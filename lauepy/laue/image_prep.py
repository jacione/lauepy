"""
Image cleanup to prepare for the Laue analysis
"""
import copy
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndi
import tifffile
from scipy.ndimage.measurements import label, find_objects, center_of_mass
from skimage import restoration
from progressbar import progressbar as pbar


def extract_substrate(config):
    """
    Takes a stack of Laue diffraction images and extracts only the stationary peaks. If the median background has
    the issue of uneven features, we can adopt the cleanup_images method to remove uneven pattern. The rolling-ball
    algorithm estimates the background intensity of a grayscale image in case of uneven exposure.


    :param config: configuration parameters
    :type config: dict
    """
    print('Extracting substrate-only image...')
    files = sorted(Path(config['data_dir']).iterdir())
    img_stack = np.array([tifffile.imread(f'{f}') for f in files], dtype='i')

    # The median of all frames removes short-lived features (such as Laue peaks!)
    rb_img = np.median(img_stack, axis=0)
    rb_img = remove_badpixel(rb_img)

    # Rolling ball filter removes large-scale features (such as vignette)
    rb_img = rb_img - restoration.rolling_ball(rb_img, radius=config['prep_rolling_ball_radius'])
    
    tifffile.imsave(f"{config['working_dir']}/substrate_peaks.tiff", np.array(rb_img, dtype='i'))

    return


def cleanup_images(config):
    """
    Takes the image stack and cleans them up so that the Laue peaks are easier to see.

    :param config: configuration parameters
    :type config: dict
    """
    print('Cleaning up Laue images...')
    output_dir = f"{config['working_dir']}/clean_images"
    if not Path(output_dir).exists():
        Path(output_dir).mkdir(parents=True)

    files = sorted(Path(config['data_dir']).iterdir())
    img_stack = np.array([tifffile.imread(f'{f}') for f in files], dtype='i')

    for i, img in enumerate(pbar(img_stack)):
        # Remove bad pixels
        img = remove_badpixel(img)

        # Normalize the intensity of each image
        img = np.clip(img, 0, np.quantile(img, 0.999))
        img = (img - np.quantile(img, config['prep_zero_fraction'])) / np.mean(img) * config['prep_coefficient']

        # Subtract the very broad features
        img = img - ndi.gaussian_filter(img, config['prep_gaussian_sigma'])

        # Remove negative values
        img[img < 0] = 0

        tifffile.imsave(f'{output_dir}/img_{i:05}.tiff', np.array(img, dtype='i'))

    if config['show_plots']:
        files = sorted(Path(output_dir).iterdir())
        img_stack = np.array([tifffile.imread(f'{f}') for f in files], dtype='i')

        plt.figure()
        plt.plot(np.median(img_stack, axis=(1, 2)))
        plt.title('median values of background-subtracted images')
        plt.show()

    return


def find_outlier_pixels(data, tolerance=10e7, worry_about_edges=True):
    # TODO: Unused function - find use or delete
    # This function finds the hot or dead pixels in a 2D dataset.
    # tolerance is the number of standard deviations used to cutoff the hot pixels
    # If you want to ignore the edges and greatly speed up the code, then set
    # worry_about_edges to False.
    #
    # The function returns a list of hot pixels and also an image with with hot pixels removed

    blurred = ndi.median_filter(data, size=2)
    difference = data - blurred
    threshold = tolerance

    # find the hot pixels, but ignore the edges
    hot_pixels = np.nonzero((np.abs(difference[1:-1, 1:-1]) > threshold))
    hot_pixels = np.array(hot_pixels) + 1  # because we ignored the first row and first column

    fixed_image = np.copy(data)  # This is the image with the hot pixels removed
    for y, x in zip(hot_pixels[0], hot_pixels[1]):
        fixed_image[y, x] = blurred[y, x]

    if worry_about_edges:
        height, width = np.shape(data)

        # Now get the pixels on the edges (but not the corners)

        # left and right sides
        for index in range(1, height - 1):
            # left side:
            med = np.median(data[index - 1:index + 2, 0:2])
            diff = np.abs(data[index, 0] - med)
            if diff > threshold:
                hot_pixels = np.hstack((hot_pixels, [[index], [0]]))
                fixed_image[index, 0] = med

            # right side:
            med = np.median(data[index - 1:index + 2, -2:])
            diff = np.abs(data[index, -1] - med)
            if diff > threshold:
                hot_pixels = np.hstack((hot_pixels, [[index], [width - 1]]))
                fixed_image[index, -1] = med

        # Then the top and bottom
        for index in range(1, width - 1):
            # bottom:
            med = np.median(data[0:2, index - 1:index + 2])
            diff = np.abs(data[0, index] - med)
            if diff > threshold:
                hot_pixels = np.hstack((hot_pixels, [[0], [index]]))
                fixed_image[0, index] = med

            # top:
            med = np.median(data[-2:, index - 1:index + 2])
            diff = np.abs(data[-1, index] - med)
            if diff > threshold:
                hot_pixels = np.hstack((hot_pixels, [[height - 1], [index]]))
                fixed_image[-1, index] = med

        # Then the corners

        # bottom left
        med = np.median(data[0:2, 0:2])
        diff = np.abs(data[0, 0] - med)
        if diff > threshold:
            hot_pixels = np.hstack((hot_pixels, [[0], [0]]))
            fixed_image[0, 0] = med

        # bottom right
        med = np.median(data[0:2, -2:])
        diff = np.abs(data[0, -1] - med)
        if diff > threshold:
            hot_pixels = np.hstack((hot_pixels, [[0], [width - 1]]))
            fixed_image[0, -1] = med

        # top left
        med = np.median(data[-2:, 0:2])
        diff = np.abs(data[-1, 0] - med)
        if diff > threshold:
            hot_pixels = np.hstack((hot_pixels, [[height - 1], [0]]))
            fixed_image[-1, 0] = med

        # top right
        med = np.median(data[-2:, -2:])
        diff = np.abs(data[-1, -1] - med)
        if diff > threshold:
            hot_pixels = np.hstack((hot_pixels, [[height - 1], [width - 1]]))
            fixed_image[-1, -1] = med

    return hot_pixels, fixed_image


def remove_badpixel(img):
    """

    :param img: Laue image to analyze
    :type img: ndarray
    """
    corrected = copy.copy(img)
    dead_pixel = np.where(img < 0)
    hot_pixel = np.where(img > 10e7)
    bad_pixel = (np.concatenate((dead_pixel[0], hot_pixel[0])), np.concatenate((dead_pixel[1], hot_pixel[1])))

    height, width = np.shape(img)
    for idx in range(len(bad_pixel[0])):
        y = bad_pixel[0][idx]
        x = bad_pixel[1][idx]
        if (5 < x < width - 5) and (5 < y < height - 5):
            # deal with no-edge bad pixels
            median = np.median(img[y - 5: y + 5, x - 5: x + 5].ravel())

        # Now get the pixels on the edges (but not the corners)
        elif x < 5 < y < height - 5:
            median = np.median(img[y - 5: y + 5, x: x + 5].ravel())

        # right sides
        elif x > width - 5 and 5 < y < height - 5:
            median = np.median(img[y - 5: y + 5, x - 5:x].ravel())

        # top sides
        elif width - 5 > x > 5 > y:
            median = np.median(img[y: y + 5, x - 5: x + 5].ravel())

        # bottom sides
        elif 5 < x < width - 5 and y > height - 5:
            median = np.median(img[y - 5: y, x - 5: x + 5].ravel())

        # Now get the pixels on the corners
        # top left corner
        elif x < 5 and y < 5:
            median = np.median(img[y: y + 5, x: x + 5].ravel())

        # top right corner
        elif x > width - 5 and y < 5:
            median = np.median(img[y: y + 5, x - 5: x].ravel())

        # bottom left corner
        elif x < 5 and y > height - 5:
            median = np.median(img[y - 5: y, x: x + 5].ravel())

        # bottom right corner
        elif x > width - 5 and y > height - 5:
            median = np.median(img[y - 5: y, x - 5: x].ravel())
        else:
            raise ValueError(f'median was not defined for point (x={x}, y={y})')

        # print("median", median)
        if np.isnan(corrected[y, x]):
            corrected[y, x] = np.nan_to_num(corrected[y, x])
            corrected[y, x] = median
        else:
            corrected[y, x] = median

    return corrected


def get_peaklist(image, threshold_max=1000, threshold=3, minimal_pixel=5):
    pre_processed = image
    # Blob finding
    seg_img = copy.copy(pre_processed)
    seg_img[seg_img > threshold_max] = 0

    avg = np.median(seg_img)
    sigma = np.std(seg_img)
    threshold = 3
    seg_img[seg_img < avg + threshold * sigma] = 0
    # fig = plt.figure(figsize=(10, 10), dpi=50)
    # ax = fig.add_subplot(111)
    # ax.imshow(filename[:, :],vmax=100)
    # gc.collect()
    peak = 1
    peak_dict = {}

    min_pix = minimal_pixel
    frame = seg_img[:, :]
    labels, num_feature = label(frame)  # find all clusters that could be peaks
    c = Counter(list(labels.ravel())).items()
    c = [cc[0] for cc in c if cc[1] > min_pix]

    mask = np.isin(labels, c)
    labels[~mask] = 0
    seg_img[~mask] = 0
    frame[~mask] = 0
    num_feature = len(c) - 1

    locations = find_objects(labels)  # find the slices where these peaks exist

    # c = np.arange(0,num_feature+1)

    lPos = center_of_mass(frame, labels=labels, index=c[1:])

    # lPos = [(pos[1],pos[0]) for pos in lPos]
    lpos_cor = []
    for i, pos in enumerate(lPos):
        peak_x = pos[1]
        peak_y = pos[0]
        #             if (peak_x)>256:
        #                 peak_x += 4
        #             if (peak_y)>256:
        #                 peak_y += 5

        lpos_cor.append((peak_x, peak_y))

        # lPos = np.array(lPos).reshape([-1,2])
    lPos = np.array(lpos_cor).reshape([-1, 2])

    for center in lPos:
        peak_dict['peak_%s' % peak] = {'XY': list(center), 'Center_Frame': 0}
        peak += 1

    #     end = time.time()

    #     print('Features:',num_feature)
    #     print('Number of Peaks:',len(list(peak_dict)))

    real_xys = [(peak_dict[pk]['XY']) for pk in peak_dict]
    FullPeakList = {'x0': [], 'y0': []}
    for peak in real_xys:
        FullPeakList['x0'].append(peak[0])
        FullPeakList['y0'].append(peak[1])

    # fig = plt.figure(figsize=(10, 10), dpi=50)
    # ax = fig.add_subplot(111)
    # ax.imshow(filename[:, :], vmax=100)
    # ax.scatter(FullPeakList['x0'], FullPeakList['y0'], alpha=0.7, edgecolor='red', facecolor='None', s=160)

    # open('peak_sapphire.txt', 'w') as json_file:
    #     json.dump(peak_dict, json_file)    
    # tfile.imsave(seg_img_path,np.int32(data)) 
    return peak_dict


