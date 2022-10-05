"""
Image cleanup to prepare for the Laue analysis
"""
import copy
from collections import Counter
from pathlib import Path
import time

import matplotlib.pyplot as plt
import cupy as cp
import cupyx.scipy.ndimage as ndi
import tifffile
from scipy.ndimage.measurements import label, find_objects, center_of_mass
from progressbar import progressbar as pbar


def extract_substrate(working_dir=None, data_dir=None, prep_substrate_quantile=None, prep_substrate_sigma=None,
                      prep_substrate_radii=None, show_plots=None, **kwargs):
    """
    Takes a stack of Laue diffraction images and extracts only the stationary peaks. If the median background has
    the issue of uneven features, we can adopt the cleanup_images method to remove uneven pattern. The rolling-ball
    algorithm estimates the background intensity of a grayscale image in case of uneven exposure.

    :param working_dir: (str) Where the processed images will be saved
    :param data_dir: (str) The location of the raw images to be analyzed
    :param prep_substrate_quantile: (float) Determines the pixel-wise quantile to use when generating a substrate-only
        image.
    :param prep_substrate_sigma: (float) Width used in the Gaussian pre-filter. Generally, a value below 0.5 will
        reduce the effectiveness of the rolling-ball filter, while a value above 1.0 will reduce the visibility of
        dimmer peaks.
    :param prep_substrate_radii: (list of floats) Radii used for iterative rolling-ball background subtraction. If
        length > 1, will be applied in the order given.
    :param show_plots: (bool) If true, show a plot of before and after processing.
    :param kwargs: Any additional keyword arguments will be ignored. Required for compatibility with the main LauePy
        workflow, which is based on configuration dictionaries.
    :return: None
    """
    print('Extracting substrate-only image...')
    files = sorted(Path(data_dir).iterdir())
    img_stack = cp.array([tifffile.imread(f'{f}') for f in files], dtype='i')

    # A quantile filter over all frames removes short-lived features (such as Laue peaks!)
    # q=0.5 is a median filter
    # q=0.0 is a minimum filter
    # q=1.0 is a maximum filter
    raw_img = cp.quantile(img_stack, prep_substrate_quantile, axis=0)
    img = raw_img
    filt_img = ndi.median_filter(img, size=5)
    img[img < cp.quantile(img, 0.1)] = filt_img[img < cp.quantile(img, 0.1)]
    img[img > 10**7] = filt_img[img > 10**7]

    # Rolling ball filter removes large-scale features (such as vignette)
    sigma = prep_substrate_sigma
    kernels = make_rb_kernels(prep_substrate_radii)
    img = remove_background(img, sigma, kernels)

    if show_plots:
        plt.figure(tight_layout=True)
        ax = plt.subplot(211, xticks=[], yticks=[], title='Raw image')
        plt.imshow(raw_img.get(), vmax=cp.quantile(raw_img, 0.999))
        plt.subplot(212, xticks=[], yticks=[], title='Background removed', sharex=ax, sharey=ax)
        plt.imshow(img.get(), vmax=cp.quantile(img, 0.999))
        plt.show()
    
    tifffile.imsave(f"{working_dir}/substrate/substrate_peaks.tiff", cp.array(img, dtype='i').get())

    return


def cleanup_images(working_dir=None, data_dir=None, prep_sample_sigma=None, prep_sample_radii=None, **kwargs):
    """
    Applies a series of filters to all images found in `data_dir`.
    TODO describe the filters

    :param working_dir: (str) Where the processed images will be saved
    :param data_dir: (str) The location of the raw images to be analyzed
    :param prep_sample_sigma: (float) Width used in the Gaussian pre-filter. Generally, a value below 0.5 will
        reduce the effectiveness of the rolling-ball filter, while a value above 1.0 will reduce the visibility of
        dimmer peaks.
    :param prep_sample_radii: (list of floats) Radii used for iterative rolling-ball background subtraction. If
        length > 1, will be applied in the order given.
    :param kwargs: Any additional keyword arguments will be ignored. Required for compatibility with the main LauePy
        workflow, which is based on configuration dictionaries.
    :return:
    """
    print('Cleaning up Laue images...')
    output_dir = f"{working_dir}/clean_images"
    if not Path(output_dir).exists():
        Path(output_dir).mkdir(parents=True)

    files = sorted(Path(data_dir).iterdir())
    img_stack = cp.array([tifffile.imread(f'{f}') for f in files], dtype='i')

    t0 = time.perf_counter()

    print('Removing bad pixels...')
    filt_stack = ndi.median_filter(img_stack, size=[1, 5, 5])
    img_stack[img_stack < 0] = filt_stack[img_stack < 0]
    img_stack[img_stack > 10**7] = filt_stack[img_stack < 0]

    # Subtract the broad features
    print('Subtracting background features...')
    sigma = prep_sample_sigma
    kernels = make_rb_kernels(prep_sample_radii)
    for i, img in enumerate(pbar(img_stack)):
        img_stack[i] = remove_background(img_stack[i], sigma, kernels)

    print('Saving cleaned-up images...')
    for i, img in enumerate(pbar(img_stack)):
        tifffile.imsave(f'{output_dir}/img_{i:05}.tiff', cp.array(img, dtype='i').get())
        
    t1 = time.perf_counter()
    print(f'Total time: {t1-t0}')

    return


def make_rb_kernels(radii):
    if not isinstance(radii, list):
        radii = [radii]
    kernels = []
    for r in radii:
        kernel_size = 2*r + 1
        kernel_ctr = kernel_size // 2
        coords = cp.meshgrid(*(cp.arange(kernel_size) - kernel_ctr for _ in range(2)))
        values = cp.sum(cp.array([c ** 2 for c in coords]), axis=0)
        kernels.append(values <= r**2)
    return kernels


def remove_background(img, sigma, kernels):
    img = ndi.gaussian_filter(img, sigma)
    for k in kernels:
        img = ndi.white_tophat(img, footprint=k)
    return img


def get_peaklist(image, threshold_max=1000, threshold=3, minimal_pixel=5):
    pre_processed = image
    # Blob finding
    seg_img = copy.copy(pre_processed)
    seg_img[seg_img > threshold_max] = 0

    avg = cp.median(seg_img)
    sigma = cp.std(seg_img)
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

    mask = cp.isin(labels, c)
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
    lPos = cp.array(lpos_cor).reshape([-1, 2])

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
