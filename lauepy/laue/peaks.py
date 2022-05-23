import copy
import gc
import json
import re
import time
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndi
import tifffile
from scipy.ndimage.measurements import label, center_of_mass
from scipy.spatial.distance import cdist
from skimage.feature import peak_local_max
from skimage import draw
from progressbar import progressbar as pbar

import lauepy.laue.forward_sim as fsim
from lauepy.rxlibs.xmd34 import geometry as geo
from lauepy.rxlibs.xmd34 import lattice as latt
import lauepy.laue.utils as ut


def save_peaks(config, peak_dict):
    with open(f"{config['working_dir']}/peaks/peaks.json", 'w') as f:
        json.dump(peak_dict, f)


def load_peaks(config):
    with open(f"{config['working_dir']}/peaks/peaks.json", 'r') as f:
        return json.load(f)


def find_substrate_peaks(config, peak_dict):
    """
    Find the coordinates of peaks in the substrate-only Laue image
    """
    start = time.perf_counter()

    # Set up the working directory
    working_dir = config['working_dir']

    # Load the substrate-filtered image
    img = tifffile.imread(f'{working_dir}/substrate/substrate_peaks.tiff')

    # Small median filter to remove any persistent noise
    img = ndi.median_filter(img, size=2)

    # Set the parameters for the peak-finding
    threshold = np.std(img) * config['pkid_substrate_threshold']
    min_dist = config['pkid_substrate_distance']

    # Find the peaks using skimage.feature.peak_local_max
    peak_coords = peak_local_max(img, min_distance=min_dist, threshold_abs=threshold, exclude_border=10)
    peak_coords = np.fliplr(peak_coords)  # This is so that they go more nicely into the following functions

    # Save the substrate peaks to the peak dictionary
    peak_dict['substrate'] = {'coords': peak_coords.tolist(), 'num_peaks': peak_coords.shape[0]}

    # Create a mask to block out the substrate peaks from the individual images
    structure = np.zeros((15, 15))
    structure[draw.disk((7, 7), 5.5)] = 1
    sub_mask = ndi.binary_dilation(img > threshold, structure=structure)
    sub_mask[254:258] = True
    sub_mask[:, 254:258] = True
    sub_mask[:, 512:516] = True
    sub_mask[:, 770:774] = True
    np.save(f"{working_dir}/substrate/substrate_mask.npy", np.array([sub_mask]))

    end = time.perf_counter()

    if config['verbose']:
        print()
        print('### Indexed substrate peaks ###')
        print(f'Min peak distance: {min_dist} pixels')
        print(f'Rel. peak threshold: {config["pkid_substrate_threshold"]}')
        print(f'Abs. peak threshold: {threshold:0.3}')
        print(f'Number of peaks found: {peak_coords.shape[0]}')
        print(f'Time to calculate: {end-start: 0.3} sec')

    if config['show_plots']:
        plt.figure()
        plt.imshow(img, vmax=np.quantile(img, 0.95))
        plt.scatter(peak_coords[:, 0], peak_coords[:, 1], edgecolor='red', facecolor='None', s=160)
        plt.figure()
        plt.imshow(sub_mask)
        plt.scatter(peak_coords[:, 0], peak_coords[:, 1], edgecolor='red', facecolor='None', s=160)
        plt.show()

    return peak_dict


def find_sample_peaks(config, peak_dict):
    start = time.perf_counter()

    working_dir = config['working_dir']

    # Load up the image files into a 3D numpy array
    files = sorted(Path(f"{working_dir}/clean_images").iterdir())
    img_stack = np.array([tifffile.imread(f'{f}') for f in files], dtype='i')

    # Set the peak-finding parameters
    threshold = np.std(img_stack) * config['pkid_sample_threshold']
    min_dist = config['pkid_sample_distance']

    # Load the mask from the substrate
    try:
        substrate_mask = np.load(f"{working_dir}/substrate/substrate_mask.npy")
    except FileNotFoundError:
        # This should really never happen, but if the substrate peaks haven't been indexed, then the substrate mask
        # won't exist. If that happens, then it just needs to run the substrate indexing routine before moving on.
        find_substrate_peaks(config, peak_dict)
        substrate_mask = np.load(f"{working_dir}/substrate/substrate_mask.npy")
    # Extend the substrate mask to match the shape of the image stack
    substrate_mask = np.repeat(substrate_mask, img_stack.shape[0], axis=0)

    # Apply the substrate mask onto the image stack
    img_stack = np.ma.array(img_stack, mask=substrate_mask)

    # This loop does two things: (1) save the peaks into a dictionary where they are organized by frame, and (2) find
    # the frame with the highest number of peaks, so that it can be plt.imshown later on.
    for frame, img in enumerate(img_stack):
        peak_coords = peak_local_max(img, min_distance=min_dist, threshold_abs=threshold, exclude_border=10)
        peak_dict[f'frame_{frame:05}'] = {
            'coords': np.fliplr(peak_coords).tolist(),
            'num_peaks': peak_coords.shape[0]
        }

    end = time.perf_counter()

    # Calculate the total number of peaks and average peaks per frame.
    total_peaks = np.sum([peak_dict[frame]["num_peaks"] for frame in peak_dict.keys()])
    mean_peaks = np.mean([peak_dict[frame]["num_peaks"] for frame in peak_dict.keys()])

    if config['verbose']:
        print()
        print('### Indexed sample peaks ###')
        print(f'Min peak distance: {min_dist} pixels')
        print(f'Rel. peak threshold: {config["pkid_sample_threshold"]}')
        print(f'Abs. peak threshold: {threshold:0.3}')
        print(f'Number of peaks found: {total_peaks}')
        print(f'Avg peaks per frame: {mean_peaks}')
        print(f'Time to calculate: {end-start: 0.3} sec')

    if config['show_plots']:
        print("Overlaying peaks...")
        plot_dir = Path(f"{config['working_dir']}/peaks/overlays")
        if plot_dir.exists():
            for p in plot_dir.iterdir():
                p.unlink()
        else:
            plot_dir.mkdir()
        plt.figure(tight_layout=True, figsize=(10, 5), dpi=150)
        for i, frame in enumerate(pbar(img_stack)):
            plt.cla()
            c = np.array(peak_dict[f'frame_{i:05}']['coords']).T
            if not len(c):
                continue
            plt.imshow(frame, vmax=config['prep_coefficient']*0.25)
            plt.scatter(c[0], c[1], edgecolor='red', facecolor='None', s=160)
            plt.savefig(f"{plot_dir}/frame_{i:0>5}.png")

    return peak_dict


def record_positions(config, peak_dict):
    spec_positions = ut.read_spec_log(config)

    for frame, frame_data in peak_dict.items():
        if frame.startswith('frame_'):
            frame_data['positions'] = list(spec_positions[int(frame[-3:])])
        elif frame == 'substrate':
            frame_data['positions'] = [0, 0, 0, 0, 0, 0]

    peak_dict['info'] = {'angles': ut.read_spec_init(config, 'Phi', 'Chi', 'Theta').tolist()}
    return peak_dict


# Below this line lives the old code. It's here in case we need it.


def get_peak_coords(img, min_dist, threshold, border=10):
    return peak_local_max(img, min_distance=min_dist, threshold_abs=threshold, exclude_border=border)


def isolate_peaks(input_yml, substrate_sim_peak_dict, distance=5):
    """
    Old function for non-substrate peak finding.
    """
    print("isolate peaks")
    with open(input_yml, 'r') as yml_file:
        x = json.load(yml_file)
    peak_dict_path = x['peak_dict_path']
    seg_img_path = x['seg_img_path']
    verbose = x['verbose']
    frames = x['frames']
    threshold = x['threshold']
    data_path = x['data_path']
    dir_path = x['dir_path']
    threshold_max = x['threshold_max']
    min_pix = x['min_pix']
    max_pix = x['max_pix']
    start = time.perf_counter()

    with open(substrate_sim_peak_dict, 'r') as yml_file:
        y = json.load(yml_file)
    substrate = y['substrate_peak']

    data = []
    center_frames = []

    new_bkg = tifffile.imread(data_path + dir_path + '/' + dir_path + '_bkg_rb')
    bright_bkg = new_bkg > 4 * np.median(new_bkg.ravel())
    bright_bkg = ndi.binary_dilation(bright_bkg, iterations=3)

    for i in range(frames[0], frames[1], frames[2]):
        center_frames.append(i)

        raw_img = tifffile.imread(data_path + dir_path + '/' + dir_path + '_%05d.tif' % i).astype(float)
        # fig = plt.figure(figsize=(30,30))
        # plt.title('the raw_img image is',fontsize = 12)
        # plt.imshow(raw_img, vmin=0, vmax=100)
        # segment data
        # raw = np.int32(raw_img)
        # fig = plt.figure(figsize=(40,40))
        # plt.title(f'raw image {idx}',fontsize = 12)
        # plt.imshow(raw,vmax = 200)
        # corrected_raw = remove_badpixel(raw)
        # rb_bkg = median_filter(corrected_raw, size = 50) 
        # print('median_filter mean',np.mean(rb_bkg))
        # print('median_filter standard deviation',np.std(rb_bkg))
        # fig = plt.figure(figsize=(40, 40))
        # plt.title('rolling ball image')
        # plt.imshow(rb_bkg, vmin=0, vmax=100)
        # print('corrected_raw mean',np.mean(corrected_raw))
        # print('corrected_raw standard deviation',np.std(corrected_raw))
        # corrected_img = corrected_raw - rb_bkg
        # print('corrected_img mean ',np.mean(corrected_img))
        # print('corrected_img standard deviation',np.std(corrected_img))
        seg_img = copy.copy(raw_img)

        seg_img[seg_img > threshold_max] = 0

        avg = np.median(seg_img)
        sigma = np.std(seg_img)
        # print('avg',avg)
        # print('sigma',sigma)
        seg_img[seg_img < threshold] = 0  # large threshold will give less peaks
        #         gaussian_filter_seg = filters.gaussian(seg_img, sigma=2)
        #         med_filter_seg = filters.median(seg_img, np.ones((3, 3)))
        #         tv_filter_seg = restoration.denoise_tv_chambolle(seg_img, weight=0.1)

        #         plt.figure(figsize=(16, 16))
        #         plt.subplot(221)
        #         plt.imshow(seg_img, vmin=0, vmax=100)
        #         plt.axis('off')
        #         plt.title('Image')
        #         plt.subplot(222)
        #         plt.imshow(gaussian_filter_seg, vmin=0, vmax=100)
        #         plt.axis('off')
        #         plt.title('Gaussian filter')
        #         plt.subplot(223)
        #         plt.imshow(med_filter_seg, vmin=0, vmax=100)
        #         plt.axis('off')
        #         plt.title('Median filter')
        #         plt.subplot(224)
        #         plt.imshow(tv_filter_seg, vmin=0, vmax=100)
        #         plt.axis('off')
        #         plt.title('TV filter')
        #         plt.show()
        # fig = plt.figure(figsize=(30,30))
        # plt.title('the raw seg_img image is',fontsize = 12)
        # plt.imshow(seg_img, vmin=0, vmax=100) 
        # seg_img[seg_img<threshold] = 0
        # read in segmented data

        data.append(seg_img)
    data = np.array(data)
    #     data = np.moveaxis(np.array(data),0,-1)

    gc.collect()
    peak = 1
    peak_dict = {}

    for z, center_frame in enumerate(center_frames):
        print(f'frame {center_frame}')

        frame = data[z, :, :]
        print('max value is ', max(frame.ravel()))
        print('min value is ', min(frame.ravel()))

        labels, num_feature = label(frame)  # find all clusters that could be peaks
        c = Counter(list(labels.ravel())).items()
        c = [cc[0] for cc in c if min_pix < cc[1] < max_pix]
        # abnormal_c = [cc[0] for cc in c if cc[1]<= min_pix or cc[1] >= max_pix]
        # c = [cc[0] for cc in c if cc[1]<max_pix]
        # print(abnormal_c)
        ### the following procedure get the list of peaks ###
        mask = np.isin(labels, c)
        labels[~mask] = 0
        data[z][~mask] = 0
        frame[~mask] = 0
        num_feature = len(c) - 1

        #         c = np.arange(0,num_feature+1)

        lPos = center_of_mass(frame, labels=labels, index=c[1:])
        # print(lPos)
        # lPos = [(pos[1],pos[0]) for pos in lPos]
        lpos_cor = []
        for i, pos in enumerate(lPos):

            peak_x = pos[1]
            peak_y = pos[0]
            #             if (peak_x)>256:
            #                 peak_x += 4
            #             if (peak_y)>256:
            #                 peak_y += 5
            # print('peak_x',peak_x)
            # print('peak_y',peak_y)
            # print(np.shape(bright_bkg))
            if not bright_bkg[int(peak_y), int(peak_x)]:
                lpos_cor.append((peak_x, peak_y))

            # lPos = np.array(lPos).reshape([-1,2])
        if len(lpos_cor):
            lPos = delete_substrate(peak_dict=lPos, substrate_peak=substrate, d=distance)
        lPos = np.array(lpos_cor).reshape([-1, 2])

        FullPeakList = {'x0': [], 'y0': []}
        for center in lPos:
            # print('peak',peak, 'XY', center, 'Center_Frame',center_frame)
            peak_dict['peak_%s' % peak] = {'XY': list(center), 'Center_Frame': center_frame}
            peak += 1
            # peak_dict = delete_substrate(peak_dict = peak_dict, substrate_peak = substrate, startID = start_idx)
            # print(peak_dict)

            # real_xys_frame = [(peak_dict[pk]['XY']) for pk in peak_dict]
            # FullPeakList = {}
            # FullPeakList['x0'] = []
            # FullPeakList['y0'] = []
            # for peak in real_xys:
            FullPeakList['x0'].append(center[0])
            FullPeakList['y0'].append(center[1])
        fig = plt.figure(figsize=(25, 25), dpi=50)
        # plt.title(f"threshold value is {threshold}")
        ax = fig.add_subplot(111)
        # ax.title.set_text(f'threshold value is {threshold}')
        # ax.imshow(filename[:, :],vmax=5e2)
        ax.imshow(frame[:, :], vmax=100)
        ax.scatter(FullPeakList['x0'], FullPeakList['y0'], alpha=0.7, edgecolor='yellow', facecolor='None', s=160)
        # tfile.imsave(data_path+dir_path+'/'+dir_path+'_%05d_isolated_peaks.tif'%i,np.int32(fig)) 
        fig.savefig(data_path + dir_path + '/' + dir_path + '_%05d_isolated_peaks.tif' % z)

    end = time.perf_counter()
    if verbose:
        print('Features:', num_feature)
        print('Number of Peaks:', len(list(peak_dict)))
        print('time to calculate:', end - start, 's')
    with open(peak_dict_path, 'w') as json_file:
        json.dump(peak_dict, json_file)
    tifffile.imsave(seg_img_path, np.int32(data))
    return


def extract_spec(config, frame_nums):
    with open(config['spec_file']) as f:
        file_name = f.read()
    scan = config['scan']
    find_scan = re.compile(r"#S %s (.*?)#S %s" % (scan, scan + 1), re.DOTALL)
    find_ang = re.compile(r"#P0 (.*?)\n", re.DOTALL)
    # find_samp = re.compile(r"#P2(.*?)\n",re.DOTALL)
    find_piezo = re.compile(r"(-?\d+.*?)\n")
    scan_txt = find_scan.findall(file_name)
    if len(scan_txt) == 0:
        find_scan = re.compile(r"#S %s .*" % scan, re.DOTALL)
        scan_txt = find_scan.findall(file_name)

    vals = find_ang.findall(scan_txt[0])[0].split(" ")

    theta, chi, phi = np.array(vals[1:4], dtype=float)
    print('theta', theta, 'chi', chi, 'phi', phi)
    new = scan_txt[0].split('Monitor  Detector\n')[1]
    new1 = new.split('#')[0]
    # sampvals = find_samp.findall(scan_txt[0])[0].split(" ")
    # samp = np.array(sampvals[3:5],dtype=float)

    extracted = []

    frames = find_piezo.findall(new1)
    # print(scan_txt)

    for num in frame_nums:
        p = np.array(frames[num].split(" ")[:2], dtype=float)

        extracted.append([list(p), list([phi, chi, theta])])
    return extracted


def delete_substrate(peak_dict, substrate_peak, d=5):
    xys_sample = peak_dict
    xys_substrate = substrate_peak
    # for peak in peak_dict:
    #     x = peak_dict[peak]['XY'][0]
    #     y = peak_dict[peak]['XY'][1]
    #     xys_sample.append((x,y))
    # print(peak,peak_dict[peak])
    peak_list_sample = np.array(xys_sample)
    peak_list_substrate = np.array(xys_substrate)
    # print('peak_list_sample',peak_list_sample, peak_list_sample.shape)
    # print('peak_list_substrate',peak_list_substrate,peak_list_substrate.shape)
    d = cdist(peak_list_sample, peak_list_substrate, metric='euclidean')
    close_points = np.where(d < d)
    # print('close_points[0]',close_points[0], type(close_points[0]))
    # print('set(close_points[0].tolist())',set(close_points[0].tolist()))
    # print('list(set(close_points[0].tolist()))',list(set(close_points[0].tolist())),
    # type(list(set(close_points[0].tolist()))))
    temp = sorted(list(set(close_points[0].tolist())))
    # print('')
    remove_idx = reversed(temp)
    # print('remove_idx',remove_idx)
    for idx in remove_idx:
        # print('idx', idx)
        # print('point', xys_sample[idx])
        del xys_sample[idx]
        # del peak_dict['peak_'+str(idx + startID)]

    return xys_sample


def simulate_substrate(pattern_dict, crystal_path, det_path, raw_image):
    with open(crystal_path) as f:
        crystal_params = json.load(f)
    pos_list = crystal_params['pos_list']

    a, b, c, ang1, ang2, ang3 = crystal_params['lattice_params']
    atom_list = [latt.AtomInCell(atom_pos[0], *atom_pos[1]) for atom_pos in pos_list]
    xtal = latt.Xtal(a, b, c, ang1, ang2, ang3, atomlist=atom_list)

    with open(pattern_dict, 'r') as yml_file:
        x = json.load(yml_file)
    xtal_rmat = np.array(x['pattern_1']['Rot_mat'])

    with open(det_path) as f:
        det_params = json.load(f)
    # trans_vec = det_params['trans_vec']
    # pix_x,pix_y,pitch_x,pitch_y,name = det_params['pixels']
    det_name = det_params['det_name']

    rot_vec = np.array(det_params['rot_vec'])
    trans_vec = np.array(det_params['trans_vec'])
    args = det_params['pixels']
    ccd1 = geo.Detector(*args)
    # geo_ccd = geo.DetectorGeometry('ccd1',ccd1, trans_vec, rot_vec)
    # det_params = det_params['pixels']
    # detector = geo.Detector(pix_x,pix_y,pitch_x,pitch_y,det_type)
    geo_det = geo.DetectorGeometry(det_name, ccd1, trans_vec, rot_vec)
    # # normalize vectors in rot matrix
    for i in (0, 1, 2):
        vec = xtal_rmat[i]
        vec = vec / np.linalg.norm(vec)
        xtal_rmat[i] = vec

    energy_cutoff_keV = [3.0, 20]
    xtal_simu = fsim.LaueSimulation(geo_det, xtal, xtal_rmat.T, keV_cutoffs=energy_cutoff_keV,
                                  plot_range_padding_percent=0)

    fig = plt.figure(figsize=(20, 20), dpi=50)
    ax = fig.add_subplot(111)
    # if isinstance(raw_image, str):
    #     image = tfile.imread(raw_image)
    # else:
    #     image = raw_image
    # ax.imshow(image, vmin=0, vmax=100, interpolation='nearest')
    # ax.scatter(FullPeakList['x0'], FullPeakList['y0'], alpha=0.7, edgecolor='red', facecolor='None', s=160)

    xtal_simu.gen_peaks_to_plot()
    xtal_simu.make_plots(ax)
    print('simulate')

    fig.savefig('simulated_sapphire')
    # plt.scatter(FullPeakList['x0'], FullPeakList['y0'], alpha=0.5, edgecolor='yellow', facecolor='yellow', s = 81)
    plt.show()
    xys_substrate = list(zip(xtal_simu.xs, xtal_simu.ys))

    substrate_sim_dict = {"substrate_peak": xys_substrate}

    with open('substrate_sim_peak_dict.json', 'w') as yml_file:
        json.dump(substrate_sim_dict, yml_file)
    # substrate = y['substrate_peak']

    # crystal_params = '../../Crystal_Parameters/Au_crystal_params.json'
    # crystal = {"pos_list": pos_list,
    #               "material": 'Au',
    #               "lattice_params": latts,
    #               "space_group": 225}

    # with open(crystal_params,'w') as file:
    #     json.dump(crystal,file)

    return xys_substrate
