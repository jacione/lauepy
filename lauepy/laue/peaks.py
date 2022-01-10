import copy
import gc
import json
import re
import time
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.ndimage as ndi
import tifffile
from scipy.ndimage.measurements import label, find_objects, center_of_mass
from scipy.spatial.distance import cdist

import lauepy.laue.forward_sim as fsim
from lauepy.rxlibs.xmd34 import geometry as geo
from lauepy.rxlibs.xmd34 import lattice as latt


def isolate_substrate_peaks(config):
    """ This function takes a numpy array of data and finds the center of mass of all peaks,
     thresholding by a minimum peak width in pixels

    inputs:
            data_path: path to reduced data
            weight_path: path to weights

    outputs:
            center_xy: list of the xy values of the centers of mass for all peaks
            center_z: list of the z values (frames) of the centers of mass of each peak
            z_len: the length in z that a peak exists (i.e. number of frames)

    """
    output_dir = f"{config['working_dir']}/{config['working_id']}"
    threshold = config['pkid_threshold']
    cutoff = config['pkid_cutoff']
    min_pix = config['pkid_min']
    max_pix = config['pkid_max']
    start = time.perf_counter()

    raw_img = tifffile.imread(f'{output_dir}/substrate_peaks.tiff')
    # segment data
    seg_img = np.copy(raw_img)
    seg_img[seg_img > cutoff] = 0  # FIXME does this do anything if hot pixels have already been removed?

    avg = np.median(seg_img)
    sigma = np.std(seg_img)
    seg_img[seg_img < avg + threshold * sigma] = 0  # FIXME should this be (avg + threshold) * sigma?

    # seg_img[seg_img<threshold] = 0
    # read in segmented data

    data = np.array(seg_img)
    #     data = np.moveaxis(np.array(data),0,-1)

    gc.collect()  # FIXME vestigial?
    substrate_peak_dict = {}

    labels, _ = label(data)  # find all clusters that could be peaks
    c = Counter(list(labels.ravel())).items()
    c = [cc[0] for cc in c if min_pix < cc[1] < max_pix]
    # abnormal_c = [cc[0] for cc in c if cc[1]<= min_pix or cc[1] >= max_pix]
    # c = [cc[0] for cc in c if cc[1]<max_pix]
    # print(abnormal_c)
    ### the following procedure get the list of peaks ###
    mask = np.isin(labels, c)
    labels[~mask] = 0
    data[~mask] = 0
    data[~mask] = 0
    num_feature = len(c) - 1

    locations = find_objects(labels)  # find the slices where these peaks exist

    #         c = np.arange(0,num_feature+1)

    lPos = center_of_mass(data, labels=labels, index=c[1:])
    # print(lPos)
    # lPos = [(pos[1],pos[0]) for pos in lPos]
    lpos_cor = []
    for i, pos in enumerate(lPos):
        peak_x = pos[1]
        peak_y = pos[0]
        # if (peak_x)>256:
        #     peak_x += 4
        # if (peak_y)>256:
        #     peak_y += 5

        lpos_cor.append((peak_x, peak_y))

        # lPos = np.array(lPos).reshape([-1,2])
    lPos = np.array(lpos_cor).reshape([-1, 2])

    for i, center in enumerate(lPos):
        substrate_peak_dict[f'peak_{i}'] = {'XY': list(center)}
    # substrate_peak_dict = delete_substrate(peak_dict = peak_dict, substrate_peak = substrate, startID = start_idx)

    end = time.perf_counter()

    real_xys = [(substrate_peak_dict[pk]['XY']) for pk in substrate_peak_dict]
    FullPeakList = {'x0': [], 'y0': []}
    for peak in real_xys:
        FullPeakList['x0'].append(peak[0])
        FullPeakList['y0'].append(peak[1])

    if config['verbose']:
        print('Features:', num_feature)
        print('Number of Peaks:', len(list(substrate_peak_dict)))
        print('time to calculate:', end - start, 's')
    with open(f'{output_dir}/substrate_peaks.json', 'w') as json_file:
        json.dump(substrate_peak_dict, json_file)
    tifffile.imsave(f'{output_dir}/segmented_data.tiff', np.int32(data))

    if config['show_plots']:
        fig = plt.figure(figsize=(10, 10), dpi=50)
        ax = fig.add_subplot(111)
        ax.imshow(raw_img[:, :], vmax=100)
        ax.scatter(FullPeakList['x0'], FullPeakList['y0'], alpha=0.7, edgecolor='red', facecolor='None', s=160)

    return


def isolate_peaks(input_yml, substrate_sim_peak_dict, distance=5):
    """ This function takes a numpy array of data and finds the center of mass of all peaks,
     thresholding by a minimum peak width in pixels

    inputs:
            data_path: path to reduced data
            weight_path: path to weights

    outputs:
            center_xy: list of the xy values of the centers of mass for all peaks
            center_z: list of the z values (frames) of the centers of mass of each peak
            z_len: the length in z that a peak exists (i.e. number of frames)

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


def group_peaks(input_yml, min_peaks=3, max_peaks=250):
    """ 
    This function groups peaks by the z values of their centers of mass. Peaks that 
    are centered on the same frame will be in roughly the same location. After this,
    they are sorted by their length in z because larger length means larger grain. The 
    output should be groups of grains which are roughly the same size in the same
    position on the grid.

    inputs:
            center_xy: list of xy values for all peaks
            center_z: list of z values for all peaks
            min_peaks: minimum number of peaks to be considered a group
    outputs:
            df: pandas dataframe of the groups
            XYS: List of xys for each group
    """
    print('group peaks')
    with open(input_yml, 'r') as yml_file:
        x = json.load(yml_file)
    scan = x['scan']
    peak_dict_path = x['peak_dict_path']
    group_dict_path = x['group_dict_path']
    spec_file = x['spec_file']
    plot = x['plot']

    with open(peak_dict_path) as f:
        peak_dict = json.load(f)
    center_z = [peak_dict[pk]['Center_Frame'] for pk in peak_dict]
    group_dict = {}
    d = {'Peak_ID': range(1, len(center_z) + 1), 'Center_Z': center_z}
    df = pd.DataFrame(data=d)

    #     df = df.set_index('Peak_ID')
    groups = df.groupby('Center_Z')['Peak_ID'].apply(list).reset_index(name='Peak_ID')  # group peaks by center z value
    #     groups1 = df.groupby('Center_Z')['Peak_ID'].apply(list).reset_index(name='ID')

    IDS = groups['Peak_ID'].to_list()
    df = df.groupby('Center_Z').count().reset_index()

    idx_low = df.index[df['Peak_ID'] > min_peaks].to_list()
    idx_high = df.index[df['Peak_ID'] < max_peaks].to_list()
    idx = [value for value in idx_low if value in idx_high]

    #     print(idx)
    df = df.loc[df['Peak_ID'] > min_peaks]
    df = df.loc[df['Peak_ID'] < max_peaks]
    idx = df.index.tolist()
    #     print(idx)
    IDS = [IDS[i] for i in idx]

    zs = groups['Center_Z'].to_list()
    zs = [zs[i] for i in idx]
    n = len(IDS)
    spec_vals = extract_spec(spec_file, scan, zs)
    for i in range(n):
        group_dict['group_%s' % (i + 1)] = {'Center_Frame': zs[i], 'ID_List': IDS[i], 'phichitheta': spec_vals[i][1],
                                            'Pos': spec_vals[i][0]}
    print('Number of Groups:', len([grp for grp in group_dict]))
    if plot:
        IDS = [len(id) for id in IDS]
        plt.scatter(zs, IDS)
        plt.ylabel('number of peaks')
        plt.xlabel('frame')
        plt.show()

    with open(group_dict_path, 'w') as json_file:
        json.dump(group_dict, json_file)
    return


def group_substrate_peaks(config, min_peaks=3, max_peaks=250):
    """ 
    This function groups peaks by the z values of their centers of mass. Peaks that 
    are centered on the same frame will be in roughly the same location. After this,
    they are sorted by their length in z because larger length means larger grain. The 
    output should be groups of grains which are roughly the same size in the same
    position on the grid.

    inputs:
            center_xy: list of xy values for all peaks
            center_z: list of z values for all peaks
            min_peaks: minimum number of peaks to be considered a group
    outputs:
            df: pandas dataframe of the groups
            XYS: List of xys for each group
    """
    output_dir = f"{config['working_dir']}/{config['working_id']}"
    scan = config['scan']
    spec_file = config['spec_file']

    with open(f"{output_dir}/substrate_peaks.json") as f:
        peak_dict = json.load(f)
    center_z = [peak_dict[pk]['Center_Frame'] for pk in peak_dict]
    group_dict = {}
    d = {'Peak_ID': range(1, len(center_z) + 1), 'Center_Z': center_z}
    df = pd.DataFrame(data=d)

    #     df = df.set_index('Peak_ID')
    groups = df.groupby('Center_Z')['Peak_ID'].apply(list).reset_index(name='Peak_ID')  # group peaks by center z value
    #     groups1 = df.groupby('Center_Z')['Peak_ID'].apply(list).reset_index(name='ID')

    IDS = groups['Peak_ID'].to_list()
    df = df.groupby('Center_Z').count().reset_index()

    idx_low = df.index[df['Peak_ID'] > min_peaks].to_list()
    idx_high = df.index[df['Peak_ID'] < max_peaks].to_list()
    idx = [value for value in idx_low if value in idx_high]

    #     print(idx)
    df = df.loc[df['Peak_ID'] > min_peaks]
    df = df.loc[df['Peak_ID'] < max_peaks]
    idx = df.index.tolist()
    #     print(idx)
    IDS = [IDS[i] for i in idx]

    zs = groups['Center_Z'].to_list()
    zs = [zs[i] for i in idx]
    n = len(IDS)
    spec_vals = extract_spec(spec_file, scan, zs)
    for i in range(n):
        group_dict['group_%s' % (i + 1)] = {'Center_Frame': zs[i], 'ID_List': IDS[i], 'phichitheta': spec_vals[i][1],
                                            'Pos': spec_vals[i][0]}
    print('Number of Groups:', len([grp for grp in group_dict]))
    if config['show_plots']:
        IDS = [len(i) for i in IDS]
        plt.scatter(zs, IDS)
        plt.ylabel('number of peaks')
        plt.xlabel('frame')
        plt.show()

    with open(f"{output_dir}/substrate_groups.json", 'w') as json_file:
        json.dump(group_dict, json_file)
    return


def extract_spec(file, scan, frame_nums):
    file_name = open(file).read()
    find_scan = re.compile(r"#S %s (.*?)#S %s" % (scan, scan + 1), re.DOTALL)
    find_ang = re.compile(r"#P0 (.*?)\n", re.DOTALL)
    #     find_samp = re.compile(r"#P2(.*?)\n",re.DOTALL)
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
    #     sampvals = find_samp.findall(scan_txt[0])[0].split(" ")
    #     samp = np.array(sampvals[3:5],dtype=float)

    extracted = []

    frames = find_piezo.findall(new1)
    #     print(scan_txt)

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
