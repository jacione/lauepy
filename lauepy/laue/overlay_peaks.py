import json

import matplotlib.pyplot as plt
import numpy as np
import tifffile as tiff
import yaml


def overlay(outpath, pattern, grain, threshold=300, reduction='reduction.yml', isolation='isolate_peaks.yml',
            pat_dict='pattern_dict.json', det_params='det_params.json', cryst_params='crystal_params.json', plot=True,
            frame_step=1):
    ### patter ####
    with open(reduction, 'r') as yml_file:
        x = yaml.safe_load(yml_file)
    print(grain)
    data_path = x['dataDirectory']
    img_path = x['identifier']
    with open(pat_dict) as f:
        pattern_dict = json.load(f)
        # print('Frame Goodness Pix_Dist  RMS Num_Peaks \t   Spec_Orientation \t\t Samplex_Sampley')
    patt = pattern_dict[pattern]

    ### frames #####
    #     img_path = 'LANLPUP421a_S%04d'%scan
    frame = patt['Center_Frame']
    #     print("frame is",frame)
    xys = patt['xys']

    with open(isolation, 'r') as yml_file:
        y = yaml.safe_load(yml_file)
    ind = frame - 1 // frame_step  # list(np.unique([pattern_dict[pat]['Center_Frame'] for pat in pattern_dict])).index(frame)
    start_frame = y['frames'][0]
    xs, ys = zip(*xys)
    ####### 
    p = patt
    original_path = '%s/%s/%s_%05d.tif' % (data_path, img_path, img_path, frame)
    image = tiff.imread(original_path)

    reduce_img = tiff.imread('seg_stack.tif')[ind - start_frame]

    if plot:
        fig = plt.figure(figsize=[16, 10])
        plt.subplot(1, 2, 1)
        plt.title('raw_img_scan' + '_%05d' % (p['Center_Frame']) + '%d %d %.2f %.2f %d %s %s' % (
        p['Center_Frame'], p['Goodness'], p['Dist'], p['RMS'], p['Num_Peaks'],
        np.int32(np.round(p['Spec_Orientation'], 0)).tolist(), np.round(p['Pos'], 3)))

        #         data[data<threshold] = 0
        plt.imshow(image, vmin=0, vmax=threshold, origin='upper')
        #     xtal_simu.gen_peaks_to_plot()
        red_path = '%s/%s_reduced/%s_reduced_%05d.tif' % (data_path, img_path, img_path, frame)
        plt.scatter(xs, ys, s=81, facecolors='none', edgecolors='r')
        plt.subplot(1, 2, 2)
        plt.title('red_img_scan' + '_%05d' % frame)
        plt.imshow(reduce_img, vmin=0, vmax=1, origin='upper')

        # xtal_simu.make_plots(ax)

        # print(xtal_simu.xs)
        plt.scatter(xs, ys, s=81, facecolors='none', edgecolors='r')
        plt.show()
    # plt.show()
    # print('  %d \t %d \t %.2f \t %.2f \t %d \t %s \t %s'%(p['Center_Frame'],p['Goodness'],p['Dist'],p['RMS'],p['Num_Peaks'],np.int32(np.round(p['Spec_Orientation'],0)).tolist(),np.round(p['Pos'],3)))

    tiff.imsave('%s/%s' % (outpath, grain), image.astype(np.float32))
    return
