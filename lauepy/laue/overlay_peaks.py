import matplotlib.pyplot as plt
import numpy as np
import tifffile as tifffile


def overlay(config, pattern, grain):
    ### pattern ####
    print(grain)

    ### frames #####
    #     img_path = 'LANLPUP421a_S%04d'%scan
    frame = pattern['Center_Frame']
    #     print("frame is",frame)
    xys = pattern['xys']

    xs, ys = zip(*xys)
    ####### 
    raw_img = tifffile.imread(f"{config['data_dir']}/"
                              f"{config['exp_id']}{config['spec_seq']}_S{config['scan']:04}_{frame:05}.tif")
    clean_img = tifffile.imread(f"{config['working_dir']}/clean_images/img_{frame:05}.tiff")

    if config['show_plots']:
        plt.figure(figsize=[16, 10])
        plt.subplot(1, 2, 1)
        plt.title(f"raw_img_scan_{pattern:05}, {pattern['Goodness']}, {pattern['Dist']}, {pattern['RMS']}, "
                  f"{pattern['Num_Peaks']}, {np.int32(np.round(pattern['Spec_Orientation'], 0)).tolist()}, "
                  f"{np.round(pattern['Pos'], 3)}")
        plt.imshow(raw_img, vmin=0, vmax=config['grain_threshold'], origin='upper')
        plt.scatter(xs, ys, s=81, facecolors='none', edgecolors='r')
        plt.subplot(1, 2, 2)
        plt.title('red_img_scan' + '_%05d' % frame)
        plt.imshow(clean_img, vmin=0, vmax=1, origin='upper')
        plt.scatter(xs, ys, s=81, facecolors='none', edgecolors='r')
        plt.show()

    tifffile.imsave(f"{config['working_dir']}/grains/{grain}.tiff", raw_img.astype(np.float32))
    return
