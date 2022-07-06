import matplotlib.pyplot as plt
import numpy as np
import tifffile as tifffile


def overlay(config, pattern, grain):
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

    tifffile.imsave(f"{config['working_dir']}/grains/{grain}.tiff", raw_img.astype(np.float32))
    return
