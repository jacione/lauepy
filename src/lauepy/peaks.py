import json
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndi
import tifffile
from skimage.feature import peak_local_max
from skimage import draw
from progressbar import progressbar as pbar

import src.lauepy.utils as ut


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

    if config['verbose']:
        print()
        print('### Indexed substrate peaks ###')
        print(f'Min peak distance: {min_dist} pixels')
        print(f'Rel. peak threshold: {config["pkid_substrate_threshold"]}')
        print(f'Abs. peak threshold: {threshold:.3f}')

    # Find the peaks using skimage.feature.peak_local_max
    peak_coords = peak_local_max(img, min_distance=min_dist, threshold_abs=threshold, exclude_border=10)
    peak_coords = np.fliplr(peak_coords)  # This is so that they go more nicely into the following functions

    # Save the substrate peaks to the peak dictionary
    peak_dict['substrate'] = {'coords': peak_coords.tolist(), 'num_peaks': peak_coords.shape[0]}

    # Create a mask to block out the substrate peaks from the individual images
    mask_threshold = np.std(img) * config['pkid_mask_threshold']
    r = config["pkid_mask_dilation"]
    s = int(r+2)
    structure = np.zeros((2*s+1, 2*s+1))
    structure[draw.disk((s, s), r)] = 1
    sub_mask = ndi.binary_dilation(img > mask_threshold, structure=structure)

    # Add the joining lines on the detector to the mask
    sub_mask[254:258] = True
    sub_mask[:, 254:258] = True
    sub_mask[:, 512:516] = True
    sub_mask[:, 770:774] = True

    # Save the substrate mask
    np.save(f"{working_dir}/substrate/substrate_mask.npy", np.array([sub_mask]))

    end = time.perf_counter()

    if config['verbose']:
        print(f'Number of peaks found: {peak_coords.shape[0]}')
        print(f'Time to calculate: {end-start:.3f} sec')

    if config['show_plots']:
        plt.figure()
        plt.imshow(img, vmax=np.quantile(img, 0.999))
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
    threshold = np.std(img_stack, axis=(1, 2)) * config['pkid_sample_threshold']
    min_dist = config['pkid_sample_distance']

    if config['verbose']:
        print()
        print('### Indexed sample peaks ###')
        print(f'Min peak distance: {min_dist} pixels')
        print(f'Threshold rel: {config["pkid_sample_threshold"]}')
        print(f'Threshold mean: {np.mean(threshold):.1f}')
        print(f'Threshold stdv: {np.std(threshold):.1f}')

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

    # This loop saves the peaks into a dictionary where they are organized by frame
    for frame, img in enumerate(pbar(img_stack)):
        peak_coords = peak_local_max(img, min_distance=min_dist, threshold_abs=threshold[frame], exclude_border=1)
        peak_dict[f'frame_{frame:05}'] = {
            'coords': np.fliplr(peak_coords).tolist(),
            'num_peaks': peak_coords.shape[0]
        }

    end = time.perf_counter()

    # Calculate the total number of peaks and average peaks per frame.
    total_peaks = np.sum([peak_dict[frame]["num_peaks"] for frame in peak_dict.keys()])
    mean_peaks = np.mean([peak_dict[frame]["num_peaks"] for frame in peak_dict.keys()])

    if config['verbose']:
        print(f'Number of peaks found: {total_peaks}')
        print(f'Avg peaks per frame: {mean_peaks:.3f}')
        print(f'Time to calculate: {end-start:.3f} sec')

    # if config['show_plots']:
    #     print()
    #     print("Overlaying peaks...")
    #     print("This is a very slow process. If you're confident in your peakfinding parameters, you may want to"
    #           "disable the show_plots parameter.")
    #     plot_dir = Path(f"{config['working_dir']}/peaks/overlays")
    #     if plot_dir.exists():
    #         for p in plot_dir.iterdir():
    #             p.unlink()
    #     else:
    #         plot_dir.mkdir()
    #     plt.figure(tight_layout=True, figsize=(10, 5), dpi=150)
    #     for i, frame in enumerate(pbar(img_stack)):
    #         plt.cla()
    #         c = np.array(peak_dict[f'frame_{i:05}']['coords']).T
    #         if not len(c):
    #             continue
    #         plt.imshow(frame, vmax=np.quantile(frame, 0.999))
    #         plt.scatter(c[0], c[1], edgecolor='red', facecolor='None', s=160)
    #         plt.savefig(f"{plot_dir}/frame_{i:0>5}.png")

    return peak_dict


def record_positions(config, peak_dict):
    axes, spec_positions = ut.read_spec_log(config)
    peak_dict['info'] = {'angles': ut.read_spec_init(config, 'Phi', 'Chi', 'Theta').tolist()}
    peak_dict['info']['axes'] = axes  # TODO make sure that the axes labels get propagated to the grain dict.

    for frame, frame_data in peak_dict.items():
        if frame.startswith('frame_'):
            frame_data['positions'] = list(spec_positions[int(frame[-3:])])
        elif frame == 'substrate':
            frame_data['positions'] = [0 for _ in axes]

    return peak_dict
