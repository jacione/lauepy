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


def save_peaks(peak_dict, working_dir=None, **kwargs):
    with open(f"{working_dir}/peaks/peaks.json", 'w') as f:
        json.dump(peak_dict, f)


def load_peaks(working_dir=None, **kwargs):
    try:
        with open(f"{working_dir}/peaks/peaks.json", 'r') as f:
            return json.load(f)
    except FileNotFoundError:
        return {}


def find_substrate_peaks(peak_dict, working_dir=None, pkid_substrate_threshold=None, pkid_substrate_distance=None,
                         pkid_mask_threshold=None, pkid_mask_dilation=None, show_plots=True, **kwargs):
    """
    Find the coordinates of peaks in the substrate-only Laue image
    Todo update docstring
    """
    start = time.perf_counter()

    # Set up the working directory
    working_dir = working_dir

    # Load the substrate-filtered image
    img = tifffile.imread(f'{working_dir}/substrate/substrate_peaks.tiff')

    # Small median filter to remove any persistent noise
    img = ndi.median_filter(img, size=2)

    # Set the parameters for the peak-finding
    threshold = np.std(img) * pkid_substrate_threshold
    min_dist = pkid_substrate_distance

    print()
    print('### Indexed substrate peaks ###')
    print(f'Min peak distance: {min_dist} pixels')
    print(f'Rel. peak threshold: {pkid_substrate_threshold}')
    print(f'Abs. peak threshold: {threshold:.3f}')

    # Find the peaks using skimage.feature.peak_local_max
    peak_coords = peak_local_max(img, min_distance=min_dist, threshold_abs=threshold, exclude_border=10)
    peak_coords = np.fliplr(peak_coords)  # This is so that they go more nicely into the following functions

    # Save the substrate peaks to the peak dictionary
    peak_dict['substrate'] = {'coords': peak_coords.tolist(), 'num_peaks': peak_coords.shape[0]}

    # Create a mask to block out the substrate peaks from the individual images
    mask_threshold = np.std(img) * pkid_mask_threshold
    r = pkid_mask_dilation
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

    print(f'Number of peaks found: {peak_coords.shape[0]}')
    print(f'Time to calculate: {end-start:.3f} sec')

    if show_plots:
        plt.figure()
        plt.imshow(img, vmax=np.quantile(img, 0.999))
        plt.scatter(peak_coords[:, 0], peak_coords[:, 1], edgecolor='red', facecolor='None', s=160)
        plt.figure()
        plt.imshow(sub_mask)
        plt.scatter(peak_coords[:, 0], peak_coords[:, 1], edgecolor='red', facecolor='None', s=160)
        plt.show()

    return peak_dict


def find_sample_peaks(peak_dict, working_dir=None, pkid_sample_threshold=None, pkid_sample_distance=None, **kwargs):
    # todo docstring
    start = time.perf_counter()

    working_dir = working_dir

    # Load up the image files into a 3D numpy array
    files = sorted(Path(f"{working_dir}/clean_images").iterdir())
    img_stack = np.array([tifffile.imread(f'{f}') for f in files], dtype='i')

    # Set the peak-finding parameters
    threshold = np.std(img_stack, axis=(1, 2)) * pkid_sample_threshold
    min_dist = pkid_sample_distance

    print()
    print('### Indexed sample peaks ###')
    print(f'Min peak distance: {min_dist} pixels')
    print(f'Threshold rel: {pkid_sample_threshold}')
    print(f'Threshold mean: {np.mean(threshold):.1f}')
    print(f'Threshold stdv: {np.std(threshold):.1f}')

    # Load the mask from the substrate
    try:
        substrate_mask = np.load(f"{working_dir}/substrate/substrate_mask.npy")
    except FileNotFoundError:
        # This should really never happen, but if the substrate peaks haven't been indexed, then the substrate mask
        # won't exist. If that happens, then it just needs to run the substrate indexing routine before moving on.
        find_substrate_peaks(peak_dict, working_dir=working_dir, **kwargs)
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
    total_peaks = np.sum([peak_dict[frame]["num_peaks"] for frame in peak_dict.keys() if frame.startswith("frame")])
    mean_peaks = np.mean([peak_dict[frame]["num_peaks"] for frame in peak_dict.keys() if frame.startswith("frame")])

    print(f'Number of peaks found: {total_peaks}')
    print(f'Avg peaks per frame: {mean_peaks:.1f}')
    print(f'Time to calculate: {end-start:.1f} sec')

    return peak_dict


def overlay_peaks(config):
    # todo docstring
    working_dir = config["working_dir"]
    try:
        # Load up the image files into a 3D numpy array
        files = sorted(Path(f"{working_dir}/clean_images").iterdir())
        img_stack = np.array([tifffile.imread(f'{f}') for f in files], dtype='i')

        # Apply the substrate mask to the images
        substrate_mask = np.load(f"{working_dir}/substrate/substrate_mask.npy")
        substrate_mask = np.repeat(substrate_mask, img_stack.shape[0], axis=0)
        img_stack = np.ma.array(img_stack, mask=substrate_mask)

        # Load up the peak dictionary
        peak_dict = load_peaks(config)
    except FileNotFoundError:
        print("Could not overlay peaks! Make sure all the peak-finding routines have been run first.")
        return
    if config['show_plots']:
        print()
        print("Overlaying peaks. Be warned, this is a very slow process...")
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
            plt.imshow(frame, vmax=np.std(img_stack)*config["pkid_sample_threshold"])
            plt.scatter(c[0], c[1], edgecolor='red', facecolor='None', s=160)
            plt.savefig(f"{plot_dir}/frame_{i:0>5}.png")


def record_positions(peak_dict, spec_file=None, scan=None, **kwargs):
    # todo docstring
    axes, spec_positions = ut.read_spec_log(spec_file, scan)
    peak_dict['info'] = {'angles': ut.read_spec_init(spec_file, scan, 'Phi', 'Chi', 'Theta').tolist()}
    peak_dict['info']['axes'] = axes  # TODO make sure that the axes labels get propagated to the grain dict.

    for frame, frame_data in peak_dict.items():
        if frame.startswith('frame_'):
            frame_data['positions'] = list(spec_positions[int(frame[-3:])])
        elif frame == 'substrate':
            frame_data['positions'] = [0 for _ in axes]

    return peak_dict
