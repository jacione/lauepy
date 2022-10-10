# LauePy

LauePy is a python-based code library that analyzes Laue diffraction patterns collected from many locations on a sample---either a polycrystal or separated nanocrystals---to obtain a crystal orientation map. The Laue diffraction data is collected on a large area detector while scanning a wide-band (or "pink") x-ray beam over a sample. As a complement to Bragg coherent diffraction imaging (BCDI), it has a distinct advantage over the more common orientation-mapping technique of electron backscatter diffraction (EBSD). Because scanning Laue analysis and BCDI use almost identical setups, one can be performed directly after the other simply by inserting or removing a monochrometor.

## Installation & setup

### APS beamline 34-ID-C
If you're using LauePy, there's a decent chance you're doing so at APS beamline 34-ID-C. In this case, all of the dependencies should already be installed on Sayre, and you just need to run the following commands for the initial setup:
```
ssh -Y cxduser@sayre
cd 34idc-work/2022/lauepy_dev
conda activate lauepy
```
Once this is done, you can (and should) skip the rest of this section.

### Requirements
LauePy makes certain assumptions about the setup, process, and data management of your experiment. These assumptions are based on APS beamline 34-ID-C, where the code was developed.

Current requirements:
* Linux operating system (for the compiled Euler binary)
* Python 3.7 (other versions probably work but haven't been tested)
* NVidia graphics card (for CUDA-accelerated calculation)
* The experiment must be managed in [SPEC](https://certif.com/content/spec/).

While a limited framework has been laid for the expansion of the software beyond this beamline, it is not currently at a plug-and-play level.

### Install
Clone the git repository [here](https://github.com/jacione/lauepy)! Once you've downloaded it, you'll need to install it, preferably in its own environment. Navigate your terminal to the repository and run
```
pip install --upgrade pip
pip install build
python -m build
pip install .
```

## Basic Usage

The easiest way to use LauePy is by running
```
python src/lauepy_scripts/gui.py
```
> Note: The path must be relative to the toplevel LauePy directory for the submodules to import properly.

Each tab on the application window has a list of parameters which are linked to a configuration file. When the app is first opened, it loads an example configuration from `lauepy/config_example/config.yml`. Changes made to the configuration are automatically saved when any command is run.

The following sections describe each tab in the application. They also include descriptions of each editable parameter, along with their corresponding name in the `config.yml` file.

### General Parameters
The "General" tab defines the experiment. Editing any of the first four parameters (and saving) will prompt LauePy to create a new working directory containing your `config.yml` file. The working directory, data, and spec file are located based on the conventions used at APS beamline 34-ID-C, as well as on the parameters provided by the user:
```
Working directory:
/home/beams/CXDUSER/{beamline}-work/{year}/{exp_id}/Analysis/lauepy_output/scan_{scan}{alt_id}
Data directory:
/home/beams/CXDUSER/{beamline}-data/{year}/{exp_id}/AD34idcLaue_{exp_id}{spec_seq}/{exp_id}{spec_seq}_S{scan:04}
SPEC file:
/home/beams/CXDUSER/{beamline}-data/{year}/{exp_id}/{exp_id}{spec_seq}.spec
```
From this tab, the user can run the full set of Laue analysis routines with a single click, or proceed to the other tabs for a more step-by-step approach.

| Parameter     | Config name | Description                                                                                                                                                                                            |
|---------------|-------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Scan          | scan        | Scan number in the experiment's SPEC file                                                                                                                                                              |
| Alternate ID  | alt_id      | Alternative identifier to differentiate multiple analyses of the same scan. This string will appear as a suffix to the working directory name. To leave off, enter "None".                             |
| Year          | year        | Year in which the data was taken                                                                                                                                                                       |
| Experiment ID | exp_id      | Working directory name for this beamtime (e.g. "LauePUP422")                                                                                                                                           |
| SPEC sequence | spec_seq    | Differentiates between multiple SPEC files created during a single beamtime. Appears as a suffix on the SPEC file name.                                                                                |
| Beamline      | beamline    | Identifies the beamline where this experiment took place (e.g. "34idc")                                                                                                                                |
| Calibration   | calibration | Detector calibration filepath. Unlike other parameters, clicking on the entry field will open a load-file dialog, which must be used to locate the calibration file.                                   |
| Substrate     | substrate   | Chemical formula for the substrate material. Must have a corresponding file in `lauepy/src/crystals`. If there is no crystal file for the material used in your experiment, you will need to make one. |
| Sample        | sample      | Chemical formula for the sample material (see above).                                                                                                                                                  |
| Show plots    | show_plots  | If true, produce and save plots at various stages of analysis                                                                                                                                          |
| Verbose       | verbose     | If true, print verbose output while running                                                                                                                                                            |


### Image prep tab

The image prep routine makes it easier to find Laue peaks. It first extracts a substrate-only image by taking a pixel-wise quantile of the image stack. This keeps only those peaks which persist through the entire image stack. All images (substrate-only included) are then passed through several filters:
1. A selective median filter replaces "dead" and "hot" pixels with more locally appropriate values. This filter is not modifiable.
2. A Gaussian filter smooths out pixel-to-pixel noise. The rolling-ball (RB) method for background subtraction (used next) is highly sensitive to this type of noise. If not removed, it will often cause the RB filter to leave much of the original background, as well as many circular artifacts.
3. A series of rolling-ball filters remove large-scale background features, isolating the Laue peaks as bright points on a dark background.

| Parameter          | Config name                                  | Description                                                                                                                                                                                                               |
|--------------------|----------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Quantile filter    | prep_substrate_quantile                      | Determines the pixel-wise quantile to use when generating a substrate-only image. Recommended value: 0.5                                                                                                                  |
| Gaussian sigma     | prep_substrate_sigma <br/> prep_sample_sigma | Width used in the Gaussian pre-filter. Generally, a value below 0.5 will reduce the effectiveness of the rolling-ball filter, while a value above 1.0 will reduce the visibility of dimmer peaks. Recommended value: 0.75 |
| Rolling ball radii | prep_substrate_radii <br/> prep_sample_radii | Radii used for iterative rolling-ball background subtraction. Must be given in [square brackets]. Recommended value: [45, 15, 5]                                                                                          |

### Peak finding tab

The peak-finding routine uses `skimage.feature.peak_local_max()` to return the coordinates of local maxima within each image. During this process, the substrate peaks are indexed first, using the substrate-only composite image. Then a mask is generated by taking a binary dilation of the set of all pixels whose values were above the substrate peak threshold. Applying this mask to the original image stack effectively removes the substrate peaks from the data, allowing the routine to focus instead on the (typically much dimmer) sample peaks. The mask also includes three columns and one row (each two pixels wide) which have been observed to respond differently from the rest of the detector.

| Parameter     | Config name                                        | Description                                                                                                                                                                                                           |
|---------------|----------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Threshold     | pkid_substrate_threshold<br/>pkid_sample_threshold | Minimum value of a valid peak, relative to the standard deviation intensity of the image. A good starting value is 0.2, but the ideal value varies by experiment. Decrease to find more peaks, increase to find less. |
| Min. distance | pkid_substrate_distance<br/>pkid_sample_distance   | Minimum distance (in pixels) between valid peaks. If two otherwise valid peaks are closer together, only the brighter of the two will be returned.                                                                    |

### Laue indexing tab

This tab deals with the actual Laue analysis. At the heart of this analysis is a compiled binary called Euler, taken from [LaueGo](https://github.com/34IDE/LaueGo), a software library developed at APS beamline 34-ID-E. Euler determines the most likely crystal orientation based on the location of Laue peaks in each frame, the lattice parameters of the diffracting material, and several tolerance factors defined in this tab.

Once the Laue diffraction patterns have been indexed, patterns that appear over multiple frames are grouped into "grains". For each grain, the frames' corresponding beam positions are averaged, weighted with the number of grain-specific peaks that appear in that frame, to determine the most likely position of the grain itself. The grain orientations are also compared to determine which pairs of grains are likely to be twin-related.

| Parameter       | Config name                                        | Description                                                                                                                               |
|-----------------|----------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------|
| Misorientation  | laue_substrate_mis_err<br/>laue_sample_mis_err     | Maximum misorientation for two Laue diffraction patterns to be considered the same.                                                       |
| Error tolerance | laue_substrate_tolerance<br/>laue_sample_tolerance | Maximum angular error (peak-wise average) for indexed Laue patterns.                                                                      |
| Grain tolerance | grain_tolerance                                    | Maximum misorientation for two Laue diffraction patterns to be considered the same.                                                       |
| Grain threshold | grain_threshold                                    | Maximum number of frames that a grain can appear in, used to filter out incorrectly indexed substrate peaks that didn't get filtered out. |
| Twin tolerance  | twin_tolerance                                     | Maximum crystal misorientation when determining whether two grains are twins                                                              |

### Output
After everything has run, the working directory (defined above) will contain the following:

```
scan_####
├── config.yml
├── clean_images
│   ├── img_00000.tiff
│   ├── ...
│   └── img_#####.tiff
├── grains
│   ├── grains.json
│   ├── grain_1.tiff
│   ├── ...
│   └── grain_#.tiff
├── macros
│   ├── grain_1.mac
│   ├── ...
│   └── grain_#.mac
├── peaks
│   ├── Index.txt
│   ├── patterns.json
│   ├── peaks.json
│   ├── Peaks.txt
│   └── overlays
│       ├── frame_00000.png
│       ├── ...
│       └── frame_#####.png
├── substrate
│   ├── rb_background.tiff
│   ├── substrate_mask.npy
│   └── substrate_peaks.tiff
└── clean_images
    ├── angle_list.txt
    ├── hiconf_twins.txt
    └── possible_twins.txt
```
Most of these are primarily for internal use, but can be useful when debugging. The most important output files are the grain macros. These essentially program a specific grain's position, orientation, and lattice parameters into SPEC. With this information, SPEC can then calculate where to position both the sample and the diffractometer arm in order to measure a specific Bragg peak.

## Alternative Usage
There are several command-line scripts in `src/lauepy_scripts` that can be run outside of the context provided by the primary LauePy workflow. These include:

| File                    | Purpose                                                          |
|-------------------------|------------------------------------------------------------------|
| image_prep_no_config.py | Run image processing routines on any stack of images             |
| index_to_macro.py       | Convert any `Index.txt` file (output by Euler) into a SPEC macro |

If you would like additional scripts to be implemented (I'm looking at you, Ross), please leave a note on my [issue tracker](https://github.com/jacione/lauepy/issues), since that's essentially my running to-do list.

## How to contribute

If you run into a bug, or if you have an idea for improvement, please-oh-please create an issue on the [issue tracker](https://github.com/jacione/lauepy/issues). If it's not there, it probably won't get done.

As for contributions to the code itself, we would prefer to keep the development within our research team for now. There is enough high-speed, low-level development happening that it would be difficult to manage contributions from those unfamiliar with the underlying algorithms as well as how the project fits into our overall research goals. Once the code is a bit more stable and robust, we will happily change this policy.
