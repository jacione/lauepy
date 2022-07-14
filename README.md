# LauePy

LauePy is a python-based code library that analyzes Laue diffraction patterns collected from a sample---either a polycrystal or separated nanocrystals---to obtain a crystal orientation map. This is a direct alternative to the more common technique of electron backscatter diffraction (EBSD). Unlike EBSD, which requires the sample to be mounted in an electron microscope, Laue mapping is performed in a coherent x-ray (CXR) beamline, making it an excellent complement to other CXR techniques such as Bragg coherent diffraction imaging (BCDI). Furthermore, by using coherent x rays instead of electrons, Laue mapping is able to penetrate deep into the sample, unlike EBSD, which only provides information about the surface orientation.

## Installation

Clone the git repository [here](https://github.com/jacione/lauepy)! Once you've downloaded it, you'll need to install it, preferably in its own environment. Navigate your terminal to the repository and run
```
pip install --upgrade pip
pip install build
python -m build
pip install .
```
Assuming I've done my job, you should be good.

## Usage

The easiest way to use LauePy is by running
```
python lauepy/scripts/gui.py
```
Each tab on the application window has a list of parameters which are linked to a configuration file. When the app is first opened, it loads an example configuration from `lauepy/config_example/config.yml`. Whenever any of the parameters is edited, the buttons to run the software are disabled until you either save or revert your changes.

### General Parameters
The `General` tab defines the experiment. Editing any of the first four parameters (and saving) will prompt LauePy to create a new working directory containing your `config.yml` file. The working directory structure was defined based on the conventions used at APS beamline 34-ID-C:
```
/home/beams/CXDUSER/{beamline}-work/{year}/{exp_id}/Analysis/lauepy_output/scan_{scan}{alt_id}
```
While this should work for any experiment done at that beamline, adjustments may be required for LauePy to work with other file systems.

| Parameter     | Config name | Description                                                                                                                                                                                              |
|---------------|-------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Scan          | scan        | Scan number in the experiment's SPEC file                                                                                                                                                                |
| Alternate ID  | alt_id      | Alternative identifier to differentiate multiple analyses of the same scan. This string will appear as a suffix to the working directory name. To leave off, enter "None".                               |
| Year          | year        | Year in which the data was taken                                                                                                                                                                         |
| Experiment ID | exp_id      | Working directory name for this beamtime (e.g. "LauePUP422")                                                                                                                                             |
| SPEC sequence | spec_seq    | Differentiates between multiple SPEC files created during a single beamtime. Appears as a suffix on the SPEC file name.                                                                                  |
| Beamline      | beamline    | Identifies the beamline where this experiment took place (e.g. "34idc")                                                                                                                                  |
| Calibration   | calibration | Detector calibration filename. This file will be loaded from the                                                                                                                                         |
| Substrate     | substrate   | Chemical formula for the substrate material. Must have a corresponding file in `lauepy/src/crystals`. If there is no crystal file for the material used in your experiment, you will need to create one. |
| Sample        | sample      | Chemical formula for the sample material (see above).                                                                                                                                                    |
| Show plots    | show_plots  | If true, produce and save plots at various stages of analysis                                                                                                                                            |
| Verbose       | verbose     | If true, print verbose output while running                                                                                                                                                              |

### Image prep tab

| Parameter          | Config name                                  | Description                                                                                                                                                                                                               |
|--------------------|----------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Quantile filter    | prep_substrate_quantile                      | Determines the pixel-wise quantile to use when generating a substrate-only image. Recommended value: 0.5                                                                                                                  |
| Gaussian sigma     | prep_substrate_sigma <br/> prep_sample_sigma | Width used in the Gaussian pre-filter. Generally, a value below 0.5 will reduce the effectiveness of the rolling ball filter, while a value above 1.0 will reduce the visibility of dimmer peaks. Recommended value: 0.75 |
| Rolling ball radii | prep_substrate_radii <br/> prep_sample_radii | Radii used for iterative rolling-ball background subtraction. Must be given in [square brackets]. Recommended value: [30, 10, 3]                                                                                          |

### Peak finding tab
| Parameter     | Config name                                        | Description                                                                                                                                                                                                          |
|---------------|----------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Threshold     | pkid_substrate_threshold<br/>pkid_sample_threshold | Minimum value of a valid peak, relative to the standard deviation intensity of the mage. A good starting value is 0.2, but the ideal value varies by experiment. Decrease to find more peaks, increase to find less. |
| Min. distance | pkid_substrate_distance<br/>pkid_sample_distance   | Minimum distance (in pixels) between valid peaks. If two otherwise valid peaks are closer together, only the brighter of the two will be returned.                                                                   |

### Laue indexing tab

| Parameter       | Config name                                        | Description                                                                         |
|-----------------|----------------------------------------------------|-------------------------------------------------------------------------------------|
| Misorientation  | laue_substrate_mis_err<br/>laue_sample_mis_err     | Maximum misorientation for two Laue diffraction patterns to be considered the same. |
| Error tolerance | laue_substrate_tolerance<br/>laue_sample_tolerance | Maximum angular error (peak-wise average) for indexed Laue patterns.                |
| Grain tolerance | grain_tolerance                                    | Maximum misorientation for two Laue diffraction patterns to be considered the same. |
| Grain threshold | grain_threshold                                    | I don't think this actually does anything...                                        |
| Twin tolerance  | twin_tolerance                                     | Maximum crystal misorientation when determining whether two grains are twins        |

If you're more comfortable with a command-line interface, there are other files in the `scripts` directory that do all the same things, and can be run one after another.

## How to contribute

If you run into a bug, or if you have an idea for improvement, please-oh-please create an issue on the [issue tracker](https://github.com/jacione/lauepy/issues). If it's not there, it probably won't get done.

As for contributions to the code itself, we would prefer to keep the development within our research team for now. There is enough high-speed, low-level development happening that it would be difficult to manage contributions from those unfamiliar with the underlying algorithms as well as how the project fits into our overall research goals. Once the code is a bit more stable and robust, we will happily change this policy.
