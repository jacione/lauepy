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

## Example

The easiest way to use LauePy (particularly if you don't care how it works) is by running
```
python lauepy/scripts/gui.py
```
That will give you this interface:

![image](https://user-images.githubusercontent.com/34489968/170314681-c63eef50-e269-441d-990d-db94823dbd72.png)

You can then adjust the parameters to your needs and execute the commands from top to bottom. If you're more comfortable with a command-line interface, there are other files in the `scripts` directory that do all of the same things, and can be run one after another.

## How to contribute

If you run into a bug, or if you have an idea for improvement, please-oh-please create an issue on the [issue tracker](https://github.com/jacione/lauepy/issues). If it's not there, it probably won't get done.

As for contributions to the code itself, we would prefer to keep the development within our research team for now. There is enough high-speed, low-level development happening that it would be difficult to manage contributions from those unfamiliar with the underlying algorithms as well as how the project fits into our overall research goals. Once the code is a bit more stable and robust, we will happily change this policy.
