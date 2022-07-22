import tkinter as tk
from tkinter import ttk, filedialog
from pathlib import Path

import yaml

from src.lauepy import utils as ut
from src.lauepy import image_prep as prep
from src.lauepy import peaks as pk
import src.lauepy.auto_laue as al
import src.lauepy.make_grain_dict as grain
import src.lauepy.find_twins as twins


PARAMS = {
    "scan": "Scan number",
    "alt_id": "Alternate ID",
    "year": "Year",
    "exp_id": "Experiment ID",
    "spec_seq": "SPEC sequence",
    "beamline": "Beamline",
    "calibration": "Calibration",
    "substrate": "Substrate",
    "sample": "Sample",
    "show_plots": "Show plots",
    "verbose": "Verbose",
    "prep_substrate_quantile": "Quantile filter",
    "prep_substrate_sigma": "Gaussian filter sigma",
    "prep_substrate_radii": "Rolling ball radii",
    "prep_sample_sigma": "Gaussian filter sigma",
    "prep_sample_radii": "Rolling ball radii",
    "pkid_substrate_threshold": "Peak threshold",
    "pkid_substrate_distance": "Min. distance",
    "pkid_mask_threshold": "Mask threshold",
    "pkid_mask_dilation": "Mask dilation",
    "pkid_sample_threshold": "Peak threshold",
    "pkid_sample_distance": "Min. distance",
    # "laue_substrate_goodness": "Goodness",
    "laue_substrate_mis_err": "Misorientation",
    "laue_substrate_tolerance": "Error tolerance",
    # "laue_substrate_frequency": "Frequency",
    # "laue_substrate_comb_sub": "Comb sub",
    # "laue_substrate_times": "Times",
    # "laue_sample_goodness": "Goodness",
    "laue_sample_mis_err": "Misorientation",
    "laue_sample_tolerance": "Error tolerance",
    # "laue_sample_frequency": "Frequency",
    # "laue_sample_comb_sub": "Comb sub",
    # "laue_sample_times": "Times",
    "grain_tolerance": "Grain tolerance",
    "grain_threshold": "Grain threshold",
    "twin_tolerance": "Twin tolerance"
}

TIPS = {
    "scan": "Scan number",
    "alt_id": "Alternative identifier to differentiate multiple analyses of the same scan. To ignore, enter \"None\".",
    "year": "Year",
    "exp_id": "Working directory for this beamtime (e.g. \"LauePUP422\")",
    "spec_seq": "Differentiates between multiple SPEC files created during a single beamtime. Appears as a suffix on "
                "the SPEC file name.",
    "beamline": "Identifies the beamline where this experiment took place (e.g. \"34idc\")",
    "calibration": "Detector calibration file. Must include *.json extention",
    "substrate": "Chemical formula for the substrate material. Must have a corresponding file in lauepy/src/crystals.",
    "sample": "Chemical formula for the sample material. Must have a corresponding file in lauepy/src/crystals.",
    "show_plots": "Show plots",
    "verbose": "Verbose",
    "prep_substrate_quantile": "Determines the pixel-wise quantile to use when generating a substrate-only image."
                               "\nRecommended value: 0.5",
    "prep_substrate_sigma": "Width used in the Gaussian pre-filter. Generally, a value below 0.5 will reduce the "
                            "effectiveness of the rolling ball filter, while a value above 1.0 will reduce the "
                            "visibility of dimmer peaks.\nRecommended value: 0.75",
    "prep_substrate_radii": "Radii used for iterative rolling-ball background subtraction. Must be given in [square "
                            "brackets].\nRecommended value: [100, 30, 10]",
    "prep_sample_sigma": "Width used in the Gaussian pre-filter. Generally, a value below 0.5 will reduce the "
                         "effectiveness of the rolling ball filter, while a value above 1.0 will reduce the "
                         "visibility of dimmer peaks. \nRecommended value: 0.75",
    "prep_sample_radii": "Radii used for iterative rolling-ball background subtraction. Must be given in [square "
                         "brackets].\nRecommended value: [100, 30, 10]",
    "pkid_substrate_threshold": "Minimum value of a valid peak, relative to the standard deviation intensity of the "
                                "image. A good starting value is 0.2, but the ideal value varies by experiment. "
                                "Decrease to find more peaks, increase to find less.",
    "pkid_substrate_distance": "Minimum distance (in pixels) between valid peaks. If two otherwise valid peaks are "
                               "closer together, only the brighter of the two will be returned.",
    "pkid_mask_threshold": "Threshold used for masking substrate peaks. Setting this different from the substrate peak"
                           "threshold may lead to unexpected results.",
    "pkid_mask_dilation": "Radius of binary dilation applied to the substrate mask \nRecommended value: 3.0",
    "pkid_sample_threshold": "Minimum value of a valid peak, relative to the standard deviation intensity of the "
                             "image. A good starting value is 0.2, but the ideal value varies by experiment. "
                             "Decrease to find more peaks, increase to find less.",
    "pkid_sample_distance": "Minimum distance (in pixels) between valid peaks. If two otherwise valid peaks are "
                            "closer together, only the brighter of the two will be returned.",
    "laue_substrate_mis_err": "Maximum misorientation for two Laue diffraction patterns to be considered the same.",
    "laue_substrate_tolerance": "Maximum angular error (peak-wise average) for indexed Laue patterns.",
    "laue_sample_mis_err": "Maximum misorientation for two Laue diffraction patterns to be considered the same.",
    "laue_sample_tolerance": "Maximum angular error (peak-wise average) for indexed Laue patterns.",
    "grain_tolerance": "Maximum misorientation for two Laue diffraction patterns to be considered the same.",
    "grain_threshold": "Maximum number of frames that a grain can appear in, used to filter out incorrectly indexed "
                       "substrate peaks that didn't get filtered out.",
    "twin_tolerance": "Maximum crystal misorientation when determining whether two grains are twins."
}

GEN_TEXT = ["beamline", "year", "exp_id", "scan", "alt_id", "spec_seq", "calibration", "substrate", "sample"]
GEN_BOOL = ["show_plots", "verbose"]
PREP_SUB = ["prep_substrate_quantile", "prep_substrate_sigma", "prep_substrate_radii"]
PREP_SAM = ["prep_sample_sigma", "prep_sample_radii"]
PKID_SUB = ["pkid_substrate_threshold", "pkid_substrate_distance", "pkid_mask_threshold", "pkid_mask_dilation"]
PKID_SAM = ["pkid_sample_threshold", "pkid_sample_distance"]
# LAUE_SUB = ["laue_substrate_goodness", "laue_substrate_mis_err", "laue_substrate_tolerance", "laue_substrate_frequency",
#             "laue_substrate_comb_sub", "laue_substrate_times"]
# LAUE_SAM = ["laue_sample_goodness", "laue_sample_mis_err", "laue_sample_tolerance", "laue_sample_frequency",
#             "laue_sample_comb_sub", "laue_sample_times"]
LAUE_SUB = ["laue_substrate_mis_err", "laue_substrate_tolerance"]
LAUE_SAM = ["laue_sample_mis_err", "laue_sample_tolerance"]
LAUE_OTH = ["grain_tolerance", "grain_threshold", "twin_tolerance"]


class LaueApp:
    def __init__(self, conf_path=None):
        self.root = tk.Tk()
        self.root.title("LauePy")

        self.run_buttons = []
        self.run_frame = ttk.Frame(self.root)
        self.run_frame['padding'] = 10
        self.run_frame.grid(column=1, row=0)
        run_fcns = {
            "Prepare images": self.run_prep,
            "Find peaks": self.run_peaks,
            "Index peaks": self.run_auto,
            "Find grains": self.gen_macro,
        }
        for name, fcn in run_fcns.items():
            button = ttk.Button(self.run_frame, text=name, command=fcn)
            button.pack(side="top", fill='x', pady=5)
            self.run_buttons.append(button)

        self.cfg_buttons = []
        self.cfg_frame = ttk.Frame(self.root)
        self.cfg_frame['padding'] = 5
        self.cfg_frame.grid(column=0, row=1, columnspan=2, sticky="EW")
        cfg_fcns = {
            "Load": self.load_config,
            "Revert": self.revert_config,
            "Save": self.save_config,
        }
        for name, fcn in cfg_fcns.items():
            button = ttk.Button(self.cfg_frame, text=name, command=fcn)
            button.pack(side="left", fill='x', padx=5, expand=True)
            self.cfg_buttons.append(button)

        if conf_path is None:
            conf_path = f"{Path(__file__).parents[2]}/config_example/config.yml"
        self.conf_path = conf_path
        self.config = {key: tk.StringVar(self.root, value="") for key in PARAMS.keys()}
        for val in self.config.values():
            val.trace_add("write", self.disable_buttons)
        self.load_config(conf_path)

        self.input_nb = ttk.Notebook(self.root)
        self.input_nb.grid(column=0, row=0)
        tab_names = ["General", "Image prep", "Peak finding", "Laue indexing"]

        # Create the tabs on the main window
        for j, name in enumerate(tab_names):
            frame = ttk.Frame(self.input_nb)
            frame['padding'] = 10
            frame.grid(column=0, row=0)
            frame.columnconfigure(0, weight=1)
            frame.columnconfigure(1, weight=1)
            line = 0

            if j == 0:  # Create the general tab
                line = self.entry_section(frame, None, GEN_TEXT, line)
                for i, key in enumerate(GEN_BOOL):
                    val = self.config[key]
                    ttk.Checkbutton(
                        frame, variable=val, onvalue="True", offvalue="False", text=PARAMS[key]
                    ).grid(column=i, row=line, pady=4)

            if j == 1:  # Create the image prep tab
                line = self.entry_section(frame, "Substrate", PREP_SUB, line)
                self.entry_section(frame, "Sample", PREP_SAM, line)

            if j == 2:  # Create the peak finding tab
                line = self.entry_section(frame, "Substrate", PKID_SUB, line)
                self.entry_section(frame, "Sample", PKID_SAM, line)

            if j == 3:  # Create the Laue indexing tab
                line = self.entry_section(frame, "Substrate", LAUE_SUB, line)
                line = self.entry_section(frame, "Sample", LAUE_SAM, line)
                self.entry_section(frame, "Other", LAUE_OTH, line)

            self.input_nb.add(frame, text=name)

        self.root.mainloop()

    def disable_buttons(self, *_):
        for btn in self.run_buttons:
            btn.state(['disabled'])
        for btn in self.cfg_buttons[1:]:
            btn.state(['!disabled'])

    def save_config(self):
        [year, exp_id, scan, alt_id] = [self.config[key].get() for key in ("year", "exp_id", "scan", "alt_id")]
        if alt_id == "None":
            alt_id = ""
        self.conf_path = f"/home/beams/CXDUSER/34idc-work/{year}/{exp_id}/Analysis/lauepy_output/scan_{int(scan):04}" \
                         f"{alt_id}/config.yml"
        try:
            ut.save_config({key: str(val.get()) for key, val in self.config.items()}, self.conf_path)
            ut.read_config(self.conf_path)
        except FileNotFoundError as error:
            print()
            print(str(error))
            return
        for btn in self.run_buttons:
            btn.state(['!disabled'])
        for btn in self.cfg_buttons[1:]:
            btn.state(['disabled'])

    def load_config(self, conf_path=None):
        if conf_path is None:
            conf_path = filedialog.askopenfilename()
        self.conf_path = conf_path
        with open(conf_path, 'r') as f:
            config = yaml.safe_load(f)
        for key in self.config.keys():
            self.config[key].set(f"{config[key]}")
        self.save_config()

    def revert_config(self):
        self.load_config(self.conf_path)

    def run_prep(self):
        cfg = ut.read_config(self.conf_path)
        prep.extract_substrate(cfg)
        prep.cleanup_images(cfg)

    def run_peaks(self):
        cfg = ut.read_config(self.conf_path)
        peak_dict = pk.find_substrate_peaks(cfg, {})
        peak_dict = pk.find_sample_peaks(cfg, peak_dict)
        peak_dict = pk.record_positions(cfg, peak_dict)
        pk.save_peaks(cfg, peak_dict)

    def run_auto(self):
        cfg = ut.read_config(self.conf_path)
        for p in Path(f"{cfg['working_dir']}/grains").iterdir():
            p.unlink()
        sim = al.AutoLaue(cfg)
        sim.index()

    def gen_macro(self):
        cfg = ut.read_config(self.conf_path)
        grain.make_grain_dict(cfg)
        cfg = ut.read_config(self.conf_path)
        if twins.find_possible_twins(cfg):
            twins.find_twins(cfg)
            twins.cleanup_directory()

    def get_calibration(self, _):
        f = filedialog.askopenfilename()
        print(f)
        self.config["calibration"].set(f)
        self.root.focus_set()

    def entry_section(self, parent, heading, entries, line):
        if heading is not None:
            ttk.Label(parent, text=heading, font=("Arial", 14)).grid(column=0, row=line, sticky=tk.SW, pady=4)
            line += 1
        for key in entries:
            val = self.config[key]
            label = ttk.Label(parent, text=f"{PARAMS[key]}")
            label.grid(column=0, row=line, sticky=tk.W)
            CreateToolTip(label, TIPS[key])
            enter = ttk.Entry(parent, textvariable=val)
            enter.grid(column=1, row=line, sticky=tk.E, pady=4)
            if key == "calibration":
                enter.bind("<Button>", self.get_calibration)
            line += 1
        return line


class CreateToolTip(object):
    """
    create a tooltip for a given widget
    """
    def __init__(self, widget, text='widget info'):
        self.waittime = 500     #miliseconds
        self.wraplength = 300   #pixels
        self.widget = widget
        self.text = text
        self.widget.bind("<Enter>", self.enter)
        self.widget.bind("<Leave>", self.leave)
        self.widget.bind("<ButtonPress>", self.leave)
        self.id = None
        self.tw = None

    def enter(self, event=None):
        self.schedule()

    def leave(self, event=None):
        self.unschedule()
        self.hidetip()

    def schedule(self):
        self.unschedule()
        self.id = self.widget.after(self.waittime, self.showtip)

    def unschedule(self):
        id = self.id
        self.id = None
        if id:
            self.widget.after_cancel(id)

    def showtip(self, event=None):
        x = y = 0
        x, y, cx, cy = self.widget.bbox("insert")
        x += self.widget.winfo_rootx() + 25
        y += self.widget.winfo_rooty() + 20
        # creates a toplevel window
        self.tw = tk.Toplevel(self.widget)
        # Leaves only the label and removes the app window
        self.tw.wm_overrideredirect(True)
        self.tw.wm_geometry("+%d+%d" % (x, y))
        label = tk.Label(self.tw, text=self.text, justify='left', background="#ffffff", relief='solid', borderwidth=1,
                         wraplength=self.wraplength)
        label.pack(ipadx=1)

    def hidetip(self):
        tw = self.tw
        self.tw= None
        if tw:
            tw.destroy()


if __name__ == '__main__':
    LaueApp("/home/beams/CXDUSER/34idc-work/2022/LauePUP422/Analysis/lauepy_output/scan_0939/config.yml")
