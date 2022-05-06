import tkinter as tk
from tkinter import ttk, filedialog
from pathlib import Path

import yaml

from lauepy.laue import utils as ut
from lauepy.laue import image_prep as prep
from lauepy.laue import peaks as pk
import lauepy.laue.auto_laue as al
import lauepy.laue.make_grain_dict as grain
import lauepy.laue.find_twins as twins


CONFIG_ENTRIES = {
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
    "prep_quantile": "Quantile filter",
    "prep_gamma": "Gamma",
    "prep_rb_radius": "Rolling ball radius",
    "prep_zero_fraction": "Zero point quantile",
    "prep_coefficient": "Overall brightness adjust",
    "prep_gaussian_sigma": "Gaussian filter sigma",
    "pkid_substrate_threshold": "Threshold",
    "pkid_substrate_distance": "Min. distance",
    "pkid_sample_threshold": "Threshold",
    "pkid_sample_distance": "Min. distance",
    "laue_substrate_goodness": "Goodness",
    "laue_substrate_mis_err": "Misorientation",
    "laue_substrate_tolerance": "Error tolerance",
    "laue_substrate_frequency": "Frequency",
    "laue_substrate_comb_sub": "Comb sub",
    "laue_substrate_times": "Times",
    "laue_sample_goodness": "Goodness",
    "laue_sample_mis_err": "Misorientation",
    "laue_sample_tolerance": "Error tolerance",
    "laue_sample_frequency": "Frequency",
    "laue_sample_comb_sub": "Comb sub",
    "laue_sample_times": "Times",
    "grain_tolerance": "Tolerance",
    "grain_threshold": "Threshold",
    "twin_tolerance": "Tolerance"
}

GENERAL = ["scan", "alt_id", "year", "exp_id", "spec_seq", "beamline", "calibration", "substrate", "sample",
           "show_plots", "verbose"]
PREP_SUB = ["prep_quantile", "prep_rb_radius"]
PREP_SAM = ["prep_zero_fraction", "prep_gamma", "prep_coefficient", "prep_gaussian_sigma"]
PKID_SUB = ["pkid_substrate_threshold", "pkid_substrate_distance"]
PKID_SAM = ["pkid_sample_threshold", "pkid_sample_distance"]
LAUE_SUB = ["laue_substrate_goodness", "laue_substrate_mis_err", "laue_substrate_tolerance", "laue_substrate_frequency",
            "laue_substrate_comb_sub", "laue_substrate_times"]
LAUE_SAM = ["laue_sample_goodness", "laue_sample_mis_err", "laue_sample_tolerance", "laue_sample_frequency",
            "laue_sample_comb_sub", "laue_sample_times"]
LAUE_OTH = ["grain_tolerance", "grain_threshold", "twin_tolerance"]

BOOL_VARS = ["show_plots", "verbose"]
STR_VARS = ["scan", "alt_id", "year", "exp_id", "spec_seq", "beamline", "calibration", "substrate", "sample"]
NUM_VARS = ["prep_quantile", "prep_rb_radius", "prep_zero_fraction", "prep_gamma", "prep_coefficient",
            "prep_gaussian_sigma", "pkid_substrate_threshold", "pkid_substrate_distance", "pkid_sample_threshold",
            "pkid_sample_distance", "laue_substrate_goodness", "laue_substrate_mis_err", "laue_substrate_tolerance",
            "laue_substrate_frequency", "laue_substrate_comb_sub", "laue_substrate_times", "laue_sample_goodness",
            "laue_sample_mis_err", "laue_sample_tolerance", "laue_sample_frequency", "laue_sample_comb_sub",
            "laue_sample_times", "grain_tolerance", "grain_threshold", "twin_tolerance"]


def init_entries(root):
    frame = ttk.Frame(root)
    frame['padding'] = 5, 10

    frame.columnconfigure(0, weight=1)
    frame.columnconfigure(1, weight=3)

    year = tk.StringVar()
    ttk.Label(frame, text="Year:").grid(column=0, row=0, sticky=tk.W)
    year_entry = ttk.Entry(frame, textvariable=year)
    year_entry.grid(column=1, row=0, sticky=tk.W)

    exp_id = tk.StringVar()
    ttk.Label(frame, text="Exp. directory:").grid(column=0, row=1, sticky=tk.W)
    exp_id_entry = ttk.Entry(frame, textvariable=exp_id)
    exp_id_entry.grid(column=1, row=1, sticky=tk.W)

    scan_num = tk.StringVar()
    ttk.Label(frame, text="Scan number:").grid(column=0, row=2, sticky=tk.W)
    exp_id_entry = ttk.Entry(frame, textvariable=scan_num)
    exp_id_entry.grid(column=1, row=2, sticky=tk.W)

    return frame, [year, exp_id, scan_num]


def init():
    root = tk.Tk()
    root.title("LAUEPY")

    root.rowconfigure(0)
    root.rowconfigure(1)

    entry_frame, entries = init_entries(root)
    entry_frame.grid(row=0, column=0)

    button_frame = ttk.Frame(root)
    button_frame['padding'] = 10
    button_frame.grid(row=1)

    button_frame.columnconfigure(0, weight=1)
    button_frame.columnconfigure(1, weight=1)

    def load_config():
        ret_val = filedialog.askopenfilename()
        root.quit()

    def new_config():
        ret_val = ut.new_analyis([x.get() for x in entries])
        root.quit()

    new_button = ttk.Button(button_frame, text="New", command=new_config)
    new_button.grid(column=0, row=0, padx=10)

    load_button = ttk.Button(button_frame, text="Load", command=load_config)
    load_button.grid(column=1, row=0, padx=10)

    root.mainloop()
    return


def main(conf_path=None):
    if conf_path is None:
        window = tk.Tk()
        window.withdraw()
        conf_path = filedialog.askopenfilename()
        window.destroy()

    with open(conf_path, 'r') as f:
        conf_orig = yaml.safe_load(f)
    buttons = []

    root = tk.Tk()
    root.title("LAUEPY")

    input_frame = ttk.Frame(root)
    input_frame['padding'] = 10
    input_frame.grid(column=0, row=0)
    input_frame.columnconfigure(0, weight=1)
    input_frame.columnconfigure(1, weight=3)

    config = {key: tk.StringVar(value=f"{val}") for key, val in conf_orig.items()}

    for i, (key, val) in enumerate(config.items()):
        ttk.Label(input_frame, text=f"{key}").grid(column=0, row=i, sticky=tk.W)
        if key in BOOL_VARS:
            entry_box = ttk.Checkbutton(input_frame, variable=val, onvalue="True", offvalue="False")
        else:
            entry_box = ttk.Entry(input_frame, textvariable=val)
        entry_box.grid(column=1, row=i)
        val.trace_add("write", lambda a, b, c: [btn.state(['disabled']) for btn in buttons[1:]])

    button_frame = ttk.Frame(root)
    button_frame['padding'] = 10
    button_frame.grid(column=1, row=0)

    def save_config():
        ut.save_config({key: val.get() for key, val in config.items()}, conf_path)
        for btn in buttons[1:]:
            btn.state(['!disabled'])

    def run_prep():
        cfg = ut.read_config(conf_path)
        prep.extract_substrate(cfg)
        prep.cleanup_images(cfg)

    def run_peaks():
        cfg = ut.read_config(conf_path)
        peak_dict = pk.find_substrate_peaks(cfg, {})
        peak_dict = pk.find_sample_peaks(cfg, peak_dict)
        peak_dict = pk.record_positions(cfg, peak_dict)
        pk.save_peaks(cfg, peak_dict)

    def run_auto():
        cfg = ut.read_config(conf_path)
        for p in Path(f"{cfg['working_dir']}/grains").iterdir():
            p.unlink()
        sim = al.AutoLaue(cfg)
        sim.index()

    def gen_macro():
        cfg = ut.read_config(conf_path)
        grain.make_grain_dict(cfg)

    def find_twins():
        cfg = ut.read_config(conf_path)
        if twins.find_possible_twins(cfg):
            twins.find_twins(cfg)
            twins.cleanup_directory()

    fcns = {
        "Save configuration": save_config,
        "Prepare images": run_prep,
        "Find peaks": run_peaks,
        "Index peaks": run_auto,
        "Generate macros": gen_macro,
        "Find twins": find_twins
    }

    for name, fcn in fcns.items():
        button = ttk.Button(button_frame, text=name, command=fcn)
        button.pack(fill='x', pady=20)
        buttons.append(button)

    root.mainloop()


if __name__ == '__main__':
    main()
