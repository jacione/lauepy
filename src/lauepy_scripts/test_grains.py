from src.lauepy import utils as ut
from src.lauepy import make_grain_dict as mgd


if __name__ == "__main__":
    cfg = ut.read_config("/home/beams/CXDUSER/34idc-work/2022/LauePUP422/Analysis/lauepy_output/scan_0939/config.yml")
    mgd.cluster_grains(cfg)
