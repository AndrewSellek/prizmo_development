import shutil
import numpy as np
from prizmo_commons import data_dir


def prepare():
    shutil.copyfile("../data/CO_cooling/cooling.dat", data_dir+"CO_cooling.dat")

    # data = np.loadtxt(data_dir+"CO_cooling.dat").T
