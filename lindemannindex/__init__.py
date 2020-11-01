# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 15:46:59 2020

@author: Alexander van Teijlingen
"""

import numpy as np
import MDAnalysis
from scipy.spatial.distance import cdist 


def calc_over_time(u, name):
    mols_of_interest = u.select_atoms("name "+name)
    times = np.ndarray((u.trajectory.n_frames,), dtype=np.float16)
    LindemannIndexes = []
    i = 0
    for ts in u.trajectory:
        times[i] = ts.time
        x = cdist(mols_of_interest.positions, mols_of_interest.positions)
        print(mols_of_interest.positions.shape)
    return x