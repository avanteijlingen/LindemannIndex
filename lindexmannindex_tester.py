# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 21:24:56 2020

@author: avtei
"""
import numpy as np
import lindemannindex as li
import matplotlib.pyplot as plt

#li.calc_xtc("../gromacs/traj_comp.xtc")

import MDAnalysis

def calcindex(positions, distances):
    square_average = 0
    for a in distances:
        square_average += a**2                                #<rij2>
    square_average = square_average / positions.shape[0]
    average = np.sum(distances)/distances.shape[0]
    average_squared = average ** 2                            #<rij>2
    RMS = square_average - average_squared
    RMS = np.sqrt(RMS)
    right = RMS/average
    left = 2 / (positions.shape[0] * (positions.shape[0] - 1))
    Lindex = left * right
    return Lindex

def CalculateXYZ(xyzfilepath):
    pass

def CalculateOverTraj(xtcfilepath, verbose=False):
    xtc = MDAnalysis.coordinates.XTC.XTCReader(xtcfilepath)
    i=0
    indexes = np.ndarray((xtc.n_frames, 1))
    for ts in xtc:
        atom1, atom2 = 0, 1
        distances = np.ndarray((int(((ts.positions.shape[0]-1)*(ts.positions.shape[0]))/2),))
        i=0
        while atom1 < ts.positions.shape[0]-1:
            while atom2 < ts.positions.shape[0]:
                p1 = ts.positions[atom1]
                p2 = ts.positions[atom2]
                squared_dist = np.sum((p1-p2)**2, axis=0)
                dist = np.sqrt(squared_dist)
                distances[i] = dist
                atom2 += 1
                i+=1
            atom1 += 1
            atom2 = atom1 + 1

        Lindex = calcindex(ts.positions, distances)
        indexes[ts.frame] = Lindex
        print(Lindex, indexes[ts.frame])
    return indexes

Lindexes = CalculateOverTraj("../gromacs/traj_comp.xtc")
plt.plot([x/2 for x in range(0,200)], Lindexes)
plt.xlabel("Temperature (K)")
plt.ylabel("Lindemann Index")
plt.title("887 Argon atoms undergoing heating")