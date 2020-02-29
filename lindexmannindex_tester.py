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

def calcindex(positions):
    pass

if __name__ == "__main__":
    xtc = MDAnalysis.coordinates.XTC.XTCReader("../gromacs/traj_comp.xtc")
    #xyz = np.ndarray((xtc.n_frames, 3), dtype=np.float32)
    i=0
    #indexes = np.ndarray((xtc.n_frames-1, 2))
    to = 5000
    indexes = np.zeros((3000,2))
    for ts in xtc:
        if ts.frame < 2000:
            continue
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
            
        square_average = 0
        for a in distances:
            square_average += a**2                                #<rij2>
        square_average = square_average / ts.positions.shape[0]
        average = np.sum(distances)/distances.shape[0]
        average_squared = average ** 2                            #<rij>2
        RMS = square_average - average_squared
        RMS = np.sqrt(RMS)
        right = RMS/average
        left = 2 / (ts.positions.shape[0] * (ts.positions.shape[0] - 1))
        LindemannIndex = left * right
        
        indexes[ts.frame-2000] = [ts.frame-2000, LindemannIndex]
        #print(LindemannIndex, indexes[ts.frame])
        if ts.frame == to-1:
            break
    
plt.plot(indexes[:,0], indexes[:,1])
plt.xlabel("Time (ns)")
plt.ylabel("Lindemann Index")
#plt.yscale("log")
plt.title("100 Argon atoms undergoing heating")