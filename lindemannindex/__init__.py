# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 15:46:59 2020

@author: Alexander van Teijlingen
"""

import numpy as np
import MDAnalysis

def calc_xtc(path):
    xtc = MDAnalysis.coordinates.XTC.XTCReader(path)
    xyz = []
    i=0
    for ts in xtc: #we really only need the first and last
       print(i, int(xtc.n_frames))
       print("time:", int(ts.frame), 'ns')
       for coord in ts.positions:
           xyz.append([float(coord[0]), float(coord[1]), float(coord[2])])
    #print(xyz)
    a="""LindemannIndex = self.calculator.calcindex(xyz, len(ts.positions))
    csv.write("{},{}\n".format(int(ts.frame), LindemannIndex))
    csv.flush()
    print((int(ts.frame), LindemannIndex))
    i+= 1
    frame +=1"""