#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 13:58:22 2023

@author: nathan
"""

import tarfile
import os

os.chdir("/media/nathan/T7/path_integral/simulations/out/results")
tarballs = os.listdir()

for file in tarballs: 
    print(file)
    tar = tarfile.open(file)
    member = tar.getmembers()
    
    matching = [member[s].name for s in range(len(member)) if "trees" in member[s].name][0]
    
    tar.extract(matching, path = "../trees")