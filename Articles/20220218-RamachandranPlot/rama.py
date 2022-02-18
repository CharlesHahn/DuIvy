#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Gábor Erdős, 2017
Contact info: gerdos[at]caesar.elte.hu
The preferences were calculated from the following artice:
Lovell et al. Structure validation by Calpha geometry: phi,psi and Cbeta deviation. 2003
DOI: 10.1002/prot.10286
"""

import sys
from pyrama import calc_ramachandran, plot_ramachandran
if len(sys.argv) < 2:
    sys.exit("Usage: pyrama my_pdb_file.pdb")

normals, outliers = calc_ramachandran(sys.argv[1:])
with open("rama.log", "w") as fo: 
	fo.write("File : " + " ".join(sys.argv[1:]) + "\n")
	print(" === normals === ")
	fo.write(" === normals === \n")
	for key, values in normals.items():
		print(key, len(values["x"]), len(values["y"]))
		fo.write("{} {} {}".format(key, len(values["x"]), len(values["y"])) + "\n")
	print(" === outliers === ")
	fo.write(" === outliers === \n")
	for key, values in outliers.items():
		print(key, len(values["x"]), len(values["y"]))
		fo.write("{} {} {}".format(key, len(values["x"]), len(values["y"])) + "\n")
plot_ramachandran(normals, outliers)
