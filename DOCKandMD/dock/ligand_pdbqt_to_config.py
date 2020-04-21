# author : charlie

import os 
import subprocess

path = os.getcwd()
print("-> CWD: " + path)

file_list = os.listdir() # indead, i should use shell code to do this stupid tash!
for file in file_list:
	if (".pdbqt" in file) and ("protein" not in file):
		ligand_name = file.split('.')[0]
		with open(ligand_name + ".conf", 'w', encoding = 'UTF-8') as fo:
			fo.write("receptor = protein.pdbqt\nligand = " + ligand_name + ".pdbqt\n")
			fo.write("out = result_" + ligand_name + ".pdbqt\n")
			fo.write("center_x = 0\ncenter_y = 0\ncenter_z = -3\nsize_x = 45\nsize_y = 30\nsize_z = 30\n")
			fo.write("cpu = 4\nenergy_range = 4\nexhaustiveness = 8\nnum_modes = 9")

		print("-> " + ligand_name + " from pbdqt to config FINE!")
		status, output = subprocess.getstatusoutput("vina --config " + ligand_name + ".conf " + "--log " + ligand_name + "_vina.log")
		print("-> vina running on " + ligand_name + " is " + str(status))
		print(output)