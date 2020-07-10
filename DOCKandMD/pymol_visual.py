# author : charlie

import subprocess

command = "protein.pdb"

for i in range(int(input(" ligand range => "))):
    command += " python_extract_model_"+str(i)+".pdb"

# print(command)

subprocess.getstatusoutput("pymol " + command)


