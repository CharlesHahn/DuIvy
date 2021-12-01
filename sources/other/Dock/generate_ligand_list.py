# author : charlie

import os

path = os.getcwd()
print('-> CWD: ' + path)
file_list = os.listdir()
for file in file_list:
	if '.mol2' in file:
		with open("ligand_list.txt", 'a', encoding = 'UTF-8') as fo:
			fo.write(file + '\n')

print('-> ligand_list.txt generated fine! ')