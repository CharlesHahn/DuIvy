# author : charlie

import os 
import time
import matplotlib.pyplot as plt


def write_aminos_log(line):
	with open("bond_per_amino.aminolog", 'a') as fo:
		fo.write(line+'\n')


def write_ligands_log(line):
	with open("bond_per_ligand.ligandlog", 'a') as fo:
		fo.write(line+'\n')


def load_data(filename):
	with open(filename, 'r') as fo:
		file_content = fo.read()
	data = []
	file_lines = file_content.strip().strip('\n').split('\n')
	for line in file_lines:
		data.append( line.split(',') )
	return data


def average_per_amino(data):
	aminos = {}
	for line in data[1:]:
		key = line[1] + line[2]
		keys_list = aminos.keys()
		if key not in keys_list:
			aminos[key] = [1] + [float(i) for i in line[6: ]]

		elif key in keys_list:
			value = aminos[key]
			final_value = []
			final_value.append(value[0] + 1)
			if len(value[1:]) == len(line[6:]):
				for i in range(1, len(value)):
					final_value.append(value[i] + float(line[i+5]) )
				aminos[key] = final_value
			else:
				print("length of value[0:] if not equal to length of line[6:] ")
				write_aminos_log("length of value[0:] if not equal to length of line[6:] ")
		else:
			print("WRONG in average_per_amino function ! ")
			write_aminos_log("WRONG in average_per_amino function ! ")

	# print("amino -> [ count , "+ ", ".join(data[0][4:]) + " ]")
	write_aminos_log("amino -> [ count , "+ ", ".join(data[0][6:]) + " ]")

	for key, value in aminos.items():
		aminos[key] = [aminos[key][0]] + [float(i)/int(aminos[key][0]) for i in aminos[key][1: ]]
		# print(key, end = ' -> ')
		# print(["{:.4f}".format(i) for i in aminos[key]])
		write_aminos_log(key + ' -> ' + ',  '.join(['{:.4f}'.format(i) for i in aminos[key]] ) )
	return aminos 


def average_per_ligand(data):
	ligands = {}
	for line in data[1:]:
		key = line[4]
		keys_list = ligands.keys()
		if key not in keys_list:
			ligands[key] = [1] + [float(i) for i in line[6: ]]

		elif key in keys_list:
			value = ligands[key]
			final_value = []
			final_value.append(value[0] + 1)
			if len(value[1:]) == len(line[6:]):
				for i in range(1, len(value)):
					final_value.append(value[i] + float(line[i+5]) )
				ligands[key] = final_value
			else:
				print("length of value[0:] if not equal to length of line[6:] ")
				write_ligands_log("length of value[0:] if not equal to length of line[6:] ")
		else:
			print("WRONG in average_per_ligand function ! ")
			write_ligands_log("WRONG in average_per_ligand function ! ")

	# print("ligand -> [ count , "+ ", ".join(data[0][4:]) + " ]")
	write_ligands_log("ligand -> [ count , "+ ", ".join(data[0][6:]) + " ]")

	for key, value in ligands.items():
		ligands[key] = [ligands[key][0]] + [float(i)/int(ligands[key][0]) for i in ligands[key][1: ]]
		# print(key, end = ' -> ')
		# print(["{:.4f}".format(i) for i in ligands[key]])
		write_ligands_log(key + ' -> ' + ',  '.join(['{:.4f}'.format(i) for i in ligands[key]] ) )
	return ligands


def draw_data(aminos, picturename= "picturename"):
	x_names = ['17LEU', '18VAL', '19PHE', '20PHE', '21ALA', '22GLU', '23ASP', '24VAL', '25GLY', '26SER',\
			'27ASN', '28LYS', '29GLY', '30ALA', '31ILE', '32ILE', '33GLY', '34LEU', '35MET', '36VAL',\
			'37GLY', '38GLY', '39VAL', '40VAL','41ILE', '42ALA' ]
	draw_dic = {}
	blank = [0 for i in aminos[[k for k in aminos.keys()][0] ] ]
	for xname in x_names:
		draw_dic[xname] = blank

	for key, value in aminos.items():
		if  key in draw_dic.keys():
			draw_dic[key] = aminos[key]
		else: 
			print("draw_dic KEY WRONG !!! ")
			return 1

	x_data = [i for i in draw_dic.keys()]
	y_data = [i[0] for i in draw_dic.values()]
	plt.bar( x_data, y_data)
	
	for x, li in zip(aminos.keys(), aminos.values()):
		tag = ''
		for l in li[1:]:
			tag += '\n{:.2f}'.format(l)
		plt.text(x, 0.5*li[0], tag, ha='center', va='bottom', fontsize=8)

	plt.xticks(rotation=90)
	plt.title(picturename)
	plt.ylabel("num")
	plt.show()


def draw_ligand_data(ligands, picturename= "picturename"):
	x_data = [i for i in ligands.keys()]
	y_data = [i[0] for i in ligands.values()]
	plt.bar( x_data, y_data)
	
	for x, li in zip(ligands.keys(), ligands.values()):
		tag = ''
		for l in li[1:]:
			tag += '\n{:.2f}'.format(l)
		plt.text(x, 0.5*li[0], tag, ha='center', va='bottom', fontsize=8)

	plt.xticks(rotation=90)
	plt.title(picturename)
	plt.ylabel("num")
	plt.show()


def main():
	"""
	print("running analysis_deep.py")
	status = os.system("python analysis_deep.py")
	print("run analysis_deep.py -> "+ str(status))
	"""
	pwd = os.getcwd()
	with open("bond_per_amino.aminolog", 'w') as fo:
		fo.write("#### "+ pwd + " #### \n")
	with open("bond_per_ligand.ligandlog", 'w') as fo:
		fo.write("#### "+ pwd + " #### \n")

	print("========== hydrophobic_data.charlielog ============")
	write_aminos_log("========== hydrophobic_data.charlielog ============")
	write_ligands_log("========== hydrophobic_data.charlielog ============")
	data = load_data("hydrophobic_data.charlielog")
	# print("count = " + str(len(data)))
	write_aminos_log("count = " + str(len(data)))
	write_ligands_log("count = " + str(len(data)))
	aminos = average_per_amino(data)
	draw_data(aminos, "hydrophobic_data.charlielog")
	ligands = average_per_ligand(data)
	draw_ligand_data(ligands, "hydrophobic_data.charlielog, ligand")
	
	print("========== hydrogen_bonds_data.charlielog =========")
	write_aminos_log("========== hydrogen_bonds_data.charlielog =========")
	write_ligands_log("========== hydrogen_bonds_data.charlielog =========")
	data = load_data("hydrogen_bonds_data.charlielog")
	# print("count = " + str(len(data)))
	write_aminos_log("count = " + str(len(data)))
	write_ligands_log("count = " + str(len(data)))
	aminos = average_per_amino(data)
	draw_data(aminos, "hydrogen_bonds_data.charlielog")
	ligands = average_per_ligand(data)
	draw_ligand_data(ligands, "hydrogen_bonds_data.charlielog, ligand")

	print("========== pi_stacks_data.charlielog ==============")
	write_aminos_log("========== pi_stacks_data.charlielog ==============")
	write_ligands_log("========== pi_stacks_data.charlielog ==============")
	data = load_data("pi_stacks_data.charlielog")
	# print("count = " + str(len(data)))
	write_aminos_log("count = " + str(len(data)))
	write_ligands_log("count = " + str(len(data)))
	aminos = average_per_amino(data)
	draw_data(aminos, "pi_stacks_data.charlielog")
	ligands = average_per_ligand(data)
	draw_ligand_data(ligands, "pi_stacks_data.charlielog, ligand")


	print("run average_perAA over !")

	time.sleep(6)



if __name__ == '__main__':
	main()

