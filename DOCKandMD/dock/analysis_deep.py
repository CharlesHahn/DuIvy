# author: charlie

import os
import subprocess
import xml.etree.ElementTree as ET
from matplotlib import pyplot as plt
import seaborn as sns
sns.set_style("darkgrid")


def write_logfile(line):
	with open("analysis_deep.runlog", 'a') as fo:
		fo.write(line + '\n')
	return 0


def hydrophobic_yield():
	for i in range(10000):
		yield i 

def hydrogen_yield():
	for i in range(10000):
		yield i 

def pi_stack_yield():
	for i in range(1000):
		yield i 


def data_dump(lis_data, filename):
	lines = ''
	for lis in lis_data:
		lines += ','.join([str(li) for li in lis])
		lines += '\n'
	with open(filename, 'w') as fo:
		fo.write(lines.strip('\n'))
	return 0


def make_complex(protein_file, ligand_list):
	complex_files = []
	for ligand_file in ligand_list:
		with open(ligand_file, 'r') as fo:
			contents = fo.read()
		ligand_content = []
		for line in contents.strip().split('\n'):
			if line[:4] == 'ATOM':
				ligand_content.append(line)
		with open(protein_file, 'r') as fo:
			protein_content = fo.read()
		complex_content = protein_content + '\n' + '\n'.join(ligand_content)
		complex_file = 'complex_' + ligand_file
		with open(complex_file, 'w') as fo:
			fo.write(complex_content)
		complex_files.append(complex_file)
	return complex_files 


def plip(complexs_list):
	xml_reports = []
	progro = 0
	allwork = len(complexs_list)
	for complex_file in complexs_list:
		plip_command = "D:\\CharlieAPP\\plip-stable\\plip\\plipcmd.py -f " + complex_file + " -yxv"
		status, output = subprocess.getstatusoutput(plip_command)
		## print("plip status -> " + str(status))
		# print("plip output -> " + output)
		status, output = subprocess.getstatusoutput("ren report.xml " + complex_file[: -4] + '.xml')
		## print("rename status -> " + str(status))
		# print("rename output -> " + output)
		if status == 1:
			del_status, del_output = subprocess.getstatusoutput("del " + complex_file[: -4] + '.xml')
			if del_status == 1:
				print("Trouble in delete " + complex_file[: -4] + '.xml, CHECK IT manually !' )
				print("=====================================")
			elif del_status == 0:
				## print("try to delete this motherfucker !")
				re_status, re_output = subprocess.getstatusoutput("ren report.xml " + complex_file[: -4] + '.xml')
				## print("RE rename -----> " + str(re_status))
				if re_status == 1:
					print("Trouble in re_rename report.xml" + complex_file)
					print("=====================================")

		xml_reports.append(complex_file[:-4] + '.xml')

		progro += 1 
		progrocess = progro/allwork
		print("PLIP progrocess : {: <20}{:.2%}".format( "*"*int( progrocess*20 ), progrocess ), end = '\r')

	return xml_reports


def xml_analysis(xml_reports):
	hydrophobic_data = []
	hydrogen_bonds_data = []
	pi_stacks_data = []
	hydrophobic_data.append( 
		['num', 'resnr', 'restype', 'reschain', 'ligcarbonidx','protcarbonidx','dist'] )
	hydrogen_bonds_data.append( 
		['num','resnr','restype','reschain','lig_hbidx', 'lig_hbtype','dist_h-a','dist_d-a','don_angle'])
	pi_stacks_data.append( 
		['num', 'resnr', 'restypre', 'reschain','lig_pi_idx1', 'lig_pi_idx2','centdist','angle','offset'])

	hydrophobic_yie = hydrophobic_yield()
	hydrogen_yie = hydrogen_yield()
	pi_stack_yie = pi_stack_yield()

	hydrophobic_count = 0
	hydrogen_bonds_count = 0
	pi_stacks_count = 0

	for xml_file in xml_reports:
		tree = ET.parse(xml_file)
		root = tree.getroot()
		hydrophobic_interactions = root[1][4][0] #hydrophobic_interactions
		hydrogen_bonds = root[1][4][1] # hydrogen_bonds
		pi_stacks = root[1][4][4] # pi_stacks
		## print("=====> in " + xml_file + " <=====")
		## print('==> ' + hydrophobic_interactions.tag + " -> " + str(len(hydrophobic_interactions)))
		## print('==> ' + hydrogen_bonds.tag + " -> " + str(len(hydrogen_bonds)))
		## print('==> ' + pi_stacks.tag + ' -> ' + str(len(pi_stacks)))
		hydrophobic_count += len(hydrophobic_interactions)
		hydrogen_bonds_count += len(hydrogen_bonds)
		pi_stacks_count += len(pi_stacks)

		if len(hydrophobic_interactions) > 0:
			for hydrophobic in hydrophobic_interactions:
				resnr = hydrophobic.find('resnr').text
				restype = hydrophobic.find('restype').text
				reschain = hydrophobic.find('reschain').text
				dist = hydrophobic.find('dist').text
				ligcarbonidx = hydrophobic.find('ligcarbonidx').text
				protcarbonidx = hydrophobic.find('protcarbonidx').text
				# print(resnr, restype, reschain, dist)
				hydrophobic_data.append(
					[next(hydrophobic_yie), resnr, restype, reschain, ligcarbonidx,protcarbonidx, dist])

		if len(hydrogen_bonds) > 0:
			for hydrogen_bond in hydrogen_bonds:
				resnr = hydrogen_bond.find('resnr').text
				restype = hydrogen_bond.find('restype').text
				reschain = hydrogen_bond.find('reschain').text
				dist_h_a = hydrogen_bond.find('dist_h-a').text
				dist_d_a = hydrogen_bond.find('dist_d-a').text
				don_angle = hydrogen_bond.find('don_angle').text

				if hydrogen_bond.find("protisdon").text == "True":
					lig_hbidx = hydrogen_bond.find("acceptoridx").text
					lig_hbtype = hydrogen_bond.find("acceptortype").text + '_acceptor'
				elif hydrogen_bond.find("protisdon").text == 'False':
					lig_hbidx = hydrogen_bond.find("donoridx").text
					lig_hbtype = hydrogen_bond.find("donortype").text + '_donor'
				# print(resnr, restype, reschain,dist_h_a, dist_d_a, don_angle)
				hydrogen_bonds_data.append(
					[next(hydrogen_yie),resnr,restype,reschain,lig_hbidx, lig_hbtype,dist_h_a,dist_d_a,don_angle])

		if len(pi_stacks) > 0:
			for pi_stack in pi_stacks:
				resnr = pi_stack.find('resnr').text
				restype = pi_stack.find('restype').text
				reschain = pi_stack.find('reschain').text
				centdist = pi_stack.find('centdist').text
				angle = pi_stack.find('angle').text
				offset = pi_stack.find("offset").text
				lig_pi_idx1 = pi_stack.find("lig_idx_list")[0].text
				lig_pi_idx2 = pi_stack.find("lig_idx_list")[1].text
				pi_stacks_data.append( 
					[next(pi_stack_yie), resnr, restype, reschain, lig_pi_idx1, lig_pi_idx2, centdist, angle, offset ])

	status = data_dump(hydrophobic_data, 'hydrophobic_data.charlielog')
	if status == 0:
		print("hydrophobic data dump -> " + str(hydrophobic_count))
		write_logfile("hydrophobic data dump -> " + str(hydrophobic_count))
	else:
		print("hydrophobic data dump ERROR ! ")
	status = data_dump(hydrogen_bonds_data, 'hydrogen_bonds_data.charlielog')
	if status == 0:
		print("hydrogen_bonds data dump -> " + str(hydrogen_bonds_count))
		write_logfile("hydrogen_bonds data dump -> " + str(hydrogen_bonds_count))
	else:
		print("hydrogen_bonds data dump ERROR !")
	status = data_dump(pi_stacks_data, 'pi_stacks_data.charlielog')
	if status == 0:
		print("pi_stacks data dump -> " + str(pi_stacks_count))
		write_logfile("pi_stacks data dump -> " + str(pi_stacks_count))
	else:
		print("pi_stacks data dump ERROR !")

	return hydrophobic_data, hydrogen_bonds_data, pi_stacks_data 


def file_cp():
	status, output = subprocess.getstatusoutput("copy ..\\protein.pdb protein.pdb")
	## print("protein copy -> " + str(status))

	ligand_list = []
	for i in range(10):
		cp_command = "copy ..\\python_extract_model_" + str(i) + ".pdb ligand_" + str(i) + ".pdb"
		status, output = subprocess.getstatusoutput(cp_command)
		## print("ligand_" + str(i) + " copy -> " + str(status))
		ligand_list.append('ligand_' + str(i) + '.pdb')
	return ligand_list


def draw_data(bond_data, picturename):
	x_names = ['17LEU', '18VAL', '19PHE', '20PHE', '21ALA', '22GLU', '23ASP', '24VAL', '25GLY', '26SER',\
			'27ASN', '28LYS', '29GLY', '30ALA', '31ILE', '32ILE', '33GLY', '34LEU', '35MET', '36VAL',\
			'37GLY', '38GLY', '39VAL', '40VAL','41ILE', '42ALA' ]
	draw_dic = {}
	for xname in x_names:
		draw_dic[xname] = 0
	for lis in bond_data[1:]:
		dic_key = str(lis[1].strip()) + str(lis[2].strip())
		if  dic_key in draw_dic.keys():
			draw_dic[dic_key] += 1
		else: 
			print("draw_dic KEY WRONG !!! ")
			return 1
	x_data = [i for i in draw_dic.keys()]
	y_data = [i for i in draw_dic.values()]
	plt.bar( x_data, y_data)
	plt.xticks(rotation=90)
	plt.title(picturename)
	plt.ylabel("num")
	plt.show()

	lig_dic = {}
	for lis in bond_data[1:]:
		lig_key = str(lis[4])
		if lig_key in lig_dic.keys():
			lig_dic[lig_key] += 1
		elif lig_key not in lig_dic.keys():
			lig_dic[lig_key] = 1 
	lig_xdata = [i for i in lig_dic.keys()]
	lig_ydata = [i for i in lig_dic.values()]
	plt.bar( lig_xdata, lig_ydata)
	plt.xticks(rotation=90)
	plt.title(picturename + " ligand")
	plt.ylabel("num")
	plt.show()



def everage_data(bond_data, bond_name):
	## print("========> " + bond_name)
	write_logfile("========> " + bond_name )
	length = len(bond_data)
	if length > 1:
		for i in range(6, len(bond_data[0])):
			dataname = bond_data[0][i]
			total = 0 
			for j in range(1, len(bond_data)):
				total += float(bond_data[j][i])
			ave = total / ( length - 1 )
			## print(dataname + " --> " + str(ave))
			write_logfile(dataname + " --> " + str(ave))
	else:
		## print("length <= 1, unable to computer average ")
		write_logfile("length <= 1, unable to computer average ")
	return 0



def main():
	# initial the runlog file 
	with open("analysis_deep.runlog", 'w') as fo:
		fo.write('')

	ligand_list = file_cp()
	print("In analysis_deep.py -> cp over !")
	complexs_list = make_complex('protein.pdb', ligand_list)
	print("In analysis_deep.py -> make complexs over !")
	xml_reports = plip(complexs_list)
	print("In analysis_deep.py -> generate report.xml over ! ")
	hydrophobic_data, hydrogen_bonds_data, pi_stacks_data = xml_analysis( xml_reports )
	print("In analysis_deep.py -> generate bond_data finished ! ")
	## draw_data(hydrophobic_data, 'hydrophobic_interactions')
	## draw_data(hydrogen_bonds_data, 'hydrogen_bonds')
	## draw_data(pi_stacks_data, 'pi_stacks')
	everage_data(hydrophobic_data, 'hydrophobic_interactions')
	everage_data(hydrogen_bonds_data, 'hydrogen_bonds')
	everage_data(pi_stacks_data, 'pi_stacks')
	print("In analysis_deep.py -> generate everage_data over !")



if __name__ == '__main__':
	main()

		

