# author : charlie

import os
import re
import subprocess


def fetch_EFEB(ligand):
	energy = re.findall(r"Estimated Free Energy of Binding    =.*?kcal/mol  \[=", ligand )
	binding_energy = energy[0][37: -4]
	#print(binding_energy)
	binding_energy_num = float( binding_energy[:-8] )
	#print( binding_energy_num )
	return (binding_energy_num, binding_energy)


def fetch_ki(ligand):
	ki = re.findall(r"Estimated Inhibition Constant, Ki   =.*?Temperature.*?\]", ligand )
	if len(ki) == 0:
		inhibition_constant_num = 0
		inhibition_constant = "none"
	else:
		inhibition_constant = ki[0][37:-38]
		#print(inhibition_constant)
		inhibition_constant_num = float( inhibition_constant[ :9] )
		#print(inhibition_constant_num)
	return (inhibition_constant_num, inhibition_constant)


def write_data(data_lis):
	# data_lis should be a list containing the EFEB, ki
	# data_lis = [model_num, binding_energy_num,binding_energy,inhibition_constant_num,inhibition_constant]
	if data_lis[0] == -1:
		with open("models_data.csv", "w") as fo:
			fo.write(
				"name,binding_energy_num,binding_energy,inhibition_constant_num(um),inhibition_constant\n")
	else:
		with open("models_data.csv", "a") as fo:
			fo.write("model_" + ",".join([str(item) for item in data_lis] ) + "\n")


def write_pdb(content, num): 
	with open("python_extract_model_" + str(num) + ".pdbqt", 'w') as fo:
		fo.write(content)
	status, output = subprocess.getstatusoutput("cut -c-66 python_extract_model_"
		+str(num)+".pdbqt > python_extract_model_"+str(num)+".pdb")
	print("saving pdb status: " + str(status) )
	print("saving pdb output: " + output)
	print("rm python_extract_model_" + str(num) + ".pdbqt : ", end = "")
	os.remove("python_extract_model_" + str(num) + ".pdbqt")
	print("done")


def main():
	# convert dlg to pdbqt
	status, output = subprocess.getstatusoutput("grep '^DOCKED' dock.dlg | cut -c9- > results_dock.pdbqt")
	print("converting dlg to pdbqt : " + str(status))
	print(output)
	
	with open('results_dock.pdbqt', 'r') as fo:
	    content = fo.read()
	    
	models = re.findall(r"MODEL[\s\S]*?ENDMDL", content)
	# print( models[0] )

	models_lis = []

	for model in models:
		model_dic = {}
		binding_energy_num, binding_energy = fetch_EFEB(model)
		inhibition_constant_num, inhibition_constant = fetch_ki(model)
		model_dic["binding_energy_num"] = binding_energy_num
		model_dic["binding_energy"] = binding_energy 
		model_dic["inhibition_constant_num"] = (inhibition_constant_num,
			inhibition_constant_num*1000)["mM" in inhibition_constant]
		model_dic["inhibition_constant"] = inhibition_constant 
		model_dic["model_content"] = model 
		models_lis.append(model_dic)

	sorted_models_lis = sorted( models_lis, key = lambda item : item["binding_energy_num"])
	
	# initial the write_data()
	write_data([-1])

	for i in range(len(sorted_models_lis)):
		data_lis = [str(i), sorted_models_lis[i]["binding_energy_num"]]
		data_lis.append(sorted_models_lis[i]["binding_energy"])
		data_lis.append(sorted_models_lis[i]["inhibition_constant_num"])
		data_lis.append(sorted_models_lis[i]["inhibition_constant"])
		write_data(data_lis)
		# print(data_lis)
		write_pdb(sorted_models_lis[i]["model_content"], i)




main()
