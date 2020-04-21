# author : charlie

#################################################
# 通过 E = prolig - pro - lig 求得能量
# 要求三个文件的数据列数相同，且每列元素相同
#################################################

import os
import sys
from xvgshow import picture_subplot, xvg_deal, picture_oneplot

def energy_compute():
	try:
		pro_lig_file = sys.argv[1]
		pro_file = sys.argv[2]
		lig_file = sys.argv[3]
	except:
		print("No input filename! ")
		return 
	else:
		plotmode = '-s'
		try:
			if sys.argv[4] == '-o' or sys.argv[4] == '-s':
				plotmode = sys.argv[4]
			else:
				print("wrong plotmode, using '-s' ")
		except:
			pass 

		if os.path.exists(pro_lig_file):
			pro_lig_title, pro_lig_xlabel, pro_lig_ylabel, pro_lig_data = xvg_deal(pro_lig_file)
		else:
			print("pro_lig_file not exists in this directory ")
			return 
		if os.path.exists(pro_file):
			pro_title, pro_xlabel, pro_ylabel, pro_data = xvg_deal(pro_file)
		else:
			print("pro_file not exists in this directory ")
			return 
		if os.path.exists(lig_file):
			lig_title, lig_xlabel, lig_ylabel, lig_data = xvg_deal(lig_file)
		else:
			print("lig_file not exists in this directory ")
			return 

		final_data = pro_lig_data
		for i in range(1, len( final_data)):
			for j in range(1, len(final_data[i])):
				final_data[i][j] -= pro_data[i][j] + lig_data[i][j]
				# print(final_data[i][j], pro_data[i][j], lig_data[i][j])

		final_data.append([0 for x in range(len(final_data[0]))])
		final_data[-1][0] = "ETOTAL"
		for j in range(1, len(final_data[0])):
			for i in range(1, len(final_data)-1):
				final_data[-1][j] += final_data[i][j]

		with open("energy_compute_results.xvg", 'w') as fo:
			fo.write("##energy_compute_results.xvg generated from ")
			fo.write(pro_lig_file + ', ' + pro_file + ' and ' + lig_file + '\n')
			for j in range(len(final_data[0])):
				line_str = ' '
				for i in range(len(final_data)):
					line_str += str(final_data[i][j]) + '  ' 
				fo.write(line_str + '\n')

		cols = [ i for i in range(1, len(final_data))]
		pro_lig_title += " generated from energy_compute.py "
		if plotmode == '-s':
			picture_subplot(pro_lig_title, pro_lig_xlabel, pro_lig_ylabel, final_data, cols)
		elif plotmode == '-o':
			picture_oneplot(pro_lig_title, pro_lig_xlabel, pro_lig_ylabel, final_data, cols)




if __name__ == '__main__':
	energy_compute()


