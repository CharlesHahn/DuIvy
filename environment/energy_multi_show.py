# author : charlie

import sys
import matplotlib
from matplotlib import pyplot as plt 
matplotlib.style.use('ggplot')



def xvg_deal(filename):
	with open(filename, 'r') as fo:
		content = fo.read()

	line_list = content.strip('\n').split('\n')
	column_num = len(line_list[-1].split())
	data = []
	for i in range(column_num):
		data.append([])

	title_lis = []
	for line in line_list:
		if line[0] != "#":
			if 'time' in line: 
				for title in line.split(' '):
					if title == "":
						continue
					if '(' == title[0]:
						title_lis[-1] += title
					else:
						title_lis.append(title) 

				for i in range(column_num):
					data[i].append( title_lis[i] )

			else:
				line_split = line.split()
				if len(line_split) == column_num:
					for i in range( column_num ):
						data[i].append( float(line_split[i]) )
				else:
					print("len(line_split) != column_num")

	return data



def multi_plot(data_lis, select, filename_lis):
	number_list = []
	for i in range(2, len(select)):
		number_list.append(int( select[i] ))

	for i in range(len(number_list)):
		ax = plt.subplot(len(number_list), 1, i+1 )
		ax_legend = []
		for data in data_lis:
			ax.plot( data[0][1:], data[ number_list[i] ][1:] ) 
			ax_legend.append( data[number_list[i]][0] )
			# print( data[ number_list[i] ][0:10] )
		for f in range(len(filename_lis)):
			filename = filename_lis[f].split(".")[0]
			ax_legend[f] += ' of ' + filename
		ax.legend(ax_legend)
		plt.ylabel("energy (kJ/mol)") 
		if i == 0:
			plt.title("multi energy comparison")

	plt.xlabel("time (ps)")
	plt.show()




def main():
	# energy_multi_show.py pro.xvg
	cmds = [ i for i in sys.argv[1:] ]

	data_lis = []
	filename_lis = []
	for cmd in cmds:
		if '-n' in cmd:
			column_select = cmd
		else:
			data_lis.append( xvg_deal(cmd) )
			filename_lis.append( cmd )
	
	multi_plot(data_lis, column_select, filename_lis)
	print(" Over ~ ")



if __name__ == '__main__':
	main()
