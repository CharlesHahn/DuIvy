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



def multi_plot(data_lis, select, filename_lis, title, ylabel, showlegend):
	number_list = []
	for i in range(2, len(select)):
		number_list.append(int( select[i] ))

	for i in range(len(number_list)):
		ax = plt.subplot(len(number_list), 1, i+1 )
		ax_legend = []
		for data in data_lis:
			ax.plot( [ d/1000 for d in data[0][1:] ], data[ number_list[i] ][1:] ) 
			ax_legend.append( data[number_list[i]][0] )
			# print( data[ number_list[i] ][0:10] )
		for f in range(len(filename_lis)):
			filename = filename_lis[f].split(".")[0]
			ax_legend[f] += ' of ' + filename
		if showlegend != 0:
			ax_legend = showlegend
		ax.legend(ax_legend)
		plt.ylabel(ylabel) 
		if i == 0:
			plt.title(title)


	plt.xlabel("time (ns)")
	plt.show()




def main():
	# energy_multi_show.py pro.xvg
	cmds = [ i for i in sys.argv[1:] ]

	print("energy_compare.py fileA, fileB, fileC, ... -n165, -ttitle_for_plot (optional) ")

	data_lis = []
	filename_lis = []
	title = "Multi Energy Comparison"
	ylabel = "energy (kJ/mol)"
	showlegend = 0
	for cmd in cmds:
		if '-n' == cmd[:2]:
			column_select = cmd
		elif '-t' == cmd[:2]:
			title = cmd[2: ].replace('_', ' ')
		elif '-y' == cmd[:2]:
			ylabel = cmd[2:].replace('_', ' ')
		elif '-l' == cmd[:2]:
			showlegend = cmd[2:].split('_')
		else:
			data_lis.append( xvg_deal(cmd) )
			filename_lis.append( cmd )
	
	multi_plot(data_lis, column_select, filename_lis, title, ylabel, showlegend)
	print(" Over ~ ")



if __name__ == '__main__':
	main()
