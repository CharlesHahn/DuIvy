#####################
# author : charlie
# time : 20200709
# usage : to compare some files which contains same columns
# command : 
#   python3 xvg_compare.py file1 file2 file3 file... 
#       -y for ylabel (e.g. -ylabel_1,label_2,label_3,...)
#       -n for colomns selected for draw (e.g. -n56)
#       -t for plot title (e.g. -tEnergy_for_Protein_and_Ligands)
#       -x for xlabel (e.g. -xTime_(ns))
# notice for command:
#   all "_" will be replaced by space
#   all "," will be set as a symbel to split
#   all item formed by split() will show one by one with subplots
#####################

import sys
import matplotlib
from matplotlib import pyplot as plt 
from matplotlib import pylab as pylab

# 绘图控制参数
myparams = {
    'axes.labelsize': '12',
    'xtick.labelsize': '12',
    'ytick.labelsize': '12',
    'lines.linewidth': '1',
    'legend.fontsize': '12',
    'font.family': 'Times New Roman'
}
pylab.rcParams.update(myparams)
# matplotlib.style.use('ggplot')


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
        if line[0] != "#" and line[0] != '@':
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
                    print(" WRONG! len(line_split) != column_num")

    return data


def yield_ylabel( num ):
    for y in range(num):
        yield y


def multi_plot(data_lis, select, filename_lis, title, ylabel, xlabel, showlegend):
    number_list = []
    ylabel_list = ylabel
    ylabel_index = yield_ylabel(len(ylabel_list))
    for i in range(2, len(select)):
        number_list.append(int( select[i] ))

    for i in range(len(number_list)):
        ax = plt.subplot(len(number_list), 1, i+1 )
        ax_legend = []
        for data in data_lis:
            ax.plot( [ d for d in data[0][1:] ], data[ number_list[i] ][1:] ) 
            ax_legend.append( data[number_list[i]][0] )
            # print( data[ number_list[i] ][0:10] )
        if len(ylabel_list) == len(number_list):
            ax.set_ylabel( ylabel_list[ next( ylabel_index )] )
        for f in range(len(filename_lis)):
            filename = filename_lis[f].split(".")[0]
            ax_legend[f] = str(ax_legend[f]) +  ' of ' + filename
        if showlegend != 0:
            ax_legend = showlegend
        # print(ax_legend)
        ax.legend(ax_legend).get_frame().set_linewidth(0.0) #载入图例并设置无边框
        if len(ylabel_list) != len(number_list):
            plt.ylabel(ylabel) 
        if i == 0:
            plt.title(title)


    plt.xlabel(xlabel)
    plt.show()




def main():
    # energy_multi_show.py pro.xvg
    cmds = [ i for i in sys.argv[1:] ]

    print("""xvg_compare.py fileA, fileB, fileC, ... -n165, 
        -tTitle -yy_label,ylable_2 -xxlabel -llegend_1,legend_2 """)

    data_lis = []
    filename_lis = []
    title = ""
    ylabel = "xvg ylabel"
    xlabel = 'xvg xlabel'
    showlegend = 0
    for cmd in cmds:
        if '-n' == cmd[:2]:
            column_select = cmd
        elif '-t' == cmd[:2]:
            title = cmd[2: ].replace('_', ' ')
        elif '-y' == cmd[:2]:
            ylabel = cmd[2:].replace('_', ' ').split(",")
        elif '-x' == cmd[:2]:
            xlabel = cmd[2:].replace("_", ' ')
        elif '-l' == cmd[:2]:
            showlegend = cmd[2:].replace('_', ' ').split(',')
        else:
            data_lis.append( xvg_deal(cmd) )
            filename_lis.append( cmd )
    
    multi_plot(data_lis, column_select, filename_lis, title, ylabel, xlabel, showlegend)
    print(" ~ Over ~ ")



if __name__ == '__main__':
    main()
