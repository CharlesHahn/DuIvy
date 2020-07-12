#################################################
# author : charlie
# time : 20200712
# command : python3 energy_compute.py prolig.xvg pro.xvg lig.xvg -foutput
# usage :
# 用于计算蛋白配体之间的相互作用，计算方法参考Jerkwin博客：
#     https://jerkwin.github.io/2019/09/06/使用GROMACS计算分子间相互作用/
# 通过 E = prolig - pro - lig 求得能量
# 
# 输入一共需要三个文件：蛋白配体的能量xvg文件，蛋白的能量xvg文件，配体的能量xvg文件
# 每个输入文件需要包含四列能量数据，用于各项能量以及总库伦势能和总能量的计算
# 要求四列能量数据的顺序为 
### LJ-SR  | Disper.corr. | Coulomb-SR | Coul.-recip.
# 
# 生成数据文件的各列名如下
### LJ-SR  | Disper.corr. | Coulomb-SR | Coul.-recip. | ETOTAL | COULOMB
#################################################

import os
import sys
import seaborn as sns
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


def cols_num_gen(cols):
    for i in range(1, len(cols)+1):
        yield i


def picture_oneplot(title, xlabel, ylabel, data, cols):
    legend_lis = []
    for i in cols:
        legend_lis.append(data[i][0])
        plt.plot(data[0][1:], data[i][1:])

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(labels=legend_lis, loc='best').get_frame().set_linewidth(0.0)
    plt.show()


def picture_subplot(title, xlabel, ylabel, data, cols):
    ylabel_li = ylabel.split(',')
    ylabel_use = []
    if len(ylabel_li) != len(cols):
        ylabel_use = [ylabel for i in range(len(cols))]
    else:
        ylabel_use = ylabel_li

    plot_num = cols_num_gen(cols)
    print(len(data[0]))
    for i in cols:
        ax = plt.subplot(len(cols), 1, next(plot_num))
        ax.plot(data[0][1:], data[i][1:])
        if i == cols[-1]:
            plt.xlabel(xlabel)
        ax.set_ylabel(ylabel_use[i-1])
        if i != cols[-1]:
            ax.set_xticks([])
        legend_lis = []
        legend_lis.append(data[i][0])
        ax.legend(labels=legend_lis, loc='best').get_frame().set_linewidth(0.0)
        if i == cols[0]:
            plt.title(title)
    plt.show()


def xvg_deal(filename):
    with open(filename, 'r') as fo:
        content = fo.read()
    line_list = content.strip('\n').split('\n')
    column_num = len(line_list[-1].split())
    data = []
    for i in range(column_num):
        data.append([])
    data[0].append('time')

    for line in line_list:
        if line[0] == "@":
            if "title" in line:
                title_line = line.strip('"')
                if '"' in title_line:
                    title = title_line.split('"')[-1]
                else:
                    title = 'Null'
            elif "xaxis" in line:
                xlabel_line = line.strip('"')
                if '"' in xlabel_line:
                    xlabel = xlabel_line.split('"')[-1]
                else:
                    xlabel = 'Null'
            elif "yaxis" in line:
                ylabel_line = line.strip('"')
                if '"' in ylabel_line:
                    ylabel = ylabel_line.split('"')[-1]
                else:
                    ylabel = 'Null'
            elif "legend" in line and "@ s" in line:
                legend = line.strip('"').split('"')[-1]
                for i in range(1, column_num):
                    if len(data[i]) == 0:
                        data[i].append(legend)
                        break

        elif line[0] != '#' and line[0] != '@':
            num_list = line.strip().split()
            for i in range(column_num):
                if len(data[i]) == 0:
                    data[i].append('no-legend')
                data[i].append(float(num_list[i]))

    return title, xlabel, ylabel, data


def energy_compute():
    try:
        pro_lig_file = sys.argv[1]
        pro_file = sys.argv[2]
        lig_file = sys.argv[3]
    except:
        print("No input filename! ")
        return 
    else:
        filename_output = ""
        plotmode = '-s'
        try:
            if sys.argv[4] == '-o' or sys.argv[4] == '-s':
                plotmode = sys.argv[4]
            elif '-f' in sys.argv[4]:
                filename_output = sys.argv[4][2:]
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
        
        if len(pro_lig_data) != 5 or len(pro_data) != 5 or len(lig_data) != 5:
            print( "Wrong, data columns not equal to 5! ")
            print("Each file must contain 5 columns: ")
            print("    | LJ-SR | Disper.corr. | Coulomb-SR | Coul.-recip. |" )
            return 

        final_data = pro_lig_data
        for i in range(1, len( final_data)):
            for j in range(1, len(final_data[i])):
                final_data[i][j] -= pro_data[i][j] + lig_data[i][j]
                # print(final_data[i][j], pro_data[i][j], lig_data[i][j])

        # ETOTAL
        final_data.append([0 for x in range(len(final_data[0]))])
        final_data[-1][0] = "ETOTAL"
        for j in range(1, len(final_data[0])):
            for i in range(1, len(final_data)-1):
                final_data[-1][j] += final_data[i][j]

        # COULOMB
        final_data.append([0 for x in range(len(final_data[0]))])
        final_data[-1][0] = "COULOMB"
        for j in range(1, len(final_data[0])):
            for i in [ 3, 4 ]: 
                final_data[-1][j] += final_data[i][j]

        with open("energy_results_" + filename_output + ".xvg", 'w') as fo:
            fo.write("## energy_results_" + filename_output + ".xvg generated from ")
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
