#############################################################################
# author : charlie
# time : 20200714
# command : python xvg_average.py xvgfilename full 10 100
# usage : 
#    xvg_average.py XVG_Filename column_select start_index end_index
#        column_select -> e.g. 1,3 or full; default full
#        start_index   -> optional, e.g. 10; default 0
#        end_index     -> optional, but must be POSITIVE, NO -1 ! default End
#############################################################################

import os
import sys
from statistics import mean

def loadxvg(file):
    with open(file, 'r') as fo:
        content = fo.read()
    lines = content.strip().split('\n')
    column_num = len([i for i in lines[-1].split() if i !=' '])
    print('>> ' + file + " ->>> column_num == " + str(column_num) )
    data = [ [] for i in range(column_num)]
    for line in lines:
        if '@' in line:
            print(' ---> ' + line)
            continue
        elif line[0] == '#':
            continue
        line_clean = [i for i in line.split() if i != ' ']
        for i in range(column_num):
            data[i].append(line_clean[i])
    
    return data


def main():
    help_str = """
== command : 
    xvg_average.py XVG_Filename column_select start_index end_index
        column_select -> e.g. 1,3 or full; default full
        start_index   -> optional, e.g. 10; default 0
        end_index     -> optional, but must be POSITIVE, NO -1 ! default End
== \n"""
    cmds = [ i for i in sys.argv[1:]]

    if len(cmds) == 0:
        print(help_str)
        return
    elif len(cmds) == 1:
        filename = cmds[0]
        start_index = 0
        end_index = -1
        column_select = 'full'
    elif len(cmds) == 2:
        filename = cmds[0]
        start_index = 0
        end_index = -1
        column_select = cmds[1]
    elif len(cmds) == 3:
        filename = cmds[0]
        start_index = cmds[2]
        end_index = -1
        column_select = cmds[1]
    elif len(cmds) == 4:
        filename = cmds[0]
        start_index = cmds[2]
        end_index = cmds[3]
        column_select = cmds[1]
    else:
        print("\n --> too many arguments ! ")
        print(help_str)
        return 

    if filename not in os.listdir() or '.xvg' not in filename:
        print("\n --> wrong xvg filename ! ")
        print(help_str)
        return
    try:
        start_index = int(start_index)
        end_index = int(end_index)
    except:
        print("\n** Wrong format of start_index or end_index ! ")
        print(help_str)
        return

    data = loadxvg( filename )
    # print("data first line -> ")
    data_first_line = ''
    for i in range(len(data)):
    	data_first_line += "{: ^16}".format(data[i][0])
    # print(data_first_line)
    
    raw_max = len(data[0]) - 1
    if start_index > raw_max:
        start_index = 0 
        print("\n* start_index larger than raw_max " 
                + str(raw_max) + ", set it to be 0")
    if end_index > raw_max:
        end_index = -1
        print("\n** end_index larger than raw_max " 
                + str(raw_max) + ", set it to be -1")

    column_range = [ i for i in range(len(data))]
    column_show = []
    if column_select == 'full':
        column_show = column_range
    else:
        try:
            column_select_deal = column_select.split(',')
            column_show = [ int(i) for i in column_select_deal
                    if int(i) in column_range ]
        except:
            column_show = column_range
            print('\n** wrong column_select, set it to be full ! ')
        
    # print(start_index, end_index, column_show)
    table_title = ""
    table_content = ""
    table_seprate = ""
    table_line = ""
    table_first_line = ''
    for i in column_show:
        table_first_line += "{: ^16}".format(data[i][0])
        # if original xvg file, add first line into compute
        try:
            float(data[i][0])
        except:
            start_index = 1
            print("  ->> Unable to convert the raw 0 '{}' into float,".format(
                data[i][0]), end = '')
            print("start_index set to be 1")
        else:
            start_index = 0
        table_seprate += "{:=^16}".format("=")
        table_line += "{:-^16}".format('-')
        table_title += "{:^16}".format(i)
        if end_index == -1:
            table_content += "{:^16.4f}".format(
                mean([float(num) for num in 
                    data[i][start_index : ]]))

        else:
            table_content += "{:^16.4f}".format(
                mean([float(num) for num in 
                    data[i][start_index : end_index +1 ]]))

    print("=========" + table_seprate)
    print("title    " + table_title)
    print("firstLine" + table_first_line)
    print("---------" + table_line)
    print("average  " + table_content)
    print("=========" + table_seprate)



if __name__ == "__main__":
    main()
