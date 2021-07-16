# author : Charlie
# date : 20210716

import os
import sys
import statistics

def MovingAverage_test(data):
    MA_data = []
    for column in data:
        MA_column = []
        MA_column.append(statistics.mean(column[:1]))
        MA_column.append(statistics.mean(column[:3]))
        for i in range(2, len(column)-2):
            MA_column.append(statistics.mean(column[i-2:i+3]))
        MA_column.append(statistics.mean(column[-3:]))
        MA_column.append(statistics.mean(column[-1:]))
        MA_data.append(MA_column)

    return MA_data


def MovingAverage(data, MA):
    MA_data = []
    if MA % 2 != 1:
        MA -= 1
        print("> MA must be odd, using ", MA)
    hMA = int((MA-1)/2)
    for column in data:
        MA_column = []
        for j in range(1,MA,2):
            MA_column.append(statistics.mean(column[:j]))
        for i in range(hMA, len(column)-hMA):
            MA_column.append(statistics.mean(column[i-hMA:i+hMA+1]))
        for j in range(-1*MA+2,0,2):
            MA_column.append(statistics.mean(column[j:]))
        MA_data.append(MA_column)

    return MA_data


def main():
    print("\n==Usage : xvg_movingaverage.py inputfile outputfile MA ==")
    print("==   eg : xvg_movingaverage.py rmsd.xvg  rmsd_MA.xvg 5 ==\n")
    inputfile, outputfile = "", ""
    MA_num = 5
    if len(sys.argv) == 2:
        inputfile = sys.argv[1]
        outputfile = inputfile.split(".")[0] + "_out.xvg"
    elif len(sys.argv) == 3:
        inputfile = sys.argv[1]
        outputfile = sys.argv[2]
    elif len(sys.argv) == 4:
        inputfile = sys.argv[1]
        outputfile = sys.argv[2]
        MA_num = int(sys.argv[3])
    elif len(sys.argv) > 4:
        print("> ERROR, too many input arguments !")

    with open(inputfile, "r") as fo:
        content = fo.read()
    lines = content.strip("\n").split("\n")
    comments = []
    if len(lines[-2].strip().split()) >= 1 and len(lines[-2].strip().split()) == len(lines[-3].strip().split()):
        column_num = len(lines[-2].strip().split())
    print("Number of columns -> ", column_num)
    data = [[] for i in range(column_num)]
    for line in lines:
        if line[0] == '#' or line[0] == '@':
            comments.append(line)
        elif line[0] == ' ' and "time" in line and len(line.strip().split()) >5:
            comments.append(line)
        else:
            items = line.strip().split()
            if len(items) != column_num:
                print("> ERROR, len(line) != column_num")
                exit(0)
            for i in range(column_num):
                data[i].append(float(items[i]))
    # moving average
    data = MovingAverage(data, MA_num)
    # write outputfile
    output = "\n".join(comments)
    for j in range(len(data[0])):
        output += "\n" + " ".join(["{:>17.4f}".format(data[i][j]) for i in range(column_num)])
    
    with open(outputfile, 'w') as fo :
        fo.write(output)
    print("> Done !")



if __name__ == "__main__":
    main()
