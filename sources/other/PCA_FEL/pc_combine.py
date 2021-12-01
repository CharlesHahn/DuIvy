# author : charlie
# date : 20210904

import os
import sys

def main():
    file1, file2, outfile = "", "", ""
    cmd = sys.argv[1:]
    if len(cmd) == 0 : 
        print(" Usage : pc_combine.py pc1 pc2 output-file ")
        exit()
    elif len(cmd) == 3:
        file1 = cmd[0]
        file2 = cmd[1]
        outfile = cmd[2]
    else: 
        print("wrong input, check it")
        exit()
    # deal with file 
    with open(file1, 'r') as fo:
        content1 = fo.read()
    with open(file2, 'r') as fo:
        content2 = fo.read()
    lines1 = content1.strip().strip("\n").split("\n")
    lines2 = content2.strip().strip("\n").split("\n")
    data1 = [li.strip().split() for li in lines1 if li[0] != '@' and li[0] != '&' and li[0] != '#']
    data2 = [li.strip().split() for li in lines2 if li[0] != '@' and li[0] != '&' and li[0] != '#']
    # print(data1)
    # print(data2)
    output = ""
    if len(data1) != len(data2):
        print("wrong length of data1 and data2, check it ")
        exit()
    for i in range(len(data1)):
        if float(data1[i][0]) == float(data2[i][0]):
            output += "{:16} {:16} {:16} \n".format(data1[i][0], data1[i][1], data2[i][1])
    with open(outfile, "w") as fo:
        fo.write(output)
    print("done")



if __name__ == "__main__":
    main()
