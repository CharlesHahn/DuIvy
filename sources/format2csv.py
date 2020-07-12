# author : charlie
###
# for formatting the energy_compute_results.xvg
###

import sys


def write_standard(line, file_output):
    with open(file_output + "_formatted.xvg", 'a', encoding='utf-8') as fo:
        line = line.strip(",")
        fo.write( line + '\n')


def format_style(file_input):
    file_output = file_input.split('.')[0]
    with open(file_input, 'r', encoding='utf-8') as fo:
        content = fo.read()
    lines = content.split('\n')
    # write_standard("## standard -> " + lines[0], file_output)
    title_origin = lines[1].split()
    title_deal = []
    for title_o in title_origin:
        if '(' == title_o[0]:
            title_deal[-1] += title_o.strip()
        elif '.' == title_o.strip()[-1] and '.' == title_deal[-1][-1]:
            # print(title_o)
            title_deal[-1] += title_o.strip()
            # print(title_deal)
        else:
            title_deal.append(title_o)
    title_line = ""
    for title in title_deal:
        title_line += "{:>20},".format(title)
    write_standard(title_line, file_output)
    for line in lines[2:]:
        line_lis = line.split()
        write_content = ""
        for li in line_lis:
            write_content += "{:>20.2f},".format(float(li))
        write_standard(write_content, file_output)
    
    print(file_output + " done ~")


def main():
    cmds = [i for i in sys.argv]
    for cmd in cmds[1:]:
        format_style(cmd)
    print("ALL DONE ! ")

if __name__ == "__main__":
    main()
