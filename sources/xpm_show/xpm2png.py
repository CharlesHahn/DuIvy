# author : charlie
# data : 20210803

# /* <width/columns> <height/rows> <colors> <chars per pixel>*/

import sys
import matplotlib.pyplot as plt 


# hex color to RGB color
def hex2rgb(hex_color):
    r = int(hex_color[1:3],16)
    g = int(hex_color[3:5],16)
    b = int(hex_color[5:7], 16)
    rgb = [r, g, b]
    return rgb


def parse_show(inputfile, show, outputfile, ip):
    # parse xpm data
    xpm_title, xpm_legend, xpm_type = "", "", ""
    xpm_xlabel, xpm_ylabel = "", ""
    xpm_width, xpm_height = 0, 0 
    xpm_color_num, xpm_char_per_pixel = 0, 0
    char2color, char2note = {}, {}
    xpm_xaxis, xpm_yaxis, xpm_data = [], [], []

    # read data
    with open(inputfile, 'r') as fo:
        lines = fo.read().strip().strip("\n").split("\n")

    flag_4_code = 0 # means haven't detected yet
    for line in lines:
        if flag_4_code == 1: # means this line is code4 line
            flag_4_code = 2 # means have detected
            code4 = [int(c) for c in line.strip().strip(",").strip("\"").split()]
            xpm_width, xpm_height = code4[0], code4[1]
            xpm_color_num, xpm_char_per_pixel = code4[2], code4[3]
            continue
        elif (flag_4_code == 0) and line.startswith("static char"):
            flag_4_code = 1 # means next line is code4 line
            continue

        if line.startswith("/* x-axis"):
            xpm_xaxis += [float(n) for n in line.strip().split()[2:-1] ]
            continue
        elif line.startswith("/* y-axis"):
            xpm_yaxis += [float(n) for n in line.strip().split()[2:-1] ]
            continue
        elif line.startswith("/* title"):
            xpm_title = line.strip().split("\"")[1]
            continue
        elif line.startswith("/* legend"):
            xpm_legend = line.strip().split("\"")[1]
            continue
        elif line.startswith("/* x-label"):
            xpm_xlabel = line.strip().split("\"")[1]
            continue
        elif line.startswith("/* y-label"):
            xpm_ylabel = line.strip().split("\"")[1]
            continue
        elif line.startswith("/* type"):
            xpm_type = line.strip().split("\"")[1]
            continue

        items = line.strip().split()
        if len(items) == 7 and items[1] == "c":
            if len(items[0].strip("\"")) == xpm_char_per_pixel:
                char2color[items[0].strip("\"")] = items[2]
                char2note[items[0].strip("\"")] = items[5].strip("\"")
            if items[0].strip() == '"':
                char2color[" "] = items[2]
                char2note[" "] = items[5].strip("\"")
        if line.strip().startswith('"') and len(line.strip().strip(",").strip("\"")) == xpm_width*xpm_char_per_pixel :
            xpm_data.append(line.strip().strip(",").strip("\""))

    ## data check and transformation
    # print(xpm_title, xpm_legend, xpm_xlabel, xpm_ylabel, xpm_type)
    # print(xpm_width, xpm_height, xpm_color_num, xpm_char_per_pixel)
    # print(len(xpm_xaxis))
    # print(len(xpm_yaxis))
    # print(char2color)
    # print(char2note)
    # print(len(xpm_data))
    char2rgb = {}
    for key, value in char2color.items():
        char2rgb[key] = hex2rgb(value)
    if len(xpm_xaxis) != xpm_width:
        print("Warning -> length of x-axis is not equal to xpm width !")
        if len(xpm_xaxis) == xpm_width + 1:
            xpm_xaxis = xpm_xaxis[:-1]
            if len(xpm_xaxis) != xpm_width:
                print("ERROR -> abandoned the last value of x-axis and it still failed.")
                exit()
            print("Warning -> the last value of x-axis has been abandoned.")
        else: 
            exit()
    if len(xpm_yaxis) != xpm_height:
        print("Warning -> length of y-axis is not equal to xpm height !")
        if len(xpm_yaxis) == xpm_width + 1:
            xpm_yaxis = xpm_yaxis[:-1]
            if len(xpm_yaxis) != xpm_width:
                print("ERROR -> abandoned the last value of y-axis and it still failed.")
                exit()
            print("Warning -> the last value of y-axis has been abandoned.")
        else: 
            exit()
    if len(char2color) != xpm_color_num:
        print("ERROR -> length of char2color is not equal to xpm_color_num, check it !")
        exit()
    if len(char2rgb) != xpm_color_num:
        print("ERROR -> length of char2rgb is not equal to xpm_color_num, ", end="")
        print("wrong in change hex color to rgb color.")
        exit()
    if len(char2note) != xpm_color_num:
        print("ERROR -> length of char2note is not equal to xpm_color_num, check it !")
        exit()
    if len(xpm_data) != xpm_height:
        print("ERROR -> raws of data is not equal to xpm height, check it !")
        exit()
    for i in range(len(xpm_data)):
        if len(xpm_data[i]) != xpm_width*xpm_char_per_pixel:
            print("ERROR -> at line ", i+1, " of xpm data, length of char is ", end="")
            print("not equal to xpm width, check it !")
            exit()
    # reverse xpm_yaxis
    xpm_yaxis.reverse()
    # deal with axis values
    flag_xaxis = 0
    for x in xpm_xaxis:
        if x % 1 != 0:
            flag_xaxis = 1
    if flag_xaxis == 0:
        xpm_xaxis = [int(x) for x in xpm_xaxis]
    elif flag_xaxis == 1:
        xpm_xaxis = [int(x*100)/100 for x in xpm_xaxis]
    flag_yaxis = 0
    for y in xpm_yaxis:
        if y % 1 != 0:
            flag_yaxis = 1
    if flag_yaxis == 0:
        xpm_yaxis = [int(y) for y in xpm_yaxis]
    elif flag_yaxis == 1:
        xpm_yaxis = [int(y*100)/100 for y in xpm_yaxis]


    print("Info -> all data has been read from xpm file and checked successfully.")

    # visualization
    img = []
    for line in xpm_data:
        rgb_line = []
        for i in range(0,xpm_width*xpm_char_per_pixel, xpm_char_per_pixel):
            rgb_line.append( char2rgb[ line[i : i+xpm_char_per_pixel] ] )
        img.append(rgb_line)

    plt.figure()
    plt.imshow(img)
    plt.title(xpm_title)
    plt.xlabel(xpm_xlabel)
    plt.ylabel(xpm_ylabel)
    print("Legend of this xpm figure -> ", xpm_legend)
    # set the ticks
    if xpm_width < 100:
        x_tick = int(xpm_width/3)
    elif xpm_width >= 100 and xpm_width < 1000:
        x_tick = int(xpm_width/5)
    elif xpm_width > 500:
        x_tick = int(xpm_width/10)
    if xpm_height < 100:
        y_tick = int(xpm_height/3)
    elif xpm_height >= 100 and xpm_height < 1000:
        y_tick = int(xpm_height/5)
    elif xpm_height > 500:
        y_tick = int(xpm_height/10)
    plt.tick_params(axis='both',which='major')
    x_major_locator = plt.MultipleLocator(x_tick)
    y_major_locator = plt.MultipleLocator(y_tick)
    ax=plt.gca()
    ax.xaxis.set_major_locator(x_major_locator)
    ax.yaxis.set_major_locator(y_major_locator)
    plt.xticks([w for w in range(1,xpm_width+1, x_tick)],
            [xpm_xaxis[i] for i in range(0,len(xpm_xaxis), x_tick)])
    plt.yticks([h for h in range(1,xpm_height+1, y_tick)], 
            [xpm_yaxis[i] for i in range(0,len(xpm_yaxis), y_tick)])
    
    # to deal with commands
    if outputfile != "":
        if not (xpm_type == "Continuous" and ip == "yes"):
            plt.savefig(outputfile, dpi=600)
            print("Info -> ",outputfile," has been saved in present working directory")
    if show == "yes":
        if xpm_type == "Discrete" and ip == "yes":
            print("Warning -> Can not apply interpolation to Discrete type xpm")
            plt.show()
        elif xpm_type == "Discrete" and ip == "no":
            plt.show()
        elif xpm_type == "Continuous" and ip == "no":
            plt.show()

    # visualization, interpolation
    if xpm_type == "Continuous" and ip == "yes":
        plt.clf()
        for key, value in char2note.items():
            char2note[key] = float(value)
        img = []
        for line in xpm_data:
            value_line = []
            for i in range(0,xpm_width*xpm_char_per_pixel, xpm_char_per_pixel):
                value_line.append(char2note[ line[i: i+xpm_char_per_pixel]])
            img.append(value_line)
        # you could modify cmap and interpolation method here
        im = plt.imshow(img, cmap='jet', interpolation="bilinear")
        plt.colorbar(im, fraction=0.046, pad = 0.04)
        plt.title(xpm_title)
        plt.xlabel(xpm_xlabel)
        plt.ylabel(xpm_ylabel)
        # set the ticks
        if xpm_width < 100:
            x_tick = int(xpm_width/3)
        elif xpm_width >= 100 and xpm_width < 1000:
            x_tick = int(xpm_width/5)
        elif xpm_width > 500:
            x_tick = int(xpm_width/10)
        if xpm_height < 100:
            y_tick = int(xpm_height/3)
        elif xpm_height >= 100 and xpm_height < 1000:
            y_tick = int(xpm_height/5)
        elif xpm_height > 500:
            y_tick = int(xpm_height/10)
        plt.tick_params(axis='both',which='major')
        x_major_locator = plt.MultipleLocator(x_tick)
        y_major_locator = plt.MultipleLocator(y_tick)
        ax=plt.gca()
        ax.xaxis.set_major_locator(x_major_locator)
        ax.yaxis.set_major_locator(y_major_locator)
        plt.xticks([w for w in range(1,xpm_width+1, x_tick)],
                [xpm_xaxis[i] for i in range(0,len(xpm_xaxis), x_tick)])
        plt.yticks([h for h in range(1,xpm_height+1, y_tick)],
                [xpm_yaxis[i] for i in range(0,len(xpm_yaxis), y_tick)])

        # to deal with commands
        if outputfile != "":
            plt.savefig(outputfile, dpi=600)
            print("Info -> ",outputfile," has been saved in present working directory")
        if show == "yes":
            plt.show()

    print("Info -> All done, wish you good day !")


def main():
    Usage = "xpm2png.py : python3 script for visualization of xpm file generated by GROMACS\n"
    Usage +="        -f : [<.xpm>]      (gibbs.xpm)\n"
    Usage +="               xpm file to be inputed\n"
    Usage +="     -show : [yes/no]      (yes)(default)\n"
    Usage +="               whether to show figure directly\n"
    Usage +="        -o : [<.png>]      (output.png)   (optional)\n"
    Usage +="               png file to save figure\n"
    Usage +="       -ip : [yes/no]      (no) (default)\n"
    Usage +="               whether apply interpolation to xpm data\n"
    Usage +="               ONLY effective for xpm file of Continuous type\n"
    Usage +="               interpolation would be useful for some cases, eg. FEL\n"
    Usage +="        -h :   show this usage info\n"

    # learned this skill from jerkwin, using nonpythonic style to save space
    # hate this long paragraph, must be for idiot user
    flag_input, flag_show, flag_output, flag_ip = 0, 0, 0, 0
    inputfile, show, outputfile, ip = "", "yes", "", "no"
    cmd_list = sys.argv[1:]
    if len(cmd_list) == 0:
        print(Usage)
        exit()
    if len(cmd_list)%2 != 0 and len(cmd_list) > 1:
        print("ERROR -> wrong command, missing arguments")
        print(Usage)
        exit()
    if len(cmd_list)%2 != 0 and len(cmd_list) == 1 and cmd_list[0] != "-h":
        print("ERROR -> wrong command, missing arguments")
        print(Usage)
        exit()
    if len(cmd_list)%2 != 0 and len(cmd_list) == 1 and cmd_list[0] == "-h":
        print(Usage)
        exit()
    for cmd in cmd_list:
        if flag_input == 1:
            flag_input = 0
            if cmd.startswith("-"):
                print("ERROR -> no argument for -f ! check your command")
                print("         or maybe you are using prefix '-', change it")
                print("         type 'xpm2png.py -h' for usage message")
                exit()
            if len(cmd) < 4 or cmd[-4:] != ".xpm":
                print("ERROR -> make sure your input file is XPM file with suffix .xpm")
                print("         type 'xpm2png.py -h' for usage message")
                exit()
            inputfile = cmd
            continue
        elif flag_show == 1:
            flag_show = 0
            if cmd.startswith("-"):
                print("ERROR -> no argument for -show ! check your command")
                print("         or maybe you are using prefix '-', change it")
                print("         type 'xpm2png.py -h' for usage message")
                exit()
            if cmd != "yes" and cmd != "no":
                print("ERROR -> wrong argument for -show, use yes or no !")
                print("         type 'xpm2png.py -h' for usage message")
                exit()
            show = cmd
            continue
        elif flag_output == 1:
            flag_output = 0
            if cmd.startswith("-"):
                print("ERROR -> no argument for -o ! check your command")
                print("         or maybe you are using prefix '-', change it")
                print("         type 'xpm2png.py -h' for usage message")
                exit()
            if len(cmd) < 4 or cmd[-4:] != ".png":
                print("ERROR -> make sure your output file is PNG file with suffix .png")
                print("         type 'xpm2png.py -h' for usage message")
                exit()
            outputfile = cmd
            continue
        elif flag_ip == 1:
            flag_ip = 0 
            if cmd.startswith("-"):
                print("ERROR -> no argument for -ip ! check your command")
                print("         or maybe you are using prefix '-', change it")
                print("         type 'xpm2png.py -h' for usage message")
                exit()
            if cmd != "yes" and cmd != "no":
                print("ERROR -> wrong argument for -ip, use yes or no !")
                print("         type 'xpm2png.py -h' for usage message")
                exit()
            ip = cmd
            continue
        if cmd.startswith("-"):
            if cmd == "-f":
                flag_input = 1; continue
            elif cmd == "-show":
                flag_show = 1; continue
            elif cmd == "-o":
                flag_output = 1; continue
            elif cmd == "-ip":
                flag_ip = 1; continue
            elif cmd == "-h":
                print(Usage); exit()
            else:
                print("ERROR -> wrong command : ",cmd," ! no this option found,")
                print("         type 'xpm2png.py -h' for usage message")
                exit()
        else:
            print("ERROR -> wrong command : ",cmd," ! no this option found,")
            print("         type 'xpm2png.py -h' for usage message")
            exit()

    # parse xpm and visualization
    parse_show(inputfile, show, outputfile, ip)


if __name__ == "__main__":
    main()

