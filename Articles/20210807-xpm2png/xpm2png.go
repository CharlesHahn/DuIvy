// author : charlie
// date : 20210805

package main

import (
	"fmt"
	"image"
	"image/color"
	"image/png"
	"io/ioutil"
	"os"
	"strconv"
	"strings"
)

func hex2rgb(hex_color string) []int64 {
	rgb := make([]int64, 0)
	r, _ := strconv.ParseInt(hex_color[1:3], 16, 0)
	rgb = append(rgb, r)
	g, _ := strconv.ParseInt(hex_color[3:5], 16, 0)
	rgb = append(rgb, g)
	b, _ := strconv.ParseInt(hex_color[5:7], 16, 0)
	rgb = append(rgb, b)
	return rgb
}

func parse_save(inputfile, outputfile string) {
	// parse xpm data
	xpm_title := ""
	xpm_legend := ""
	xpm_type := ""
	xpm_xlabel := ""
	xpm_ylabel := ""
	var xpm_width int64
	var xpm_height int64
	var xpm_color_num int64
	var xpm_char_per_pixel int64
	char2color := make(map[string]string)
	char2note := make(map[string]string)
	xpm_xaxis := make([]string, 0)
	xpm_yaxis := make([]string, 0)
	xpm_data := make([]string, 0)

	buf, err := ioutil.ReadFile(inputfile)
	if err != nil {
		fmt.Println(err.Error())
		fmt.Println("ERROR -> error in reading ", inputfile)
		return
	}
	content := string(buf)
	str_buf := strings.TrimSpace(content)
	str_buf = strings.Trim(str_buf, "\n")
	lines := strings.Split(str_buf, "\n")

	flag_4code := 0 // means have not detected yet
	for _, line := range lines {
		if flag_4code == 1 { // means this line is 4code line
			flag_4code = 2
			// this is quite strange, only TrimSpace works
			line = strings.Trim(strings.TrimSpace(line), ",")
			//line = strings.ReplaceAll(line, ",", "")
			items_withblank := strings.Split(strings.Trim(line, "\""), " ")
			items := make([]string, 0)
			for _, item := range items_withblank {
				if item != "" {
					items = append(items, item)
				}
			}
			xpm_width, _ = strconv.ParseInt(items[0], 10, 0)
			xpm_height, _ = strconv.ParseInt(items[1], 10, 0)
			xpm_color_num, _ = strconv.ParseInt(items[2], 10, 0)
			xpm_char_per_pixel, _ = strconv.ParseInt(items[3], 10, 0)
			continue
		} else if flag_4code == 0 && strings.HasPrefix(line, "static char") {
			flag_4code = 1 // means next line is 4code line
			continue
		}

		if strings.HasPrefix(line, "/* x-axis") {
			items := strings.Split(strings.TrimSpace(line), " ")
			for _, item := range items[2 : len(items)-1] {
				if item != "" {
					xpm_xaxis = append(xpm_xaxis, item)
				}
			}
			continue
		} else if strings.HasPrefix(line, "/* y-axis") {
			items := strings.Split(strings.TrimSpace(line), " ")
			for _, item := range items[2 : len(items)-1] {
				if item != "" {
					xpm_yaxis = append(xpm_yaxis, item)
				}
			}
			continue
		} else if strings.HasPrefix(line, "/* title") {
			xpm_title = strings.Split(strings.TrimSpace(line), "\"")[1]
			continue
		} else if strings.HasPrefix(line, "/* legend") {
			xpm_legend = strings.Split(strings.TrimSpace(line), "\"")[1]
			continue
		} else if strings.HasPrefix(line, "/* x-label") {
			xpm_xlabel = strings.Split(strings.TrimSpace(line), "\"")[1]
			continue
		} else if strings.HasPrefix(line, "/* y-label") {
			xpm_ylabel = strings.Split(strings.TrimSpace(line), "\"")[1]
			continue
		} else if strings.HasPrefix(line, "/* type") {
			xpm_type = strings.Split(strings.TrimSpace(line), "\"")[1]
			continue
		}

		items_withblank := strings.Split(strings.TrimSpace(line), " ")
		items := make([]string, 0)
		for _, item := range items_withblank {
			if item != "" {
				items = append(items, item)
			}
		}
		if len(items) == 7 && items[1] == "c" {
			char2color[strings.Trim(items[0], "\"")] = items[2]
			char2note[strings.Trim(items[0], "\"")] = strings.Trim(items[5], "\"")
		}
		if len(items) == 1 {
			xpm_data = append(xpm_data, strings.Trim(strings.Trim(strings.TrimSpace(line), ", "), "\""))
		}
	}

	// data check and transformation
	//fmt.Println(xpm_title, xpm_legend, xpm_xlabel, xpm_ylabel, xpm_type)
	fmt.Println("Info -> title : ", xpm_title)
	fmt.Println("Info -> legend : ", xpm_legend)
	fmt.Println("Info -> xpm type : ", xpm_type)
	fmt.Println("Info -> x-label : ", xpm_xlabel)
	fmt.Println("Info -> y-label : ", xpm_ylabel)
	fmt.Println("Info -> 4 code line : ", xpm_width, xpm_height, xpm_color_num, xpm_char_per_pixel)
	fmt.Println("Info -> see x-axis and y-axis in your input xpm file")
	//fmt.Println(len(xpm_xaxis))
	//fmt.Println(len(xpm_yaxis))
	//fmt.Println(char2color)
	//fmt.Println(char2note)
	//fmt.Println(len(xpm_data))
	char2rgb := make(map[string][]int64)
	for key, value := range char2color {
		char2rgb[key] = hex2rgb(value)
	}
	//fmt.Println(char2rgb)
	if len(xpm_xaxis) != int(xpm_width) {
		fmt.Println("Warning -> length of x-axis is not equal to xpm width !")
		fmt.Println("    length of xpm_xaxis : ", len(xpm_xaxis))
		fmt.Println("    xpm_width : ", xpm_width)
	}
	if len(xpm_yaxis) != int(xpm_height) {
		fmt.Println("Warning -> length of y-axis is not equal to xpm height !")
		fmt.Println("    length of xpm_yaxis : ", len(xpm_yaxis))
		fmt.Println("    xpm_heigth : ", xpm_height)
	}
	if len(char2color) != int(xpm_color_num) {
		fmt.Println("ERROR -> number of char is not equal to number of color, check your xpm file !")
		fmt.Println("    char2color : ", char2color)
		return
	}
	if len(char2note) != int(xpm_color_num) {
		fmt.Println("ERROR -> number of char is not equal to number of note, check your xpm file !")
		fmt.Println("    char2note : ", char2note)
		return
	}
	if len(xpm_data) != int(xpm_height) {
		fmt.Println("ERROR -> raws of data is not equal to xpm height, check it !")
		fmt.Println("    raws of data : ", len(xpm_data))
		fmt.Println("    xpm_height :", xpm_height)
		return
	}
	for i := 0; i < len(xpm_data); i++ {
		if len(xpm_data[i]) != int(xpm_width*xpm_char_per_pixel) {
			fmt.Println("ERROR -> at line ", i+1, " of xpm data, length of char is not equal to xpm width, check it !")
			fmt.Println(xpm_data[i])
			fmt.Println("    length of char : ", len(xpm_data[i]))
			fmt.Println("    xpm_width*xpm_char_per_pixel : ", xpm_width*xpm_char_per_pixel)
			return
		}
	}

	// new a picture
	file, err := os.Create(outputfile)
	if err != nil {
		fmt.Println(err.Error())
		fmt.Println("ERROR -> error in creating ", outputfile)
		return
	}
	defer file.Close()
	figure := image.NewRGBA(image.Rect(0, 0, int(xpm_width), int(xpm_height)))
	for x := 0; x < int(xpm_width); x++ {
		for y := 0; y < int(xpm_height); y++ {
			char := xpm_data[y][x*int(xpm_char_per_pixel) : x*(int(xpm_char_per_pixel))+int(xpm_char_per_pixel)]
			rgb := char2rgb[char]
			figure.Set(x, y, color.RGBA{
				uint8(rgb[0]), uint8(rgb[1]), uint8(rgb[2]), 255})
		}
	}
	err = png.Encode(file, figure)
	if err != nil {
		fmt.Println(err)
	}
	fmt.Println("Info -> saving figure to ", outputfile, " successfully")
	fmt.Println("Info -> All done, enjoy your day !")

}

func main() {
	Usage := "xpm2png.go : turn xpm file generated by GROMACS into png file\n"
	Usage += "        -f : <.xpm> xpm file to visualization\n"
	Usage += "        -o : <.png> png file to save figure, (xpm2png.png)(default)\n"
	Usage += "        -h : show help message\n"
	Usage += "eg. xpm2png -f input.xpm -o xpm2png.png\n"

	cmd_list := os.Args[1:]
	if len(cmd_list) == 0 {
		fmt.Print(Usage)
		return
	} else if len(cmd_list) > 4 {
		fmt.Println("ERROR -> too many arguments, check your command")
		return
	}
	inputfile := ""
	outputfile := "xpm2png.png"
	flag_input := 0
	flag_output := 0
	for _, cmd := range cmd_list {
		if flag_input == 1 {
			flag_input = 0
			if len(cmd) < 4 || cmd[len(cmd)-4:] != ".xpm" {
				fmt.Println("ERROR -> make sure your input file is XPM file with suffix .xpm")
				return
			}
			inputfile = cmd
			continue
		} else if flag_output == 1 {
			flag_output = 0
			if len(cmd) < 4 || cmd[len(cmd)-4:] != ".png" {
				fmt.Println("ERROR -> make sure your output file is PNG file with suffix .png")
				return
			}
			outputfile = cmd
			continue
		}
		if cmd[0] == '-' {
			if cmd == "-f" {
				flag_input = 1
			} else if cmd == "-o" {
				flag_output = 1
			} else if cmd == "-h" {
				fmt.Print(Usage)
				return
			} else {
				fmt.Println("ERROR -> unknown command argument ", cmd, ", check it !")
				return
			}
		} else {
			fmt.Println("ERROR -> unknown argument ", cmd, ", check it !")
			return
		}
	}
	fmt.Println("Info -> input: ", inputfile, ", output: ", outputfile)
	// parse xpm and visualization
	parse_save(inputfile, outputfile)

}
