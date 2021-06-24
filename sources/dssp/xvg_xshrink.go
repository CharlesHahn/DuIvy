// Author : Charlie
// Date : 20210505

package main

import (
	"fmt"
	"io/ioutil"
	"os"
	"strconv"
	"strings"
)

func main() {
	inputfile := ""
	shrink := 1000
	outputfile := ""
	// os command line args
	if len(os.Args) == 1 {
		fmt.Println("No input ! Error !")
		fmt.Println("type xvg_shrink.exe -h for help info.")
		return
	} else if len(os.Args) == 2 {
		inputfile = os.Args[1]
		if inputfile == "-h" || inputfile == "help" || inputfile == "--help" {
			fmt.Println("== Usage: xvg_xshrink.exe <InputFileName> <ShrinkTimes> <OutputFileName> == ")
			fmt.Println("       eg. xvg_shrink.exe input.xvg 1000 output.xvg")
			return
		}
	} else if len(os.Args) == 3 {
		inputfile = os.Args[1]
		shrink_times, err := strconv.ParseInt(os.Args[2], 10, 64)
		if err != nil {
			fmt.Print("Error to parse ", os.Args[2], " to shrink times, use it as output filename for instead!")
			outputfile = os.Args[2]
		} else {
			shrink = int(shrink_times)
		}
	} else if len(os.Args) == 4 {
		inputfile = os.Args[1]
		shrink_times, err := strconv.ParseInt(os.Args[2], 10, 64)
		if err != nil {
			fmt.Print("Unable to parse ", os.Args[2], " to shrink times. It should be Integer !")
			fmt.Println("type xvg_shrink.exe -h for help info.")
			return
		} else {
			shrink = int(shrink_times)
		}
		outputfile = os.Args[3]
	} else if len(os.Args) > 4 {
		fmt.Println("Too many input, only the 1st was adopted !")
		fmt.Println("type xvg_shrink.exe -h for help info.")
	}

	output := ""
	// read file
	buf, err := ioutil.ReadFile(inputfile)
	if err != nil {
		panic(err.Error())
	}
	// 去掉最后可能的空格和空行
	str_buf := strings.Trim(string(buf), " ")
	str_buf = strings.Trim(str_buf, "\n")
	lines := strings.Split(str_buf, "\n")
	fmt.Println("Lines -> ", len(lines))
	// 切分会多出空元素
	//column_number := len(strings.Split(lines[len(lines)-1], " "))
	column_number := len(strings.Fields(lines[len(lines)-1]))
	column_number_1 := len(strings.Fields(lines[len(lines)-2]))
	if column_number != column_number_1 {
		fmt.Println("The columns numbers in different rows are not the same !")
		fmt.Println("-1 line-> ", column_number, "; -2 line-> ", column_number_1)
		fmt.Println("You may need to delete blank lines in the tail of file ！")
		return
	}
	fmt.Println("column_number -> ", column_number)

	// 读取数据
	for _, line := range lines {
		if line[0] == '@' || line[0] == '#' {
			output = output + line + "\n"
		}
		if line[0] != '@' && line[0] != '#' {
			for index, item := range strings.Fields(line) {
				float, err := strconv.ParseFloat(item, 62)
				//fmt.Println(float, err)
				if err != nil {
					fmt.Println("wrong in parse data string to float !")
					fmt.Println("some figure in your file may be not correct number, check it !")
				}
				if index == 0 {
					float = float / float64(shrink)
				}
				//output += strconv.FormatFloat(float, 'f', 6, 64) + "  "
				output += fmt.Sprintf("%18.6f", float)
			}
			output += "\n"
		}
	}

	fileName := strings.Split(inputfile, ".")[0] + "_xshrink1000.xvg"
	if outputfile != "" {
		fileName = outputfile
	}
	dstFile, err := os.Create(fileName)
	if err != nil {
		fmt.Println(err.Error())
		return
	}
	defer dstFile.Close()
	dstFile.WriteString(output)

	fmt.Println("Content has been written into " + fileName)
	fmt.Println("Successfully finished !")

}
