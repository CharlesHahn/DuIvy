// Author : Charlie
// Date : 20210520

package main

import (
	"fmt"
	"io/ioutil"
	"os"
	"strconv"
	"strings"
)

func main() {
	fmt.Println("== Usage : xpm_xshrink1000.exe filename.xpm ==")
	inputfile := ""
	// os command line args
	if len(os.Args) == 1 {
		fmt.Println("No input ! Error !")
		return
	} else if len(os.Args) == 2 {
		inputfile = os.Args[1]
	} else if len(os.Args) > 2 {
		inputfile = os.Args[1]
		fmt.Println("Too many input, only the 1st was adopted !")
	}
	// read file
	buf, err := ioutil.ReadFile(inputfile)
	if err != nil {
		panic(err.Error())
	}
	// 去掉最后可能的空格和空行
	str_buf := strings.Trim(string(buf), " ")
	str_buf = strings.Trim(str_buf, "\n")
	lines := strings.Split(str_buf, "\n")
	output := make([]string, 0)
	//fmt.Println("Lines -> ", len(lines))
	for _, line := range lines {
		if len(line) > 10 && line[0:10] == "/* x-axis:" {
			//fmt.Println(line)
			items := strings.Split(line, " ")
			//fmt.Println(items)
			line_modified := ""
			for ind, it := range items {
				its := strings.Trim(it, " ")
				if ind == 0 || ind == 1 || ind == len(items)-1 || len(its) == 0 {
					line_modified = line_modified + " " + its
				} else {
					number_int, err := strconv.ParseInt(its, 10, 64)
					if err != nil {
						fmt.Println(err.Error)
						fmt.Println("Wrong to read x-axis number, be sure it's interger ")
						fmt.Println("Unsuccessful !")
						return
					}
					number_float := float64(number_int) / 1000.0
					//fmt.Println(its)
					str_float := strconv.FormatFloat(float64(number_float), 'f', -1, 64)
					line_modified = line_modified + " " + str_float
				}
			}
			line_modified = strings.Trim(line_modified, " ")
			//fmt.Println(line_modified)
			output = append(output, line_modified)
		} else if len(line) > 11 && line[0:1] == "/" && line[0:11] == "/* x-label:" {
			output = append(output, "/* x-label: \"Time(ns)\" */")
		} else {
			output = append(output, line)
		}
	}
	//fmt.Print(len(output))

	fileName := strings.Split(inputfile, ".")[0] + "_xshrink1000.xpm"
	dstFile, err := os.Create(fileName)
	if err != nil {
		fmt.Println(err.Error())
		return
	}
	defer dstFile.Close()
	out2file := strings.Join(output, "\n")
	dstFile.WriteString(out2file + "\n")

	fmt.Println("Content has been written into " + fileName)
	fmt.Println("Successfully finished !")
}
