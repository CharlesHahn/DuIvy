// Author : Charlie
// Date : 20210521

package main

import (
	"fmt"
	"io/ioutil"
	"os"
	"strconv"
	"strings"
)

func main() {
	fmt.Println(">>>")
	// initial the input args
	inputfile := ""
	colormap := ""
	output_png := ""

	// initial the colormap_list
	colormap_list := make([]string, 0)
	colormap_list = append(colormap_list, "xpm")
	colormap_list = append(colormap_list, "gmx")
	colormap_list = append(colormap_list, "pdb")
	colormap_list = append(colormap_list, "vmd")
	colormap_list = append(colormap_list, "p1")
	colormap_list = append(colormap_list, "p2")
	colormap_list = append(colormap_list, "1")
	colormap_list = append(colormap_list, "2")
	colormap_list = append(colormap_list, "3")
	colormap_list = append(colormap_list, "4")
	colormap_list = append(colormap_list, "5")
	colormap_list = append(colormap_list, "6")
	colormap_list = append(colormap_list, "7")
	colormap_list = append(colormap_list, "8")
	colormap_list = append(colormap_list, "9")
	colormap_list = append(colormap_list, "10")
	colormap_list = append(colormap_list, "11")
	colormap_list = append(colormap_list, "12")
	colormap_list = append(colormap_list, "13")
	colormap_list = append(colormap_list, "14")
	colormap_list = append(colormap_list, "15")
	colormap_list = append(colormap_list, "16")
	colormap_list = append(colormap_list, "17")
	colormap_list = append(colormap_list, "18")
	colormap_list = append(colormap_list, "19")
	colormap_list = append(colormap_list, "20")

	// initial the other colarmap
	gmx_cm := make(map[string]string)
	gmx_cm["Coil"] = "#BBBBBB"
	gmx_cm["B-Sheet"] = "#FF0000"
	gmx_cm["B-Bridge"] = "#000000"
	gmx_cm["Bend"] = "#008000"
	gmx_cm["Turn"] = "#FFFF00"
	gmx_cm["A-Helix"] = "#0000FF"
	gmx_cm["5-Helix"] = "#800080"
	gmx_cm["3-Helix"] = "#808080"

	pdb_cm := make(map[string]string)
	pdb_cm["A-Helix"] = "#E53535"
	pdb_cm["3-Helix"] = "#EC6262"
	pdb_cm["5-Helix"] = "#F28E8E"
	pdb_cm["B-Sheet"] = "#FFF000"
	pdb_cm["B-Bridge"] = "#D4A800"
	pdb_cm["Turn"] = "#00B266"
	pdb_cm["Bend"] = "#00cc00"
	pdb_cm["Coil"] = "#000000"

	vmd_cm := make(map[string]string)
	vmd_cm["A-Helix"] = "#E182E1"
	vmd_cm["3-Helix"] = "#FFA0FF"
	vmd_cm["5-Helix"] = "#E11414"
	vmd_cm["B-Sheet"] = "#FFFF64"
	vmd_cm["B-Bridge"] = "#B4B400"
	vmd_cm["Turn"] = "#469696"
	vmd_cm["Bend"] = "#00CC00"
	vmd_cm["Coil"] = "#000000"

	p1_cm := make(map[string]string)
	p1_cm["A-Helix"] = "#00FF00"
	p1_cm["3-Helix"] = "#2E8857"
	p1_cm["5-Helix"] = "#008000"
	p1_cm["B-Sheet"] = "#0000FF"
	p1_cm["B-Bridge"] = "#3399FF"
	p1_cm["Turn"] = "#FF8C00"
	p1_cm["Bend"] = "#FF0000"
	p1_cm["Coil"] = "#000000"

	p2_cm := make(map[string]string)
	p2_cm["A-Helix"] = "#017186"
	p2_cm["3-Helix"] = "#D16000"
	p2_cm["5-Helix"] = "#FCAEE5"
	p2_cm["B-Sheet"] = "#019E73"
	p2_cm["B-Bridge"] = "#CA9261"
	p2_cm["Turn"] = "#C99260"
	p2_cm["Bend"] = "#F8B0E7"
	p2_cm["Coil"] = "#000000"

	cm_1 := make(map[string]string)
	cm_1["A-Helix"] = "#0000FF"
	cm_1["3-Helix"] = "#0092FF"
	cm_1["5-Helix"] = "#00FFDB"
	cm_1["B-Sheet"] = "#00FF49"
	cm_1["B-Bridge"] = "#49FF00"
	cm_1["Turn"] = "#DBFF00"
	cm_1["Bend"] = "#FF9200"
	cm_1["Coil"] = "#FF0000"

	cm_2 := make(map[string]string)
	cm_2["A-Helix"] = "#8000FF"
	cm_2["3-Helix"] = "#376FF9"
	cm_2["5-Helix"] = "#12C7E6"
	cm_2["B-Sheet"] = "#5BF9C7"
	cm_2["B-Bridge"] = "#A4F99F"
	cm_2["Turn"] = "#EDC76F"
	cm_2["Bend"] = "#FF6F39"
	cm_2["Coil"] = "#FF0000"

	cm_3 := make(map[string]string)
	cm_3["A-Helix"] = "#0000FF"
	cm_3["3-Helix"] = "#006FF9"
	cm_3["5-Helix"] = "#00C7E6"
	cm_3["B-Sheet"] = "#24F9C7"
	cm_3["B-Bridge"] = "#5BF99F"
	cm_3["Turn"] = "#92C76F"
	cm_3["Bend"] = "#C86F39"
	cm_3["Coil"] = "#FF0000"

	cm_4 := make(map[string]string)
	cm_4["A-Helix"] = "#0000D6"
	cm_4["3-Helix"] = "#006F8D"
	cm_4["5-Helix"] = "#00C744"
	cm_4["B-Sheet"] = "#49F900"
	cm_4["B-Bridge"] = "#B6F900"
	cm_4["Turn"] = "#FFC700"
	cm_4["Bend"] = "#FF6F00"
	cm_4["Coil"] = "#FF0000"

	cm_5 := make(map[string]string)
	cm_5["A-Helix"] = "#781C86"
	cm_5["3-Helix"] = "#3F43C8"
	cm_5["5-Helix"] = "#498BC4"
	cm_5["B-Sheet"] = "#6BB18E"
	cm_5["B-Bridge"] = "#9EBE59"
	cm_5["Turn"] = "#D2B240"
	cm_5["Bend"] = "#E67D33"
	cm_5["Coil"] = "#DB2122"

	cm_6 := make(map[string]string)
	cm_6["A-Helix"] = "#3C5793"
	cm_6["3-Helix"] = "#42617B"
	cm_6["5-Helix"] = "#4A7849"
	cm_6["B-Sheet"] = "#79963F"
	cm_6["B-Bridge"] = "#C2BD4A"
	cm_6["Turn"] = "#DDB353"
	cm_6["Bend"] = "#C35C43"
	cm_6["Coil"] = "#BA3D3B"

	cm_7 := make(map[string]string)
	cm_7["A-Helix"] = "#0034F6"
	cm_7["3-Helix"] = "#0B7892"
	cm_7["5-Helix"] = "#3F9E33"
	cm_7["B-Sheet"] = "#8BB718"
	cm_7["B-Bridge"] = "#DBCA28"
	cm_7["Turn"] = "#FEB122"
	cm_7["Bend"] = "#FF7C0D"
	cm_7["Coil"] = "#FC3000"

	cm_8 := make(map[string]string)
	cm_8["A-Helix"] = "#480054"
	cm_8["3-Helix"] = "#4F307F"
	cm_8["5-Helix"] = "#435A8D"
	cm_8["B-Sheet"] = "#347F8E"
	cm_8["B-Bridge"] = "#22A187"
	cm_8["Turn"] = "#3DC36B"
	cm_8["Bend"] = "#93DB35"
	cm_8["Coil"] = "#F3E81C"

	cm_9 := make(map[string]string)
	cm_9["A-Helix"] = "#1F0266"
	cm_9["3-Helix"] = "#173D7E"
	cm_9["5-Helix"] = "#186F89"
	cm_9["B-Sheet"] = "#249689"
	cm_9["B-Bridge"] = "#3EB281"
	cm_9["Turn"] = "#65C874"
	cm_9["Bend"] = "#9ED866"
	cm_9["Coil"] = "#E9E55A"

	cm_10 := make(map[string]string)
	cm_10["A-Helix"] = "#352A87"
	cm_10["3-Helix"] = "#0168E1"
	cm_10["5-Helix"] = "#108ED2"
	cm_10["B-Sheet"] = "#0FAEB9"
	cm_10["B-Bridge"] = "#65BE86"
	cm_10["Turn"] = "#C0BC60"
	cm_10["Bend"] = "#FFC337"
	cm_10["Coil"] = "#FAFB0E"

	cm_11 := make(map[string]string)
	cm_11["A-Helix"] = "#000000"
	cm_11["3-Helix"] = "#00410B"
	cm_11["5-Helix"] = "#0B7A14"
	cm_11["B-Sheet"] = "#359D19"
	cm_11["B-Bridge"] = "#66B9E1"
	cm_11["Turn"] = "#A0CF24"
	cm_11["Bend"] = "#D1E52E"
	cm_11["Coil"] = "#FFFB3B"

	cm_12 := make(map[string]string)
	cm_12["A-Helix"] = "#5648C0"
	cm_12["3-Helix"] = "#7D87EF"
	cm_12["5-Helix"] = "#A6B9FF"
	cm_12["B-Sheet"] = "#CCD7F0"
	cm_12["B-Bridge"] = "#EBD1C1"
	cm_12["Turn"] = "#F3A889"
	cm_12["Bend"] = "#DE6953"
	cm_12["Coil"] = "#B10128"

	cm_13 := make(map[string]string)
	cm_13["A-Helix"] = "#291FCB"
	cm_13["3-Helix"] = "#5772EA"
	cm_13["5-Helix"] = "#97BFF2"
	cm_13["B-Sheet"] = "#CBE3E2"
	cm_13["B-Bridge"] = "#E5DABC"
	cm_13["Turn"] = "#E0A88A"
	cm_13["Bend"] = "#C05B56"
	cm_13["Coil"] = "#88162B"

	cm_14 := make(map[string]string)
	cm_14["A-Helix"] = "#732838"
	cm_14["3-Helix"] = "#A75656"
	cm_14["5-Helix"] = "#CC977E"
	cm_14["B-Sheet"] = "#DCCAA7"
	cm_14["B-Bridge"] = "#CBDBCA"
	cm_14["Turn"] = "#97C303"
	cm_14["Bend"] = "#5596BE"
	cm_14["Coil"] = "#234F8C"

	cm_15 := make(map[string]string)
	cm_15["A-Helix"] = "#2E4EEF"
	cm_15["3-Helix"] = "#6988F2"
	cm_15["5-Helix"] = "#B7C7F8"
	cm_15["B-Sheet"] = "#F1F4FA"
	cm_15["B-Bridge"] = "#FDFDC0"
	cm_15["Turn"] = "#F8EA58"
	cm_15["Bend"] = "#E3993B"
	cm_15["Coil"] = "#01222A"

	cm_16 := make(map[string]string)
	cm_16["A-Helix"] = "#2A4800"
	cm_16["3-Helix"] = "#69A100"
	cm_16["5-Helix"] = "#A4DC00"
	cm_16["B-Sheet"] = "#D7FB00"
	cm_16["B-Bridge"] = "#FEFC00"
	cm_16["Turn"] = "#F6E300"
	cm_16["Bend"] = "#E8AD00"
	cm_16["Coil"] = "#D80000"

	cm_17 := make(map[string]string)
	cm_17["A-Helix"] = "#00204D"
	cm_17["3-Helix"] = "#163A6D"
	cm_17["5-Helix"] = "#4B546C"
	cm_17["B-Sheet"] = "#6D6E72"
	cm_17["B-Bridge"] = "#8E8A78"
	cm_17["Turn"] = "#B3A772"
	cm_17["Bend"] = "#DBC761"
	cm_17["Coil"] = "#FFEA46"

	cm_18 := make(map[string]string)
	cm_18["A-Helix"] = "#000003"
	cm_18["3-Helix"] = "#280B54"
	cm_18["5-Helix"] = "#65156F"
	cm_18["B-Sheet"] = "#9F2A63"
	cm_18["B-Bridge"] = "#D44842"
	cm_18["Turn"] = "#F67D16"
	cm_18["Bend"] = "#FBC127"
	cm_18["Coil"] = "#FCFFA4"

	cm_19 := make(map[string]string)
	cm_19["A-Helix"] = "#000003"
	cm_19["3-Helix"] = "#231251"
	cm_19["5-Helix"] = "#5E177F"
	cm_19["B-Sheet"] = "#982D80"
	cm_19["B-Bridge"] = "#D3436E"
	cm_19["Turn"] = "#F8775D"
	cm_19["Bend"] = "#FEBA80"
	cm_19["Coil"] = "#FCFDBF"

	cm_20 := make(map[string]string)
	cm_20["A-Helix"] = "#0E0787"
	cm_20["3-Helix"] = "#5402A3"
	cm_20["5-Helix"] = "#880AA5"
	cm_20["B-Sheet"] = "#B93389"
	cm_20["B-Bridge"] = "#DB5B68"
	cm_20["Turn"] = "#F38849"
	cm_20["Bend"] = "#FEBC2A"
	cm_20["Coil"] = "#EFF921"

	// os command line args
	if len(os.Args) == 1 {
		fmt.Println("ERROR, No input.")
		fmt.Println("-> type xpm2gpl.exe -h for usage infos.")
		return
	} else if len(os.Args) == 2 {
		inputfile = os.Args[1]
	} else if len(os.Args) == 3 {
		inputfile = os.Args[1]
		colormap = os.Args[2]
	} else if len(os.Args) == 4 {
		inputfile = os.Args[1]
		colormap = os.Args[2]
		output_png = os.Args[3]
	} else if len(os.Args) > 4 {
		fmt.Println("ERROR, wrong number of input, Check your commands.")
		fmt.Println("-> type xpm2gpl.exe -h for usage infos.")
		return
	}

	// initial the colormap and output_png
	if colormap == "" {
		colormap = "xpm"
	}
	if output_png == "" {
		output_png = strings.Split(inputfile, ".")[0] + "_gpl.png"
	}

	if inputfile == "-h" || inputfile == "help" || inputfile == "--help" {
		// usage infos
		fmt.Println("=== Usage : xpm2gpl.exe <filename>.xpm colormap output.png ===")
		fmt.Println("       <filename>.xpm : xpm file to convert to gpl file")
		fmt.Println("             colormap : (optional) colormap for plot")
		fmt.Println("                       xpm : default, original colormap of xpmfile")
		fmt.Println("                       gmx : colormap of gromacs")
		fmt.Println("                       pdb : colormap of pdb")
		fmt.Println("                       vmd : colormap of vmd")
		fmt.Println("                        p1 : 10.1371/journal.pone.0178333")
		fmt.Println("                        p2 : 10.1371/journal.pcbi.1007487")
		fmt.Println("                         1 : Rainbow")
		fmt.Println("                         2 : Rainbow gnuplot_33/13/10")
		fmt.Println("                         3 : Rainbow gnuplot_26/13/10")
		fmt.Println("                         4 : Rainbow gnuplot_22/13/-31")
		fmt.Println("                         5 : Rainbow Mathmatica")
		fmt.Println("                         6 : DarkRainbow Mathmatica")
		fmt.Println("                         7 : Rainbow CETR2")
		fmt.Println("                         8 : Viridis python")
		fmt.Println("                         9 : BlueGreenYellow")
		fmt.Println("                        10 : Parula matlab")
		fmt.Println("                        11 : Avocado")
		fmt.Println("                        12 : Diverging-CoolWarm")
		fmt.Println("                        13 : Thermometer")
		fmt.Println("                        14 : Redbluetones")
		fmt.Println("                        15 : Tmap")
		fmt.Println("                        16 : Ltmap")
		fmt.Println("                        17 : Cividis")
		fmt.Println("                        18 : Inferno")
		fmt.Println("                        19 : Magma")
		fmt.Println("                        20 : Plasma")
		fmt.Println("           output_png : (optional) default <filename>_gpl.png")
		return
	}

	// check inputs
	if len(strings.Split(inputfile, ".")) < 2 || strings.Split(inputfile, ".")[1] != "xpm" {
		fmt.Println("ERROR, ", inputfile, " should be a xpm filename, like ss.xpm, check it.")
		fmt.Println("-> type xpm2gpl.exe -h for usage infos.")
		return
	}
	in_flag := false
	for _, cm := range colormap_list {
		if colormap == cm {
			in_flag = true
		}
	}
	if in_flag == false {
		fmt.Println("ERROR, colormap you input is not defined.")
		fmt.Println("-> type xpm2gpl.exe -h for usage infos.")
		return
	}
	if len(strings.Split(output_png, ".")) < 2 || strings.Split(output_png, ".")[1] != "png" {
		fmt.Println("ERROR, ", output_png, " is not a correct .png filename, input like : ss.png.")
		fmt.Println("-> type xpm2gpl.exe -h for usage infos.")
		return
	}

	// read file and turn into lines
	buf, err := ioutil.ReadFile(inputfile)
	if err != nil {
		fmt.Println(err.Error())
		fmt.Println("ERROR in reading ", inputfile)
		return
	}
	str_buf := strings.Trim(string(buf), " ")
	str_buf = strings.Trim(str_buf, "\n")
	lines := strings.Split(str_buf, "\n")
	//fmt.Println(len(lines))

	// initial the data args
	x_label := ""
	y_label := ""
	xpm_cm := make(map[string]string)
	structure_note_relations := make(map[string]string)
	note_number := make(map[string]string)
	x_axis := make([]string, 0)
	y_axis := make([]string, 0)
	data := make([][]string, 0)

	// initial xpm_cm and structure_note_relations
	xpm_cm["Coil"] = "#000000"
	xpm_cm["B-Sheet"] = "#000000"
	xpm_cm["B-Bridge"] = "#000000"
	xpm_cm["Bend"] = "#000000"
	xpm_cm["Turn"] = "#000000"
	xpm_cm["A-Helix"] = "#000000"
	xpm_cm["5-Helix"] = "#000000"
	xpm_cm["3-Helix"] = "#000000"
	structure_note_relations["Coil"] = " "
	structure_note_relations["B-Sheet"] = " "
	structure_note_relations["B-Bridge"] = " "
	structure_note_relations["Bend"] = " "
	structure_note_relations["Turn"] = " "
	structure_note_relations["A-Helix"] = " "
	structure_note_relations["5-Helix"] = " "
	structure_note_relations["3-Helix"] = " "

	// deal with xpm file contents
	cm_count := 0 // count the legend line number
	for _, line := range lines {
		items := strings.Split(line, " ")
		// len(items) > 1 means not the plot data line
		if len(items) > 1 {
			// parse x_label, y_label, x-axis, y-axis
			if items[1] == "x-label:" {
				x_label = strings.Split(line, "\"")[1]
			} else if items[1] == "y-label:" {
				y_label = strings.Split(line, "\"")[1]
			} else if items[1] == "x-axis:" {
				for ind, item := range items {
					ite := strings.Trim(item, " ")
					if ind != 0 && ind != 1 && ind != len(items)-1 && len(ite) != 0 {
						x_axis = append(x_axis, ite)
					}
				}
			} else if items[1] == "y-axis:" {
				for ind, item := range items {
					ite := strings.Trim(item, " ")
					if ind != 0 && ind != 1 && ind != len(items)-1 && len(ite) != 0 {
						y_axis = append(y_axis, ite)
					}
				}
			}
			// parse xpm_cm and structure_note_relations
			if len(items) == 8 { // mind the fucking space key
				if len(items[3]) == 7 && len(items[0]) == 2 {
					// construct the relations between notes and numbers
					note_number[strings.Trim(items[0], "\"")] = strconv.Itoa(cm_count)
					cm_count += 1
					if items[6] == "\"Coil\"" {
						xpm_cm["Coil"] = items[3]
						structure_note_relations["Coil"] = strings.Trim(items[0], "\"")
					} else if items[6] == "\"B-Sheet\"" {
						xpm_cm["B-Sheet"] = items[3]
						structure_note_relations["B-Sheet"] = strings.Trim(items[0], "\"")
					} else if items[6] == "\"B-Bridge\"" {
						xpm_cm["B-Bridge"] = items[3]
						structure_note_relations["B-Bridge"] = strings.Trim(items[0], "\"")
					} else if items[6] == "\"Bend\"" {
						xpm_cm["Bend"] = items[3]
						structure_note_relations["Bend"] = strings.Trim(items[0], "\"")
					} else if items[6] == "\"Turn\"" {
						xpm_cm["Turn"] = items[3]
						structure_note_relations["Turn"] = strings.Trim(items[0], "\"")
					} else if items[6] == "\"A-Helix\"" {
						xpm_cm["A-Helix"] = items[3]
						structure_note_relations["A-Helix"] = strings.Trim(items[0], "\"")
					} else if items[6] == "\"3-Helix\"" {
						xpm_cm["3-Helix"] = items[3]
						structure_note_relations["3-Helix"] = strings.Trim(items[0], "\"")
					} else if items[6] == "\"5-Helix\"" {
						xpm_cm["5-Helix"] = items[3]
						structure_note_relations["5-Helix"] = strings.Trim(items[0], "\"")
					}
				}
			}
		}
		// parse data and convert note to number
		if len(items) == 1 && string(line[0]) == "\"" {
			data_line := make([]string, 0)
			line4data := strings.Trim(line[0:len(line)-1], "\"")
			for _, ch := range line4data {
				data_line = append(data_line, note_number[string(ch)])
			}
			data = append(data, data_line)
		}

	}
	//fmt.Println(x_label)
	//fmt.Println(y_label)
	//fmt.Println(x_axis)
	//fmt.Println(y_axis)
	fmt.Println("number of legends : ", cm_count)

	// check xpm_cm and maybe change it to gmx
	xpm_cm_value_check := 8
	for _, value := range xpm_cm {
		if value != "#000000" {
			xpm_cm_value_check -= 1
		}
	}
	if xpm_cm_value_check > 8-cm_count+1 {
		xpm_cm = gmx_cm
		fmt.Println("Warning, colarmap can't be parsed from xpm file properly, using gmx for instead.")
	}
	// apply the colormap input
	if colormap != "xpm" {
		if colormap == "gmx" {
			xpm_cm = gmx_cm
		} else if colormap == "pdb" {
			xpm_cm = pdb_cm
		} else if colormap == "vmd" {
			xpm_cm = vmd_cm
		} else if colormap == "p1" {
			xpm_cm = p1_cm
		} else if colormap == "p2" {
			xpm_cm = p2_cm
		} else if colormap == "1" {
			xpm_cm = cm_1
		} else if colormap == "2" {
			xpm_cm = cm_2
		} else if colormap == "3" {
			xpm_cm = cm_3
		} else if colormap == "4" {
			xpm_cm = cm_4
		} else if colormap == "5" {
			xpm_cm = cm_5
		} else if colormap == "6" {
			xpm_cm = cm_6
		} else if colormap == "7" {
			xpm_cm = cm_7
		} else if colormap == "8" {
			xpm_cm = cm_8
		} else if colormap == "9" {
			xpm_cm = cm_9
		} else if colormap == "10" {
			xpm_cm = cm_10
		} else if colormap == "11" {
			xpm_cm = cm_11
		} else if colormap == "12" {
			xpm_cm = cm_12
		} else if colormap == "13" {
			xpm_cm = cm_13
		} else if colormap == "14" {
			xpm_cm = cm_14
		} else if colormap == "15" {
			xpm_cm = cm_15
		} else if colormap == "16" {
			xpm_cm = cm_16
		} else if colormap == "17" {
			xpm_cm = cm_17
		} else if colormap == "18" {
			xpm_cm = cm_18
		} else if colormap == "19" {
			xpm_cm = cm_19
		} else if colormap == "20" {
			xpm_cm = cm_20
		}
		// else ......
	}
	//fmt.Println(xpm_cm)

	// check structure_note_relations
	struc_note_check := 0
	for _, value := range structure_note_relations {
		if value == " " {
			struc_note_check += 1
		}
	}
	if struc_note_check != 8-cm_count {
		fmt.Println("ERROR, can't parse structures' note correctly, check xpm file.")
		fmt.Println(structure_note_relations)
		return
	}
	//fmt.Println(note_number)
	//for _, line := range data {
	//	fmt.Println(line)
	//}

	/////////////////////////////
	// to write gpl
	/////////////////////////////
	gpl_lines := make([]string, 0)
	gpl_lines = append(gpl_lines, "set term png")
	gpl_lines = append(gpl_lines, "set output \""+output_png+"\"")
	gpl_lines = append(gpl_lines, "unset colorbox")
	pal_line := "set pal defined("
	for i := 0; i < cm_count; i++ {
		color := ""
		for note, value := range note_number {
			if value == strconv.Itoa(i) {
				for stru, not := range structure_note_relations {
					if note == not {
						color = xpm_cm[stru]
						break
					}
				}
				break
			}
		}
		if color == "" {
			fmt.Println("Warning, color parsed from xpm_cm fail !")
			color = "#000000"
		}
		pal_line = pal_line + strconv.Itoa(i) + " \"" + color + "\","
	}
	pal_line = strings.Trim(pal_line, ",") + ")"
	gpl_lines = append(gpl_lines, pal_line)
	gpl_lines = append(gpl_lines, "")
	// add data lines
	gpl_lines = append(gpl_lines, "$data <<EOD")
	if len(y_axis) != len(data) {
		fmt.Println("ERROR, plot data is not consistent with y_axis.")
	}
	if len(x_axis) != len(data[0]) {
		fmt.Println("ERROR, plot data is not consistent with x_axis.")
	}
	y_len := len(y_axis)
	for x := 0; x < len(x_axis); x++ {
		for y := 0; y < len(y_axis); y++ {
			gpl_lines = append(gpl_lines, x_axis[x]+" "+y_axis[y]+" "+data[y_len-1-y][x])
		}
	}
	gpl_lines = append(gpl_lines, "EOD")
	gpl_lines = append(gpl_lines, "")
	// add tail part of gpl file
	gpl_lines = append(gpl_lines, "#set tmargin at screen 0.95")
	gpl_lines = append(gpl_lines, "#set bmargin at screen 0.20")
	gpl_lines = append(gpl_lines, "#set rmargin at screen 0.85")
	y_posi := [10]string{"0.92", "0.82", "0.72", "0.62", "0.52", "0.42", "0.32", "0.22", "0.12", "0.02"}
	y_posi_i := 0
	for note, _ := range note_number {
		for key, not := range structure_note_relations {
			if not == note {
				label_line := "#set label \" " + key + " \" at screen 0.85,"
				label_line = label_line + y_posi[y_posi_i]
				label_line = label_line + " left textcolor rgb " + "\""
				label_line = label_line + xpm_cm[key] + "\""
				gpl_lines = append(gpl_lines, label_line)
				y_posi_i += 1
			}
		}
	}
	font_line := "set term pngcairo enhanced truecolor font \"Arial,85\" "
	font_line = font_line + "fontscale 1 linewidth 20 pointscale 5 size 5000,6000"
	gpl_lines = append(gpl_lines, font_line)
	gpl_lines = append(gpl_lines, "set tics out nomirror;")
	gpl_lines = append(gpl_lines, "set key out reverse Left spacing 2 samplen 1/2")
	gpl_lines = append(gpl_lines, "set xl\""+x_label+"\"; set yl\""+y_label+"\";")
	gpl_lines = append(gpl_lines, "plot [0:"+x_axis[len(x_axis)-1]+"] [1:"+y_axis[len(y_axis)-1]+"] $data u 1:2:3 w imag notit, \\")
	for note, _ := range note_number {
		for key, not := range structure_note_relations {
			if not == note {
				lin := "-1 w p ps 3 pt 5 lc rgb \"" + xpm_cm[key] + "\" t\"" + key + "\", \\"
				gpl_lines = append(gpl_lines, lin)
			}
		}
	}
	gpl_lines[len(gpl_lines)-1] = strings.Trim(gpl_lines[len(gpl_lines)-1], "\\")
	gpl_lines[len(gpl_lines)-1] = strings.Trim(gpl_lines[len(gpl_lines)-1], " ")
	gpl_lines[len(gpl_lines)-1] = strings.Trim(gpl_lines[len(gpl_lines)-1], ",")

	//for _, line := range gpl_lines {
	//	fmt.Println(line)
	//}
	// write to gpl file
	fileName := strings.Split(inputfile, ".")[0] + ".gpl"
	dstFile, err := os.Create(fileName)
	if err != nil {
		fmt.Println(err.Error())
		return
	}
	defer dstFile.Close()
	out2file := strings.Join(gpl_lines, "\n")
	dstFile.WriteString(out2file + "\n")
	fmt.Println("Content has been successfully written into " + fileName)
	fmt.Println("<<<")
}
