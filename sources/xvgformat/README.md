#### description 

此脚本用于格式化某些GROMACS生成的xvg文件，可以生成去掉大部分无关信息的xvg文件以及csv文件。

#### command

```shell
python xvgformat.py inputfile1.xvg inputfile2.xvg ...
```

add '-c' to generate csv file.

```shell
python xvgformat.py inputfile1.xvg inputfile2.xvg -c
```
