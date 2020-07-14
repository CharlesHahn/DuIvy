#### description
用于gromacs生成的xvg文件以及整理过后的xvg文件中的数据求平均值，可以自定义求平均的列和每列的起始行和结束行

#### command
```shell
python xvg_average.py xvgfilename full 10 100
```

# usage 
    xvg_average.py XVG_Filename column_select start_index end_index
        column_select -> e.g. 1,3 or full; default full
        start_index   -> optional, e.g. 10; default 0
        end_index     -> optional, but must be POSITIVE, NO -1 ! default End