#### command

```shell
python3 energy_compute.py prolig.xvg pro.xvg lig.xvg -foutput
```

#### usage

用于计算蛋白配体之间的相互作用，计算方法参考Jerkwin博客：

> https://jerkwin.github.io/2019/09/06/使用GROMACS计算分子间相互作用/

通过 E = prolig - pro - lig 求得能量

输入一共需要三个文件：蛋白配体的能量xvg文件，蛋白的能量xvg文件，配体的能量xvg文件
每个输入文件需要包含（且只能包含）四列能量数据，用于各项能量以及总库伦势能和总能量的计算
要求四列能量数据的顺序为 

> LJ-SR  | Disper.corr. | Coulomb-SR | Coul.-recip.

生成数据文件的各列名如下

> LJ-SR  | Disper.corr. | Coulomb-SR | Coul.-recip. | ETOTAL | COULOMB | LJ-Total

生成的结果文件为 energy_results_.xvg可以用xvgshow.py或xvg_compare.py可视化

#### dependency

1. seaborn
2. matplotlib