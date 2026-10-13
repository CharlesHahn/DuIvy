## 📢 DuIvyTools v0.6.0 发布：分子动力学模拟分析利器的再升级！

DuIvyTools，专为GROMACS结果文件设计的分子动力学模拟分析和可视化工具，最新版本v0.6.0现已发布！🎉 这次更新不仅引入了多项功能改进，还进一步提升了工具的性能和易用性，帮助研究者更高效地处理和分析分子动力学数据。

### 🚀 DuIvyTools 的用途

DuIvyTools 是一个简洁而强大的工具集，用于处理和可视化GROMACS生成的结果文件，如 `.xvg` 和 `.xpm`。该工具专为加速分子动力学模拟数据的分析而设计，支持多种操作，如数据显示、平均值计算、能量计算等。研究人员可以通过简便的命令行操作快速生成各种可视化结果，极大地提高了数据分析的效率。

### 🛠️ v0.6.0 更新亮点

在最新的 v0.6.0 版本中，DuIvyTools 迎来了以下重要更新

1. 新增参数`-xp`、`-yp`、`-zp`用于对数据做加减操作，可用于自定义坐标轴刻度标签。
2. 增加`--legend_ncol`参数用于指定图例的列数，方便绘图，仅对matplotlib绘图引擎有效。
3. 在DIT.mplstyle中新增了`figure.figsize`参数，方便用户自定义图片尺寸
4. 新增了`--x_numticks`、`--y_numticks`、`--z_numticks`参数，用于指定坐标轴刻度标签的数量，对matplotlib绘图引擎有效；至于plotly和gnuplot，则只对X和Y轴有重复数据的矩阵图绘制有效。
5. 增加了对无列名数据文件的支持；在如pmf.xvg中，包含了多列没有列名的数据，原本DIT会将之当成字符列，现在DIT可以自动处理成数字并可视化了。
6. 修复了计算置信区间的bug，`-smv`现在可以赋予参数，默认是显示原始数据作为背景，滑动平均作为前景（这应该更常用）；如果`-smv CI`，则背景为置信区间。
7. 修复了一系列已知bugs，并优化了一些小问题。



新的滑动平均显示：

```bash
 dit xvg_compare -f rmsd.xvg -c 1 -smv -t "" -l "RMSD" -xs 0.001 -x "Time(ns)"
```



![rmsd-smv](C:\Users\hhhhh\Desktop\20240818-DIT更新\rmsd-smv.png)





自定义X和Y刻度数量：

```bash
dit xpm_show -f dccm.xpm -t "" --x_numticks 5 --y_numticks 5 --z_numticks 5 -zmin -1 -z "DCCM" -cmap bwr
```



![xpm_show](C:\Users\hhhhh\Desktop\20240818-DIT更新\xpm_show.png)





### 📦 如何获取

DuIvyTools 作为一个 Python 3 库，可以通过以下命令轻松安装：

```
pip install DuIvyTools -i https://pypi.tuna.tsinghua.edu.cn/simple
```

安装完成后，只需在命令行中输入 `dit`，即可开始探索丰富的功能。

### 📚 详细文档

关于DuIvyTools的更多信息、使用教程和命令列表，欢迎访问其[官方文档](https://duivytools.readthedocs.io/en/latest/DIT.html)。

### 🧪 适用场景

无论你是进行分子动力学模拟数据分析，还是可视化模拟结果，DuIvyTools 都能为你提供便利的支持。通过直观的命令操作，你可以快速生成所需的分析图表，大幅减少手动数据处理的时间，让你更专注于科研本身。

赶快更新到 v0.6.0 版本，体验全新的功能提升吧！🎉