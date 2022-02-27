## Description

- xvgformat 
  - 用于对gromacs生成的xvg文件进行格式化，去掉很多不必要的字段，也可以转换为csv
- energy_compute 
  - 用于计算模拟体系中两种物质之间的相互作用
- xvg_average
  - 用于对xvg文件的指定部分求平均，可以输出各项参数在一定时间内的平均值
- xvg_show
  - 用于对xvg结果进行可视化，绘制各项参数随时间变化趋势
- xvg_compare 
  - 用于对不同的xvg文件及利用xvgformat之后的数据进行对比，并可视化展示
- dssp
  - 包含对gmx的`do_dssp`命令的一些结果文件进行处理的脚本，以及一些go语言程序
- xpm2png
  - 将xpm文件转换为图片
- other
  - 一些乱七八糟的脚本
- pipi_dist_ang
  - 计算两个平面环（苯环）之间的几何中心距离和平面夹角

---
本Sources目录下，所有代码、数据等皆按照GPLv2协议开源。
