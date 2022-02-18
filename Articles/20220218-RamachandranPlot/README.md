# 绘制拉式图

### Intro

拉式图（Ramachandran图）常用于表征蛋白质二级结构的情况。

对于分子模拟中的蛋白质而言，目前也有很多工具和方法可以绘制拉式图。



### rama命令

对于GROMACS来说，其自带了`rama`命令，可以计算轨迹中蛋白质的phi-psi角度并输出到xvg文件。

```bash
gmx rama -s md.tpr -f md.xtc -o rama.xvg
```

执行命令之前需要注意做周期性校正，然后选组的时候就选蛋白质，可以适当设置`-dt`参数，以免生成过大的xvg文件。

xvg文件内部注释部分写得很详细，数据部分分为三列，Phi、Psi和氨基酸名字序号。可以直接用这些信息绘制散点图进而得到拉式图（需要自己加个背景）。



### pyrama

PyRAMA(https://github.com/gerdos/PyRAMA)是一个可以读入pdb文件并绘制拉式图的python第三方库。使用非常简单。pip安装完pyrama之后，把github仓库中bin文件夹下的pyrama.py文件下载下来并添加到环境变量就可以了。

PyRAMA可以同时绘制四幅拉式图（General、GLY、PRE-PRO以及PRO）的，基本上涵盖到了日常的需求。

PyRAMA有个好处是可以读入多帧的pdb文件，因而对于模拟轨迹导出的多帧pdb文件比较适用。

PyRAMA会把在合理范围内的店绘制成蓝色(normals)，而处在不合理位置的点则是红色(outliers)。

![fig_0](C:\Users\hhhhh\Desktop\databank\公众号\20220218\fig_0.png)

除此之外，PyRAMA的代码加起来还没有200行，自己要改改也非常容易！



### mdanalysis

MDAnalysis，一个很有名的分子模拟分析包，自然也是能做拉式图的。

官方的文档写得很详细了，可以参考 https://userguide.mdanalysis.org/1.1.1/examples/analysis/structure/dihedrals.html

只需要把gro文件、xtc文件准备一下，几行代码就可以搞定了。

![Figure_1](C:\Users\hhhhh\Desktop\20220218_RAMA\Figure_1.png)

```python
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import Ramachandran

u = mda.Universe("test.gro", "test.xtc")
protein = u.select_atoms("protein")
R = Ramachandran(protein).run()
fig, ax = plt.subplots()
R.plot(ax=ax, color='k', marker='s', ref=True)
plt.show()
```

同样你也可以绘制Gly和Pro的拉式图，选择原子的时候用残基名字做关键词选就可以了。



### Others

今天没有Others。







