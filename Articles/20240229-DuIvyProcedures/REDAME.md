# DuIvyProcedures —— 分子动力学模拟分析软件

## Intro

DuIvyProcedures (DIP) ，由杜艾维团队开发的分子动力学模拟分析软件，现在正式发布啦！

DuIvyProcedures支持多种常见的分析方法， 可以批量对多条轨迹进行多种分析，一条命令出图出数据，支持Linux、Windows等平台。

DIP的文档请见：https://duivyprocedures-docs.readthedocs.io/en/latest/index.html



## DuIvyProcedures

以下介绍DIP中几种常见分析方法的分析组件，更多组件和细节请参考文档。



### 动态互相关矩阵DCCM

DCCM是一种非常常见的针对蛋白质残基之间运动相关性的分析手段。DIP中集成了两个分析模块，可以进行DCCM的计算：`gmx_DCCM`模块可以依赖GROMACS计算协方差矩阵并由DIP转换成DCCM；`DCCM`模块则包含了改进的算法，可以以极高的速度直接从轨迹中读取坐标并计算出DCCM。

DIP会保存DCCM图片如下，同时也会保存相应的矩阵数据以供用户使用DuIvyTools或其它工具绘图。

![correlation_matrix](C:\Users\hhhhh\Desktop\20240229-DIP\correlation_matrix.png)



### 二面角主成分分析dPCA

DIP集成了两个可以计算二面角主成分分析的工具，一个是`gmx_dPCA`，一个是`PCA`组件，前者可以对蛋白质骨架的二面角进行主成分分析，后者也可以通过参数设置使之计算蛋白质phi、psi两个角度的主成分分析。DIP在完成dPCA之后会默认抽提前三个主成分并绘制两两主成分之间的散点图以及相应的数据文件，可直接作为`gmx_FEL`组件的输入并绘制成自由能形貌图。除此之外，DIP中所有PCA计算的组件都会对结果进行主成分含量计算和余弦含量（cosine content）计算，帮助用户确定PCA结果。

![dpc12_scatter](C:\Users\hhhhh\Desktop\20240229-DIP\dpc12_scatter.png)



### 主成分分析 PCA

DIP自然也包含基于坐标的主成分分析组件，`gmx_PCA`和`PCA`，这两个组件会基于选中原子的坐标进行主成分分析并抽提相应的数据，绘制成散点图，也会计算主成分占比和余弦含量。

![pc12_scatter](C:\Users\hhhhh\Desktop\20240229-DIP\pc12_scatter.png)

除了PCA这种降维手段之外，DIP中还集成了tSNE、tICA、UMAP等可以对坐标或者二面角进行降维的组件。



### 自由能形貌图FEL

前面提到DIP能进行各种降维分析，正好这些组件的结果可以与`gmx_FEL`组件组合使用来绘制基于降维结果的自由能形貌图，例如基于PCA的自由能形貌图：

![gibbs](C:\Users\hhhhh\Desktop\20240229-DIP\gibbs.png)

除了绘制自由能形貌图之外，DIP还能自动找出自由能形貌图中的最低值点，并从轨迹中导出相应的帧，免去了诸多麻烦琐碎的事情。例如DIP就在下图中标识除了最低值点的位置以及其对应的时间，同时相应的帧也被输出到pdb文件中。

![Scatter_FEL](C:\Users\hhhhh\Desktop\20240229-DIP\Scatter_FEL.png)



### 蛋白质二级结构

DIP同样集成了模拟过程中蛋白质二级结构含量分析的组件，该组件调用GROMACS命令进行蛋白质二级结构的计算，最终绘制蛋白质二级结构变化图和含量图：

![dssp_Protein](C:\Users\hhhhh\Desktop\20240229-DIP\dssp_Protein.png)



### 氢键分析

DIP可以进行两个原子组之间的或者一个原子组之内的氢键分析，不仅能够分析氢键数量、氢键占有率，还能够计算氢键的平均距离和角度等数据。不仅可以绘制出漂亮的图，还会整理出数据表格，方便后续处理。

氢键的数量变化图：

![hbnum](C:\Users\hhhhh\Desktop\20240229-DIP\hbnum.png)

氢键的时间占有率图：

![top_hbmap](C:\Users\hhhhh\Desktop\20240229-DIP\top_hbmap.png)

此图和距离、角度图都只包含了占有率最高的几个氢键，用户可以自定义呈现的氢键数量。这里的纵坐标是氢键的索引，相应的氢键名称等数据也可以在相应的数据文件中找到。

氢键的距离变化图：

![hbdistall](C:\Users\hhhhh\Desktop\20240229-DIP\hbdistall.png)



### PiStacking

DIP还能计算Pi-Pi相互作用，支持让DIP自动检测可能的芳香环，也可以自定义芳香环原子索引，然后计算芳香环之间的PiStacking；支持输出芳香环距离、角度、中心偏移等性质随时间的变化，还有每个PiStacking的占有率数据；同样也会输出相应的图和数据文件。

例如所有PiStacking的时间占有率图：

![PiStacking_Existence_Map](C:\Users\hhhhh\Desktop\20240229-DIP\PiStacking_Existence_Map.png)

芳香环中心距离随时间的变化：

![PiStacking_Distances](C:\Users\hhhhh\Desktop\20240229-DIP\PiStacking_Distances.png)

此组件还能分辨PiStacking的类型，并输出相应的PiStacking类型时间占有率图：

![PiStacking_Type_Map](C:\Users\hhhhh\Desktop\20240229-DIP\PiStacking_Type_Map.png)

这里占有率图上的纵坐标也是PiStacking的索引，相应的PiStacking名称数据等也会保存在相应的数据文件中。



### 盐桥SaltBridge

同样，DIP也支持计算盐桥。DIP可以自己寻找体系中带电的原子组，也支持用户自己定义用于计算盐桥的原子组。类似于前面的氢键分析和PiStacking分析，盐桥分子模块也会计算盐桥的距离、时间占有率等数据，并输出相应的图和数据文件。![SaltBridge_Existence_Map](C:\Users\hhhhh\Desktop\20240229-DIP\SaltBridge_Existence_Map.png)



### 残基距离接触矩阵RDCM

DIP还有一个强大的组件，支持计算残基距离接触矩阵，包括残基质心距离接触矩阵、残基C-alpha原子距离接触矩阵，以及残基间最短距离接触矩阵。

DIP会按照设置的参数输出不同时间帧的RDCM、帧间RDCM。之后，该组件还会基于残基距离接触矩阵进行一系列的分析，包括基于RDCM计算体系的均方根偏差RMSD、均方根波动RMSF等数据，还会计算残基间距离与时间等变量的相关性矩阵，还可以基于RDCM进行主成分分析PCA、对残基或者时间进行聚类等等。

该组件还支持通过设定的阈值计算contact和Encounter，计算并输出相应的时间占有率等数据。

这里列举本组件的一两个图。例如平均RDCM：

![Distance_Average](C:\Users\hhhhh\Desktop\20240229-DIP\Distance_Average.png)

基于RDCM对残基做的层次聚类图：

![RDCM_Residues_dendrogram](C:\Users\hhhhh\Desktop\20240229-DIP\RDCM_Residues_dendrogram.png)



### 其它组件

当然！DIP还包括了一些更加基础的组件，例如：

- 均方根误差RDSD
- 均方根波动RMSF
- 回旋半径Gyration
- 密度Density
- 基于坐标的聚类
- 径向分布函数RDF
- 最短路径图SPM
- ......





## Price

DuIvyProcedures是一款商业软件，但**对于学术授权和国内高校学生有特别的优惠**，欢迎大家了解和购买！

具体的授权价格和软件信息请参考DIP的文档：https://duivyprocedures-docs.readthedocs.io/en/latest/DuIvyProcedures.html





## Others

DIP的设计目标是追求更丰富的分析组件、更深入的分析、以及更完善的分析流程，力图做成一站式、自动化、简单化的分子动力学模拟分析软件。

在后续的维护和升级中，DIP还会添加新的分析组件，敬请期待。











