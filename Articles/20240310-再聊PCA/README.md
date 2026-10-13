## 再聊PCA

#### Intro

Wentao大佬在bohrium平台上发的介绍PCA的文章，可谓酣畅淋漓、尽善尽美了，当真是学习分子动力学轨迹PCA的最佳教程。

> 细看成岭又成峰——分子动力学轨迹的主成分分析 https://nb.bohrium.dp.tech/detail/6316334641

珠玉在前，这篇小文就权做一点小小的补充。



#### PCA projection

大佬的文章详细讲述了PCA的数学原理，以及PCA在轨迹上的分析和应用，包括如何执行主成分分析、如何获得前几个主成分的方差占比等等。

其文中讲到：

> 为了更好的观察主成分所代表的原子运动,原文作者根据PC1和PC2特征向量重构分子的三维结构(PathReducer通过将PC转换回全维空间来生成分子几何的笛卡尔坐标，从而用PC重建结构),并将最大的两主成分根据当前帧和下一帧的变化,将原子运动方向投影回结构上

具体的将主成分投影回结构上的过程文中并未详述，这里简单做个补充。

GROMACS软件的`anaeig`命令是集成了将主成分投影回结构上的功能的，只需要用`-extr`参数输出到结构文件即可。输出文件中包含了两帧结构，分别是PC的极小值和极大值在结构上的投影，也即可以用于表示沿PC方向上蛋白质构象的变化。

MDAnalysis本身也有一个PCA的模块，可以对轨迹执行PCA，也可以将前几个PC投影到结构上，非常的方便。

我们这里可以探讨一下将PC投影回结构的这个过程，还是挺有趣的。

使用scikit-learn对轨迹执行PCA之后，我们会得到一个特征向量矩阵和轨迹不同帧的特征值。假设是一个N个原子F帧的体系，PCA降到3维，则会得到一个$3*3N$的特征向量矩阵，这个矩阵存储在`pca.components_`里面，其中第一行也就是第一个PC的特征向量，第二行对应第二个PC的特征向量。而特征值的结果，则会是$F*3$的，第一列是第一个PC的特征值，第二列是第二个PC的特征值。同时，轨迹的平均值还会被存储到`pca.mean_`中。

那么要计算PC在构象上的投影，则只需要把对应的特征值向量和特征向量去做点积，就会得到沿PC方向构象的变化，再加上平均值，即可得到某一个PC在构象上的投影。通常情况下，我们只需要一个PC的特征值的极值就行了，利用最大和最小值，组成向量和该PC的特征向量去做点积，就能得到该PC的最大和最小值相较于平均值的“偏移”。

说起来还是非常简单的，贴一些代码方便理解：

```python
from sklearn.decomposition import PCA as skpca

## 执行PCA
pca_model = skpca(n_components=3)
data = pca_model.fit_transform(protein_coors_byframe)

## 拿到PC1的极大极小值
pc1_i_s = np.argsort(data[:, 0]).tolist()
pc1_componet = pca_model.components_[0].reshape(1, -1)
pc1_minmax = np.array([data[pc1_i_s[0], 0], data[pc1_i_s[-1], 0]]).reshape(-1, 1)

## 计算投影
pc1_reconstructed = np.dot(pc1_minmax, pc1_componet) + pca_model.mean_ 

## 将投影的维度整理到类似于原子坐标的维度形式并输出到pdb
pc1_reconstructed = np.reshape(pc1_reconstructed, (2, -1, 3))
with mda.Writer("pc1_proj.pdb", multiframe=True) as fo:
    for i in range(2):
        pro.atoms.positions = pc1_reconstructed[i]
        fo.write(pro)
```

其中`protein_coors_byframe`是F帧轨迹中N个原子的坐标展平得到的，也即shape是(F, 3N)。这里只将PC1的极大极小值投影到了结构上，也可以将PC1的所有特征值从小到大投影到结构上，那就得到了构象沿着PC1方向上的变化过程。



#### Cosine Content

除了大佬文章中提到的成分占比，余弦含量也可以作为PCA分析的一个衡量标准。一般来讲，蛋白的运动也会受到扩散的影响，如果轨迹没有收敛，那么很可能PCA得到的前几个主成分非常接近于余弦函数（随机扩散表现出来就是余弦），因而通过计算前几个主成分与余弦的相似度，也即余弦含量，一方面可以表征体系是否收敛，另一方面也可以用来衡量PCA分析的质量。

要计算某一个PC的余弦含量，代码也非常简单：

```python 
from scipy.integrate import simpson # 积分
i = 0 ## for pc1 
t = np.arange(len(pc_datas))
T = len(pc_datas)
cos = np.cos(np.pi * t * (i+1) / T)
cosine_content = ((2.0/T) * (simpson(cos*pc_datas[:, i]))**2 / simpson(pc_datas[:, i]**2))
```

对于GROMACS来说，也有`analyze`命令的`-cc`参数可以计算余弦含量。

余弦含量的值越接近于1，则表明该PC可能更代表随即扩散，一般小于0.7可能就还好。

更多资料可以参考`gmx analyze`文档，以及 [doi:10.1103/PhysRevE.65.031910](https://doi.org/10.1103/PhysRevE.65.031910).



## Others

DIP集成了轨迹PCA以及骨架二面角PCA分析方法，感兴趣的朋友可以查看DIP文档：https://duivyprocedures-docs.readthedocs.io/en/latest/



科研顺利诸位~

