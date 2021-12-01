# 蛋白配体相互作用能计算

### Intro

`gmx energy`能够计算的分子间相互作用主要分为两块：范德华相互作用和库伦相互作用。其中范德华相互作用随距离的增加而衰减很快，因而可以认为是一种短程的相互作用；而库伦相互作用随距离的增加而衰减较慢，作用的距离较长，是一种长程的相互作用。

范德华力，我们通常认为它包括了分子之间的色散力、诱导力和Pauli交换互斥力等。范德华力最初是从实验中归纳出来的，定义颇为模糊。GMX中使用Lennard-Jones potential（LJ）来描述分子间的范德华力。LJ只包含了**色散力**和**Pauli互斥力**两部分，诱导力通常是非常弱的。在不存在强极性物质的情况下，LJ已经足以描述体系了。

![LJ_graph](C:\Users\hhhhh\Desktop\databank\公众号\20211121\LJ_graph.png)

Fig 1. LJ衰减很快（GMX手册）

![LJ](C:\Users\hhhhh\Desktop\databank\公众号\20211121\LJ.png)

Fig 2. LJ的计算公式（GMX手册）

除了色散力和Pauli排斥力之外，GMX还考虑的一项非键相互作用是**库伦相互作用**。

在通常的蛋白配体模拟过程中，GMX处理LJ的方式通常都是截断，对势垒累加，分解方便。截断之内的能量被命名为LJ-SR，也即短程LJ。超过截断的部分被命名为Disper.corr，以色散校正项的形式出现。而GMX处理库伦相互作用的方式主要有两种：截断和PME。对于蛋白配体体系而言，我想应当是PME使用更多，也更方便。如果使用截断计算库伦相互作用的话，分解也是简单的。如果使用PME的话，截断距离内的库伦相互作用被记为Coulomb-SR（实空间内的库伦），超过截断的部分被记为Coul.recip.（倒易空间内的库伦），这两者加起来才是完整的库仑相互作用。在GMX中，使用PME计算的库仑相互作用不能直接分解，需要利用相互作用的定义来进行计算，也即**两个分子间的能量等于两个分子的总能量减去两个分子各自的能量**：

$$ E_{between A and B} = E_{complex AB} - E_{A} - E_{B}$$

请大家参考更多资料：

1. https://zhuanlan.zhihu.com/p/66008689
2. http://bbs.keinsci.com/thread-13147-1-1.html
3. https://jerkwin.github.io/2019/09/06/使用GROMACS计算分子间相互作用
4. GROMACS manual

对于蛋白配体体系，下文介绍两种方法计算蛋白和配体之间的能量（或许对于其它类似体系同样适用）：一种利用**相互作用的定义**来计算蛋白配体键的相互作用，一种利用**截断**的方法来计算。

计算的案例基于一个100ns的蛋白配体模拟轨迹，正方体盒子，边长9nm；跑主体模拟的时候，vdwtype设置的是默认的cut-off，coulombtype为PME，截断都设置为1.0nm的，同时为了方便使用GPU，没有设置能量检测组。软件为GMX2019.5 。



### Method 1

根据分子间相互作用的定义，我们需要计算三个体系的能量：蛋白配体复合物的能量、蛋白的能量、配体的能量。

因为没有设置能量检测组，所以edr文件中给出的能量项（比如Disper.corr和Coulomb-SR）都是一整个体系的。比如说生成edr文件时候的轨迹里包含了蛋白配体和水，那得到的这能量项就是蛋白配体和水的体系的总的值，而不是蛋白配体复合物对应的值。即使设置了能量检测组，Disper.corr和Coul.recip.也是整个体系的。故而我们在计算上述三个体系的能量的时候，需要把多余的组分全部去掉。

计算每一个体系的能量的过程大致为：

1. 抽提只包含目标组分的轨迹
2. 生成只包含目标组分的tpr文件
3. 利用只包含目标组分的tpr文件rerun轨迹得到edr文件
4. 利用`gmx energy`命令从edr文件中得到LJ-SR、Coulomb-SR、Disper.corr.和Coul.recip.四项能量。

对三个体系都执行了上面的四步之后，就可以根据相互作用的定义，计算蛋白配体之间的相互作用了。

下面给出我常用的执行这一过程的脚本：

```bash
## shell scripts for energy extract: PME
echo " === extract the prolig energy === "
gmx make_ndx -f md.tpr -o prolig.ndx 
gmx trjconv -f md.xtc -s md.tpr -n prolig.ndx -o prolig.xtc
gmx convert-tpr -s md.tpr -n prolig.ndx -o prolig.tpr
gmx mdrun -s prolig.tpr -rerun prolig.xtc -e prolig.edr 
gmx energy -f prolig.edr -o prolig.xvg

echo " === extract the protein energy === "
gmx make_ndx -f md.tpr -o pro.ndx
gmx trjconv -f md.xtc -s md.tpr -n pro.ndx -o pro.xtc
gmx convert-tpr -s md.tpr -n pro.ndx -o pro.tpr
gmx mdrun -s pro.tpr -rerun pro.xtc -e pro.edr 
gmx energy -f pro.edr -o pro.xvg

echo " === extract the ligand energy === "
gmx make_ndx -f md.tpr -o lig.ndx
gmx trjconv -f md.xtc -s md.tpr -n lig.ndx -o lig.xtc
gmx convert-tpr -s md.tpr -n lig.ndx -o lig.tpr
gmx mdrun -s lig.tpr -rerun lig.xtc -e lig.edr 
gmx energy -f lig.edr -o lig.xvg
```

以蛋白配体复合物的能量抽提为例，详解每一步的作用：

```bash
# 生成蛋白配体复合物的索引组，选择将蛋白配体组合起来成为一个新的组
gmx make_ndx -f md.tpr -o prolig.ndx 
# 从整体轨迹中抽提出复合物的轨迹，选择新生成的蛋白配体复合物
gmx trjconv -f md.xtc -s md.tpr -n prolig.ndx -o prolig.xtc
# 生成只包含复合物的tpr文件，选择复合物组
gmx convert-tpr -s md.tpr -n prolig.ndx -o prolig.tpr
# 利用新的tpr文件rerun新的xtc文件，生成edr
gmx mdrun -s prolig.tpr -rerun prolig.xtc -e prolig.edr 
# 从edr文件中导出4项能量：选择 LJ-SR, Coulomb-SR, Disper.corr., Coul.recip. 
gmx energy -f prolig.edr -o prolig.xvg
```

在执行完上面的脚本之后，会得到3个记录了能量的xvg文件：prolig.xvg、pro.xvg、lig.xvg。每个文件里面包含了四列随模拟时间变化的能量：LJ-SR、Disper.corr.、Coulomb-SR、Coul.recip.。

这里得到的能量都是对应体系的总能量，接下来我们只需要简单计算一下就可以得到蛋白配体之间的能量了。比如说我们要计算蛋白配体之间的LJ-SR，那就用同一模拟时间的复合物的LJ-SR减去蛋白的LJ-SR和配体的LJ-SR即可。LJ-SR加上色散校正项Disper.corr.就是总的范德华力，Coulomb-SR加上Coul.recip.就是总的库仑相互作用。

当然可以使用你熟悉的工具来计算。我因为要算很多次的关系，所以写了脚本来处理这一过程：energy_compute.py (https://github.com/CharlesHahn/Scripts-for-DOCK-and-MD/tree/master/sources/energy_compute)，具体的使用说明可以参见github页面。

上述脚本中用到的md.xtc是主体模拟原始的轨迹，没有经过周期性校正的，如果你rerun的是校正过了的轨迹，得到的能量数据应该也是一摸一样的。



### Method 2

利用截断来计算蛋白配体间的相互作用能量的话，需要修改一下md.mdp文件（就是你用来做主体模拟的输入参数文件）。

在本案例中，原始md.mdp文件中的计算库伦相互作用的方法是PME（coulombtype=PME），这里咱们改成coulombtype=cut-off，还有就是库伦和范德华的截断，咱都改成4.0nm；还需要添加一下蛋白质和配体的能量检测组。

修改之后的这几项：

```md.mdp
energygrps    = Protein Ligands
cutoff-scheme = Verlet
coulombtype   = cut-off
vdwtype       = cut-off
rcoulomb       = 4.0
rvdw          = 4.0
```

GMX不同版本在非键相互作用的设定上有一些差别，上述的设置方式对于我的体系work well。

Method 1中范德华力和库仑力都被分成了短程和长程两项。在截断的方法里，超过截断值的能量，我们就不考虑了。因而我们需要尽可能把截断值设得大一些。对于范德华力来说，原来的1.0nm的截断显然偏小，可能2.0nm差不多吧（我猜的）；对于库伦这种长程作用来说，还需要把截断设置得大一些。对于Verlet截断方案来说，需要将rcoulomb和rvdw设置成一样的值，否则会报错。我们还需要考虑到盒子的大小，**截断值不能超过盒子最短边的一半**，不然就会计算到分子和它的映像分子的作用，这显然不合理。所以我们这里取的截断值为4.0nm，对于本案例来说应该是合适的。

其实这个地方有很多可以思考的。对于蛋白配体体系来说，如果蛋白比较大呢？假设蛋白处于盒子中心且在盒子的每个方向上都占据了一定的长度呢？这时候处于盒子边缘的配体分子，是不是可以在截断半径内同时看到蛋白质的两侧呢？这显然也不太合理，就像两个人见面，我不能同时看见你的前胸和后背，除非你是毕加索画笔下的人物。因而在这种情况下，截断与盒子尺寸的关系应该是：盒子的最短矢量应该超过蛋白质在该方向上的长度**再加上**两倍的截断半径。但是处于节约算力的目的，这个原则常未被遵守。咱们这里也不去细追究了。

总结一下需要修改的地方：

1. 添加蛋白和配体的能量检测组
2. 将coulombtype和vdwtype都设置为cut-off
3. 将库伦和范德华的截断值设置到一个足够大的合适的值

接下来咱们利用修改过后的mdp文件去生成新的tpr文件，然后rerun一下轨迹文件就可以了。

```bash
# 利用修改之后得到的md_cutoff.mdp文件生成新的tpr文件：md_cutoff.tpr
gmx grompp -f md_cutoff.mdp -c npt.gro -r npt.gro -t npt.cpt -p topol.top -o md_cutoff.tpr -maxwarn 1
# 利用新的tpr文件生成复合物索引组
gmx make_ndx -f md_cutoff.tpr -o prolig.ndx 
# 导出复合物的轨迹：选择复合物组
gmx trjconv -f md.xtc -s md_cutoff.tpr -n prolig.ndx -o prolig.xtc
# 生成只包含复合物的tpr文件
gmx convert-tpr -s md_cutoff.tpr -n prolig.ndx -o prolig.tpr
# rerun一下轨迹
gmx mdrun -s prolig.tpr -rerun prolig.xtc -e prolig.edr 
# 生成能量数据：选择LJ-SR和Coulomb-SR
gmx energy -f prolig.edr -o energy_results.xvg
```

我这里还是按照我个人的习惯抽提了一下复合物的轨迹和tpr文件，其实可以直接rerun主体模拟的轨迹，从执行速度上来说，应该差别不大。

直接rerun主体模拟轨迹的话，命令会少很多：

```bash
gmx grompp -f md_cutoff.mdp -c npt.gro -r npt.gro -t npt.cpt -p topol.top -o md_cutoff.tpr -maxwarn 1
gmx mdrun -s md_cutoff.tpr -rerun md.xtc -e prolig.edr 
# 生成能量数据：选择LJ-SR和Coulomb-SR
gmx energy -f prolig.edr -o energy_results.xvg
```

我之前直接rerun主体模拟轨迹，在有的体系里得到的能量数据里有很多奇怪的陡升和陡降，不明觉厉。大家自行测试看看吧。

现在得到的xvg结果文件里面就只包含两列数据了，一列是Coulomb-SR，一列是LJ-SR。如果截断设置得合适，这两列数据其实就是总的库仑相互作用和总的LJ势能了。



### Comparison

我昨儿测了一些数据，我们一起看看：

Table 1. 整个100ns模拟的相互作用能

| 计算方法         | LJ (kJ/mol) | Coulomb (kJ/mol) |
| ------------ | :---------: | :--------------: |
| Method 1     |   -453.7    |      -459.6      |
| cutoff 1.4nm |   -490.5    |      -433.5      |
| cutoff 3.0nm |   -523.8    |      -447.1      |
| cutoff 4.0nm |   -525.2    |      -459.9      |

Table 2. 模拟后20ns的相互作用能

| 计算方法         | LJ (kJ/mol) | Coulomb (kJ/mol) |
| ------------ | :---------: | :--------------: |
| Method 1     |   -488.1    |      -569.5      |
| cutoff 1.4nm |   -531.9    |      -549.5      |
| cutoff 3.0nm |   -569.1    |      -549.9      |
| cutoff 4.0nm |   -570.6    |      -569.7      |

表格中的Method 1指的是利用相互作用定义的方式计算的能量；后面的三个为cutoff方法不同截断下计算得到的能量；LJ和Coulomb都是最终的蛋白配体间的范德华能量和库伦能量。

对于Method 1的计算结果，其LJ计算的时候是1.0nm截断下的LJ-SR加上色散校正项，但是数值明显较截断计算出来的大挺多。我检查了色散项的数据，发现都不大，-20kJ/mol左右，但是这一项在整个模拟过程中几乎没有什么变化。因而或许我们可以这样认为：色散校正项用来代表远程LJ似乎有点儿马虎。因为我们主体模拟时范德华的截断半径设置为1.0nm，故而Method 1计算出来的LJ-SR也不能粗略地代表整体的范德华力。通常的蛋白配体模拟中设置的rvdw大概为1.0nm、1.2nm、1.4nm等。对照下面cutoff方法下不同截断的LJ能量，很容易看出，即使是1.4nm，LJ能量也是没有被充分计算的。而Method 1计算得到的库伦相互作用，也即PME方法下实空间库伦和倒易空间库伦的和，是非常精准的，和cutoff方法下4.0nm截断的数据非常接近。

表中没有列出这两种方法的具体的计算时长，我也没有掐表。但是从感觉上来说，虽然Method 1执行的命令更多，但是似乎速度快很多。而cutoff的方法执行起来就很慢，在我的机器上，Method 2的执行时间似乎都超过了10min（当然这和体系大小也有关系）。速度的差异可能部分来源于有没有使用GPU。Method 1当中因为没有定义能量组，部分的计算是调用GPU执行的；而Method 2则全部在CPU上执行，相较起来应该会慢很多。

**总结一下**：

1. 基于Method 1 计算得到的库伦相互作用应该是非常精准的，当然，前提是主体模拟时候coulombtype设置为PME。
2. 基于Method 1 计算得到的LJ可能有较大误差。这个误差的大小可能取决于主体模拟时候设置的rvdw的数值（前提是vdwtype是设置为cut-off的）。色散项可能也有较大误差。
3. 基于Method 2，如果截断设置合理，那得到的LJ和Coulomb应当都是相当精准的。
4. Method 1的分析速度快（不设能量检测组，可以用GPU），Method 2的速度稍慢。



### Others

上述的数据是基于我手边的单个模拟案例得到的；方法和讨论都是基于我个人对相关材料的理解和分析。诸位如果要使用上述两个方法，一定要基于自己的模拟体系，结合自己的理解分析，审慎辨别地对待。

人是在不断进步的，我想我也如此。回望一年前写的分享，深觉诸多地方尚欠妥当，因而如今越发觉得想要尽善尽美。当然，可能一年后，我又会觉得今天的分享有欠妥帖。故而只能言语中劝诫大家：**勿信一家之言，兼听则明，偏听则可能毕不了业**。

还有一样。做模拟的总免不了各处搜罗来一堆脚本和代码，不管是谁写的，一定自己仔细检查下。懂代码的就一行行读一读，不懂代码的就找个测试案例测下输入输出。人家乐意分享的多半不会是恶意挖坑埋bug，但是难保人不出错，更难保代码不出错，更何况还有个你的案例适不适用的情况。

嗯。

哎本来想着内容没多少，结果写完四千多字。看来我废话确实多，或许适合退学专职去写玄幻小说。

降温啦！钻进小被窝，明年开春儿再起床！