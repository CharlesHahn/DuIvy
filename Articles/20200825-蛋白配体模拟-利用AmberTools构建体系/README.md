##  软件和材料

本文内容使用到如下软件：

- Ambertools 20
- GROMACS 2019.5 
- Gaussian 16 

如果使用win10操作系统，可以自己编译前两个软件，也可以找别人编译好的版本。比如[Ambertools](https://liuyujie714.com/63.html)、[GROMACS2019.5](http://sobereva.com/458) 或者 [各种版本GROMACS](https://liuyujie714.com/15.html)；在此表示感谢。Gaussian 16 在本文中用于计算配体的RESP电荷，RESP电荷还可以通过其他方式得到；配体还可以采用bcc电荷(准确性不如RESP)，可以直接用Ambertools计算得到。

本文还用到一个蛋白protein.pdb和一个从drugbank下载的配体小分子ligand.pdb。蛋白将使用amber的ff14SB力场处理，而配体将使用amber的GAFF力场处理。



## 配体处理

### 配体加氢

本次下载得到的ligand.pdb文件里是没有包含非极性氢的，所以需要给它加氢。加氢的方式有很多，比如Ambertools里的reduce、或者pymol等。

```shell
reduce ligand.pdb > ligand_h.pdb
```

reduce 可能更适用于氨基酸组成的分子；不论使用什么软件进行加氢操作，最后一定要检查结果文件中是否包含了加上的氢。还需要检查加氢之后的结构是否合理，不合理的结构无法通过高斯进行计算。

之后通过Ambertools的antechamber处理成Gaussian的输入文件。

```shell
antechamber -i ligand_h.pdb -fi pdb -o ligand.gjf -fo gcrt -pf y
# 当然输入的配体文件也可以是mol2等
# antechamber -i liand_h.mol2 -fi mol2 -o ligand.gjf -fo gcrt -pf y
```

然后利用Gaussian 16 计算：

```shell
g16 < ligand.gjf > ligand.out
```

在Gaussian的输出文件中就有ESP电荷数据；利用antechamber得到包含RESP电荷的mol2文件： 

```shell
antechamber -i ligand.out -fi gout -c resp -o ligand_resp.mol2 -fo mol2 -pf y
```

要在Ambertools中对配体应用GAFF力场，我们还需要生成一个输入文件：

```shell
parmchk -i ligand_resp.mol2 -f mol2 -o ligand.frcmod
```

ligand.frcmod中包含小分子的一些必要参数。

如此，我们所需要的的Ambertools输入文件就准备好了。



## Ambertools建模

### 蛋白去氢

通过reduce去除所有氢原子，tleap在载入蛋白的时候会自动添加所有氢原子的。

```bash
reduce -Trim protein.pdb > protein_removeH.pdb
```

### tleap建模

```bash
# 载入蛋白质和配体的力场
source leaprc.protein.ff14SB
source leaprc.gaff
# 载入配体的参数文件
loadamberparams ligand.frcmod
set default PBRadii mbondi2
# 载入配体
lig = loadmol2 ligand_resp.mol2
check lig
# 载入蛋白质
pro = loadpdb protein_removeH.pdb
check pro
# 将蛋白质和配体组合形成复合物
com = combine {pro lig_num}
check com
# 载入水
source leaprc.water.tip3p
# 对复合物进行溶剂化，盒子边缘距离复合物最小 1 nm
solvatebox com TIP3PBOX 10.0
# 检查体系电荷并平衡电荷
charge com
# 由于本体系带正电，所以添加氯离子知道总体电荷为0
# remind u, 这一步超级慢，建议去学门编程语言再回来
addIons2 com Cl- 0
check com 
# 保存体系的拓扑结构等数据
saveamberparm com com.prmtop com.inpcrd
```

如是，我们便得到了amber跑模拟需要的输入文件，之后利用acpype转换就可以用GROMACS跑了。

Ambertools是一件强大的工具，利用它建模有着更多的选择，譬如amber的力场以及tleap或者xleap中强大的各种命令。



### acpype转换

```bash
acpype -p com.prmtop -x com.inpcrd
```

得到com_GMX.top、com_GMX.gro



## GROMACS运行模拟

利用GROMACS运行模拟需要的mdp文件的模板可以在官方教程网站下载，根据自己的情况修改就成了。

```bash
# 改个名字
cp com_GMX.top topol.top
cp com_GMX.gro complex.gro
# 能量最小化
# 修复IM原子类型无法识别的问题：将 IM 修改为 Cl-
gmx grompp -f em_real.mdp -c complex.gro -r complex.gro -p topol.top -o em.tpr 
gmx mdrun -v -deffnm em
# nvt
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr 
gmx mdrun -deffnm nvt  
# npt
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt 
# md 
gmx grompp -f md.mdp -c npt.gro -r npt.gro -t npt.cpt -p topol.top -o md.tpr  
gmx mdrun -deffnm md 
```

## some tips

### 多配体体系构建

如果构建的体系中需要有多个配体，可以通过pymol等工具利用ligand_resp.mol2生成多个位置不同的配体文件，然后修正原子类型，使其与ligand_resp.mol2中原子类型一致。然后利用这些配体文件，在tleap中载入，构建到体系中就可以了。





