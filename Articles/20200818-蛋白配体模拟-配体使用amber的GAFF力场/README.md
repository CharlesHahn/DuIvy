## 软件与材料

本文内容使用到如下软件：

- Ambertools 20
- GROMACS 2019.5 
- Gaussian 16 

如果使用win10操作系统，可以自己编译前两个软件，也可以找别人编译好的版本。比如[Ambertools](https://liuyujie714.com/63.html)、[GROMACS2019.5](http://sobereva.com/458)；在此表示感谢。Gaussian 16 是要付费购买的，在本文中是用于计算配体的RESP电荷，RESP电荷还可以通过[Multiwfn](http://sobereva.com/441)得到；配体还可以采用bcc电荷(准确性不如RESP)，可以直接用Ambertools计算得到。

本文还用到一个蛋白protein.pdb和一个从drugbank下载的配体小分子ligand.pdb。蛋白将利用GROMACS的力场处理，而配体将用amber的GAFF力场处理。

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

### 对配体应用GAFF力场

打开Ambertools中的tleap，输入如下命令：

```shell 
source leaprc.gaff # 载入力场
loadamberparams ligand.frcmod # 载入配体参数化文件
lig = loadmol2 ligand_resp.mol2 # 载入配体
check lig # 检查配体
saveamberparm lig lig.prmtop lig.inpcrd # 保存amber输入文件
```

### 转换配体结构和拓扑文件到GROMACS格式

接下来我们需要将amber的输入文件转换为GROMACS可用的输入文件，也即top文件和gro文件。

通过Ambertools中的acpype程序进行转换：

```shell
# acpype -h # 输出帮助信息
acpype -p lig.prmtop -x lig.inpcrd
```

如此我们就得到了top文件和gro文件；配体的top文件还需要一点修改才可以被GROMACS使用。



## 构建模拟体系

### 生成蛋白质拓扑文件

我们已经有了蛋白文件protein.pdb，配体的top和gro文件；接下来首先要做的就是创建蛋白质的拓扑文件：

```shell
gmx pdb2gmx -f protein.pdb -o protein.gro -water spc -ignh
# 选择力场
# 这里选择amber99sb力场
```

然后我们在相应目录下就可以看到生成了protein.gro、topol.top、posre.itp。topol.top里面包含了选择的力场的一些参数以及蛋白质的拓扑结构参数等等。

接下来我们需要做一系列文件更改，来把配体加到我们的体系中。

### 拓扑文件修改

先说明本例子中各文件的修改步骤：

1. 复制protein.gro并重命名为complex.gro
2. 将配体的gro文件的坐标部分，复制到complex.gro的坐标部分后面，并相应修改第二行的原子数目
3. 将配体的top文件重命名为ligand.itp
4. 将pdb2gmx命令生成的topol.top文件重命名为protein.itp
5. 新建文件topol.top
6. 将protein.itp中`[ moleculetype ]`部分之前的内容完全剪切到topol.top
7. 将protein.itp中`[ dihedrals ]`部分之后的内容完全剪切到topol.top，包括各种文件的导入部分、[ system ]以及[ molecules ]。
8. 在topol.top中导入蛋白质位置限制文件的前面，导入蛋白质拓扑文件protein.itp
9. 对照ligand.itp中最后的`[ molecules ]`，在topol.top中`[ molecules ]`部分添加配体的名字和数目；在topol.top导入水拓扑文件的前面导入配体拓扑文件ligand.itp。
10. 删去ligand.itp中最后的`[ system ]`和`[ molecules ]`部分
11. 适当修改topol.top中[ system ]的名字（这个可能无所谓）
12. 对于ligand.itp中`[ moleculetype ]`前面的内容，我们也需要整合到topol.top中，一般包含两部分，`[ defaults ]`和`[ atomtypes ]`。对于这两部分内容，直接复制到topol.top文件中导入蛋白质拓扑文件部分的前面即可。
13. 现在我们topol.top文件的最开始，是导入力场参数文件，之后就是配体的`[ defaults ]`部分，由于力场文件中通常也包含这一部分，所以会出错，将`[ defaults ]`部分注释掉；然后尝试对整个体系进行溶剂化、电荷平衡操作，如果不报错，可能就可以了，如果报错，那可能还需要进一步修改topol.top的力场部分。比如用配体的`[ defaults ]`部分置换掉力场的`[ defaults ]`部分等。
14. 有关力场文件的详细信息，可参考[jerkwin的博文](https://github.com/Jerkwin/Jerkwin.github.io/blob/master/GMX/GMXman-5.md)，主要是5.7.1部分；在此表示感谢。

关于本体系中拓扑文件修改的详细情况，说明如下。

首先，在配体和蛋白质的拓扑文件（.top）中，主要包含如下几个部分的内容：

1. [ defaults ]
2. [ atomtypes ]
3. [ moleculetype ]
4. [ atoms ]
5. [ bonds ]
6. [ pairs ]
7. [ angles ]
8. [ dihedrals ]
9. [ system ]
10. [ molecules ]

前两个部分`[ defaults ]` 和`[ atomtypes ]` 主要是相应的力场参数的部分，后面的部分则定义了一些分子拓扑结构相关的参数以及系统参数。一般的itp文件的内容则主要包括3-8的部分，所以我们把ligand.top重命名为ligand.itp、pdb2gmx生成的topol.top重命名为protein.itp之后需要去掉多余的部分；但是去掉的部分并不是不需要的，这部分内容需要与体系的topol.top文件进行整合。

首先将protein.itp的1-2、9-10部分剪切到空白的topol.top中，并include蛋白质的itp文件，得到topol.top如下：

```topol.top
; Include forcefield parameters
#include "amber99sb.ff/forcefield.itp"

; Include protein topology
#include "protein.itp"

; Include Position restraint file
#ifdef POSRES
#include "posre.itp"
#endif

; Include water topology
#include "amber99sb.ff/spc.itp"

#ifdef POSRES_WATER
; Position restraint for each water oxygen
[ position_restraints ]
;  i funct       fcx        fcy        fcz
   1    1       1000       1000       1000
#endif

; Include topology for ions
#include "amber99sb.ff/ions.itp"

[ system ]
; Name
Protein

[ molecules ]
; Compound        #mols
Protein             1
```

然后我们需要将ligand.itp中的内容整合到topol.top中来。

首先是整合`[ system ]`和`[ molecules ]`部分，整合后的topol.top的这两部分如下：

```topol.top
[ system ]
; Name
Protein_and_ligand

[ molecules ]
; Compound        #mols
Protein             1
MOL                 1     
;主要就是添加了配体这一行，这里配体名为MOL
```

整合之后记得删除ligand.itp中的这两部分，并在topol.top导入水拓扑文件的前面导入ligand.itp。

然后整合ligand.itp的前两个部分到topol.top中，直接剪切并粘贴到topol.top导入蛋白质拓扑文件前面、导入力场的后面，然后注释掉`[ default ]`部分。

整合之后的topol.top如下：

```topol.top
; Include forcefield parameters
#include "amber99sb.ff/forcefield.itp"

; [ defaults ]
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
; 1               2               yes             0.5     0.8333

[ atomtypes ]
;name   bond_type     mass     charge   ptype   sigma         epsilon       Amb
 c2       c2          0.00000  0.00000   A     3.39967e-01   3.59824e-01 ; 1.91  0.0860
 c1       c1          0.00000  0.00000   A     3.39967e-01   8.78640e-01 ; 1.91  0.2100
 os       os          0.00000  0.00000   A     3.00001e-01   7.11280e-01 ; 1.68  0.1700
 o        o           0.00000  0.00000   A     2.95992e-01   8.78640e-01 ; 1.66  0.2100

; Include protein topology
#include "protein.itp"

; Include Position restraint file
#ifdef POSRES
#include "posre.itp"
#endif

; Include ligand topology
#include "ligand.itp"

; Include water topology
#include "amber99sb.ff/spc.itp"

#ifdef POSRES_WATER
; Position restraint for each water oxygen
[ position_restraints ]
;  i funct       fcx        fcy        fcz
   1    1       1000       1000       1000
#endif

; Include topology for ions
#include "amber99sb.ff/ions.itp"

[ system ]
; Name
Protein_and_ligand

[ molecules ]
; Compound        #mols
Protein             1
MOL                 1     
```

### 后续建模

即可开始通过GROMACS进行后续的建模，如果在添加离子等后续步骤中不报错的话，体系一般就是ok的。

```gmx
gmx editconf -f complex.gro -o newbox.gro -bt cubic -d 1.0  
gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro  
gmx grompp -f em.mdp -c solv.gro -r solv.gro -p topol.top -o ions.tpr -maxwarn 2  
gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -np 5  
```

### 关于力场

这种体系中，小分子的力场最好选择和蛋白质的力场比较相近的力场；也即同一个体系中，所有的物质，最好都采用同一个力场或者同一类力场。本文中，GAFF力场和amber99sb算是同一类力场吧。

如果采用不同类型的力场（如假设本文中蛋白质采用GROMOS系列的力场），可能会有一些理论上的问题，请咨询相关的资深人员。

### 关于多配体体系

如果蛋白配体体系中需要有多个位置不同的配体，那可以通过pymol等工具，载入ligand_resp.mol2文件，通过复制移动等功能实现蛋白周围的不同位置配体的生成，然后分别保存每个配体为mol2文件。这样就可以得到多个有resp电荷的配体文件，但要注意：pymol等工具导出的文件中原子类型可能和ligand_resp.mol2中定义的**原子类型**不同了，需要对照ligand_resp.mol2文件对原子类型一列进行修改，才能与ligand.frcmod配合进行后续建模。之后再通过Ambertools处理，acpype转换。最后整合和include到体系中应该就可以了。

