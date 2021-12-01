# “混用力场参数”

目前主要使用amber和gromacs这两个模拟软件（组件套装），前面已经记录了**如何配体使用amber的gaff力场**，以及**如何蛋白配体都使用amber的力场**，这两种方法（流程）都有各自的小问题，比如蛋白质还可能还是使用amber的ff14SB更精准，而整个建模过程都在ambertools中进行的话…… addIons2的速度想必诸位是知道的；跑一次这个命令，可以出门看一部电影了orz。

有什么办法可以绕开addIons2呢？既想蛋白配体离子水都使用amber的力场，又不想用这个命令。

本文将介绍**如何利用tleap对蛋白配体应用amber的力场，并利用gromacs添加离子和溶剂**。

本文可能涉及到不同力场（当然是同类力场）的混合，并且建模过程可能不常规，**实用价值不强**，建议**当作消遣来看**；如是正经研究，请一定再详细咨询相关专业人员（不包括我，我小白）。


## 材料与软件准备

本文会用到一个蛋白质pdb文件protein.pdb，以及已经计算好resp电荷的配体文件（lig.mol2和ligand.frcmod）。

使用到的程序是ambertools20（主要是tleap和acpype）以及gromacs2019。

相关的部分内容之前发文过了，本文中可能略述。


## 对蛋白应用ff14SB力场

使用tleap对蛋白质应用力场并保存为pro.prmtop, pro.inpcrd。

如果你的蛋白质里面有**半胱氨酸**或者**组氨酸**的话，一定记得根据其状态修改残基名字。可参考[tutorial](https://ambermd.org/tutorials/advanced/tutorial1_adv/index.htm)，关键段落抄录如下：Since pdb files do not distinguish between cysteine residues that are involved in bonds or other things (and hence which have no proton on the sulphur atom) we need to edit the relevant cysteine residues to correct for this. The residue name used by Leap for regular protonated cysteine residues is CYS, for deprotonated (-ve charge) it is CYM and for cysteine residues involved in disulphide bridges and other bonds it is CYX. The same is true for histidine residues which can be protonated in the 'delta' position (HID), the epsilon position (HIE) or at both (HIP).

```tleap
source leaprc.protein.ff14SB
pro = loadpdb protein.pdb
check pro
saveamberparm pro pro.prmtop pro.inpcrd
quit
```

之后使用acpype转换为gmx可用的top和gro文件格式。

```acpype
acpype -p pro.prmtop -x pro.inpcrd
```

如此就得到了只包含蛋白的坐标文件pro.gro和拓扑文件pro.top。


## 对配体应用gaff力场

跟前一步的步骤一样，先tleap处理。

```tleap
source leaprc.gaff
loadamberparams ligand.frcmod
lig = loadmol2 lig.mol2
check lig
saveamberparm lig lig.prmtop lig.inpcrd
quit
```

然后acpype转换

```acpype
acpype -p lig.prmtop -x lig.inpcrd
```

得到只包含配体的坐标文件lig.gro和拓扑文件lig.top。


## 组合力场

这一步是最重要的，我们需要把蛋白质、配体、离子和水的力场和拓扑文件组合起来。

我们得到的pro.top和lig.top文件的最前面，应该就是力场参数，后面是一大堆拓扑的信息，比如这里lig.top文件的结构是这样的(......表示省略很多行数据)：

```lig.top
[ defaults ]
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
1               2               yes             0.5     0.8333

[ atomtypes ]
;name   bond_type     mass     charge   ptype   sigma         epsilon       Amb
 c2       c2          0.00000  0.00000   A     3.39967e-01   3.59824e-01 ; 1.91  0.0860
    ......

[ moleculetype ]
;name            nrexcl
 lig               3

[ atoms ]
;   nr  type  resi  res  atom  cgnr     charge      mass       ; qtot   bond_type
     1   c2     1   MOL    C1    1     0.100372     12.01000 ; qtot 0.100
    ......

[ bonds ]
;   ai     aj funct   r             k
     1      2   1    1.3343e-01    4.7647e+05 ;     C1 - C2    
    ......

[ pairs ]
;   ai     aj    funct
     1      7      1 ;     C1 - C7    
    ......

[ angles ]
;   ai     aj     ak    funct   theta         cth
     1      2      3      1    1.2181e+02    5.7990e+02 ;     C1 - C2     - C3    
    ......

[ dihedrals ] ; propers
; for gromacs 4.5 or higher, using funct 9
;    i      j      k      l   func   phase     kd      pn
     1      2      3      4      9   180.00  27.82360   2 ;     C1-    C2-    C3-    C4
    ......

[ dihedrals ] ; impropers
; treated as propers in GROMACS to use correct AMBER analytical function
;    i      j      k      l   func   phase     kd      pn
     1      3      2     22      4   180.00   4.60240   2 ;     C1-    C3-    C2-    O7
    ......

[ system ]
 lig

[ molecules ]
; Compound        nmols
 lig              1     
 ```

那么很清晰我们可以看到它主要包含了力场、拓扑、系统三个方面的信息：

- **force field**
    + defaults
    + atomtypes
- **topol**
    + moleculetype
    + atoms
    + bonds
    + pairs
    + angles
    + dihedrals (propers)
    + dihedrals (impropers)
- **system**
    + system
    + molecules


在了解了这些之后，我们就可以开始着手整合啦，首先我们来整合蛋白质和配体的文件，步骤如下：

- 整合两者的gro文件
    1. 复制pro.gro到comlex.gro
    2. 复制lig.gro的坐标部分到comlex.gro的坐标部分后面
    3. 记得更改comlex.gro开头的原子数目了
- 整合top文件
    1. 新建一个文件forcefield.itp，这里面存放所有的力场相关参数
        1.1 将蛋白质top文件的defaults和atomtypes两部分复制到forcefield.itp中
        1.2 配体top文件中的defaults部分应该是和蛋白质一样的，那就只需要复制配体top文件的atomtypes部分到forcefield.top的atomtypes部分后面就可以了
    2. 新建一个pro.itp来存放蛋白质的拓扑数据，将蛋白质的拓扑部分从pro.top一股脑复制过去就成
    3. 新建一个lig.itp来存放配体的拓扑数据，将配体的拓扑部分从lig.top中一股脑复制过去
    4. 新建一个topol.top作为总的拓扑文件，在里面导入各物质的itp，还需要在里面整合pro.top和lig.top的系统部分，系统名字无所谓，molecules那里加一下就好啦

整理好的topol.top大致像这个样子：

```topol.top
; include force field
#include "forcefield.itp"

; include protein topology
#include "pro.itp"

; include ligand topology
#include "lig.itp"


; Include water topology


; Include topology for ions


[ system ]
 ProLig

[ molecules ]
; Compound        nmols
 pro              1     
 lig              1     
```

要注意到，用gmx建模，一般第一步`gmx pdb2gmx`的时候，我们就指定了力场和要选用的水模型，这时候生成的topol.top里面就自动导入了水和离子的itp文件。而我们是没有这一步的，因而在上面的文件中我空出来了水和离子的itp导入位置。

那么接下来我们就来尝试导入相应的离子和水的itp文件。

如果你不介意蛋白使用ff14SB力场、配体使用gaff力场、水和离子使用gromacs内置的amber99sb力场的话，请看method one；如果还是想水和离子的itp都来自于amber的话，请看method two。

#### method one

gromacs2019里面是有内置amber99sb力场的，比起ff14SB来说，肯定是老了，但是至少同类力场，各种参数的定义也差不离，混混影响不大。

amber99sb力场的水和离子的itp在哪里呢？一般是在gromacs安装目录下的share/gromacs/top这个路径下，里面有个叫amber99sb.ff的文件夹，打开就能看到ions.itp和tip3p.itp（spc.itp也在这里），直接复制到我们的工作路径。

当然，事情不可能这么简单（柯南推眼镜.jpg)。第一，水和离子的力场参数还没有，我们只是拿到了拓扑文件；第二，ions.itp里面有很多种离子，我们真的需要那么多嘛？一般就需要Na+/Cl-吧？成，咱删掉多余的，只保留Na+/Cl-两种。

至于水和离子的力场参数，我们可以根据ions.itp和tip3p.itp中的atomtypes去找；打开amber99sb.ff文件夹下的forcefield.itp看看，发现它导入了ffbonded.itp和ffnonbonded.itp，依次打开这俩文件找找。

可以在ffnonbonded.itp中发现如下字段：

```ffnonbonded.itp
[ atomtypes ]
; name      at.num  mass     charge ptype  sigma      epsilon
......
Cl          17      35.45    0.0000  A   4.40104e-01  4.18400e-01
Na          11      22.99    0.0000  A   3.32840e-01  1.15897e-02
......
HW           1       1.008   0.0000  A   0.00000e+00  0.00000e+00
......
OW           8      16.00    0.0000  A   3.15061e-01  6.36386e-01
......
```

那么我们先把这部分抄到我们工作路径下的forcefield.itp的atomtypes下面。

然后在ffbonded.itp中发现如下字段：

```ffbonded.itp
[ bondtypes ]
; i    j  func       b0          kb
  OW HW         1    0.09572   462750.4 ; P water
  HW HW         1    0.15136   462750.4 ; P water
......

[ angletypes ]
;  i    j    k  func       th0       cth
HW  OW  HW           1   104.520    836.800 ; TIP3P water
HW  HW  OW           1   127.740      0.000 ; (found in crystallographic water with 3 bonds)
......
```

其实像角度等，在tip3p.itp文件中也有定义的，但咱还是把相关的都抄进咱的forcefield.itp中。

最后，我们需要在topol.top里导入tip3o.itp和ions.itp。

#### method two

如果你之前用tleap做过某些类似体系（相同力场），用到了离子和水模型，那么那些体系里面就有水和离子的力场、拓扑参数。

可能像这样：

```ff
[ atomtypes ]
;name   bond_type     mass     charge   ptype   sigma         epsilon       Amb
 ......
 Cl-      Cl-         0.00000  0.00000   A     4.47766e-01   1.48913e-01 ; 2.51  0.0356
 OW       OW          0.00000  0.00000   A     3.15075e-01   6.35968e-01 ; 1.77  0.1520
 HW       HW          0.00000  0.00000   A     0.00000e+00   0.00000e+00 ; 0.00  0.0000
 ......
 ......
 ......

[ moleculetype ]
  ; molname       nrexcl
  CL-             1

[ atoms ]
  ; id_    at type res nr  residu name     at name  cg nr  charge   mass
    1       Cl-      1         CL-           CL-      1     -1     35.45300

[ moleculetype ]
; molname       nrexcl ; TIP3P model
  WAT             2

[ atoms ]
;   nr   type  resnr residue  atom   cgnr     charge       mass
     1     OW      1     WAT     O      1     -0.834   16.00000
     2     HW      1     WAT    H1      1      0.417    1.00800
     3     HW      1     WAT    H2      1      0.417    1.00800

#ifdef FLEXIBLE
[ bonds ]
; i j   funct   length  force.c.
1   2   1   0.09572   462750.4 0.09572   462750.4
1   3   1   0.09572   462750.4 0.09572   462750.4

[ angles ]
; i j   k   funct   angle   force.c.
2   1   3   1   104.520    836.800  104.520    836.800
#else
[ settles ]
; i j   funct   length
1   1   0.09572 0.15139

[ exclusions ]
1   2   3
2   1   3
3   1   2
#endif
```

这就很好，可以直接抄。咱直接把atomtypes部分抄到工作目录下的forcefield.itg文件中，把下面的离子的moleculetype和atoms抄到ions.itp中，把剩下的水的部分抄到tip3p.itp。

首先对离子的部分进行修改，把forcefield.itp中的两个Cl-**去掉负号**。然后把ions.tpr中的所有CL-或者Cl-**去掉负号**。

这里还需要对tip3p.itp进行一下修改，把moleculetype下的molname以及atoms下面的residue的名字改成**SOL**，因为gmx里面水叫这名儿，然后把atoms部分下的atom列的O、H1、H2依次改为**OW、HW1、HW2**，因为我们之后溶剂化时使用的spc216.gro里面的原子是这样的名字。要保证所有的名儿能对应起来。

抄好之后把水和离子的itp导入到topol.top中就行了。

这种方法可能更好，当然差别估计也不大。

#### 后续步骤

```gmx
gmx editconf -f complex.gro -o newbox.gro -bt cubic -d 1.0  
gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro  
gmx grompp -f em.mdp -c solv.gro -r solv.gro -p topol.top -o ions.tpr -maxwarn 1  
gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -nname CL -pname NA -neutral
# EM
gmx grompp -f em_real.mdp -c solv_ions.gro -r solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em  
```

如果一直到EM都没有出错，那前面的整个过程应该就是ok的。

这里注意到填充水盒子的时候使用的是spc216.gro，因为gromacs2019中没有tip3p.gro，spc也是三点电荷的水模型，gro里面也只记录了原子坐标，而其拓扑结构是由tip3p.itp决定的，所以这里是ok的。

#### Others

> 如果需要对蛋白或者某种物质进行位置限制该怎么办？

如果看过某个gmx的位置限制itp文件你就会知道，它格式很简单，举例如下：

```porse.itp
; In this topology include file, you will find position restraint
; entries for all the heavy atoms in your original pdb file.
; This means that all the protons which were added by pdb2gmx are
; not restrained.

[ position_restraints ]
; atom  type      fx      fy      fz
     1     1  1000  1000  1000
     5     1  1000  1000  1000
     7     1  1000  1000  1000
    10     1  1000  1000  1000
    13     1  1000  1000  1000
    14     1  1000  1000  1000
......
```

后面四列内容是不变的，第一列是要位置限制的物质的重原子序号。这完全可以编个简单的脚本来实现。

如果你不介意限制所有原子的话，那就更简单了，在vim里面两行python就可以搞定。

之后把限制文件include到topol.top就可以啦。

> 这样吭哧吭哧搞了一大堆，有啥优势嘛？

没。。。

正经要发文章的还是正经的来。

如果自己探索的话，这样的方式可以让你对力场、软件等更熟悉；同时你用过的tip3p.itp、ions.itp等也可以下一次用，那就会节省不少时间；听说tleap要控制溶液盐浓度好像挺不方便的，这样做或许能够比较容易利用gmx genion来设定盐浓度吧。

本文成文仓促，如有错漏，烦请不吝赐教和批评。


最后祝大家天天有灵感，并且灵感没有几年前就被人做过了吧（笔者最近苦恼这个哈哈哈）。









