### 体系说明

本文主要针对蛋白配体模拟体系的预处理过程做一个简短的记录。体系中包含一个蛋白质，六个初始位置不同的配体分子。蛋白质使用PDB编号为2BEG的Aβ<sub>17-42</sub> 五聚体，配体为某一多酚类小分子。

本文所使用的力场是根据相应文献和研究进行确定的，实际使用的时候，对于力场以及相关参数的确定需要更加费心加以确定。

### 1 设定初始构象

#### 1.1 蛋白质 PDB 文件获取与处理

首先从PDB下载编号为2BEG的pdb文件，此蛋白质结构是通过NMR方法测定得到。

由于2beg.pdb一共包含了10 个state，故而利用pymol将第一个模型另存为protein_split_from_2beg.pdb。之后需要检查蛋白质结构完整性，利用pymol打开该文件并将展示模式设置为stick，检查N端及C端是否正确；发现C端缺少结束的氧原子，因而需要利用Swiss-PdbViewer添加C端的OXT原子。利用Swiss-PdbViewer打开 protein_split_from_2beg.pdb ，通过Build选项里面的Add C-terminal Oxygen(OXT) 添加OXT原子，之后通过save选项的current layer保存当前的pdb文件为protein_modified_by_SPDBV.pdb 。这时候得到的pdb文件里面包含了很多 SPDBV的信息，这不是GROMACS需要的，需要删除。利用pymol打开protein_modified_by_SPDBV.pdb ，为蛋白质加氢，然后保存为protein.pdb；pymol在保存时会自动去掉没必要的字段。如此便得到了所需的pdb文件。

#### 1.2 配体文件处理
配体的初始文件是来源于ZINC网站的mol2文件；通过pymol打开mol2格式的配体文件，直接导出结构为pdb文件即可。

#### 1.3 初始构象处理

利用pymol同时载入前文提到的protein.pdb和ligand.pdb，利用pymol的分子移动和复制功能，构建蛋白配体的复合物构象。构建的复合物中包含一个蛋白质和6个配体，配体较为均匀和随机地分布在蛋白质周围，且利用pymol的surface模式观察，配体和蛋白质没有表面的接触和重叠。

之后分别以pdb格式单独导出每一个配体，得到6个配体的pdb文件。

### 2 拓扑结构转换

#### 2.1 蛋白质拓扑结构生成

```shell
gmx pdb2gmx -f protein.pdb -o protein_processed.gro -water spc -ignh 
```

上面的命令可以根据蛋白质的pdb文件生成所需的分子坐标文件(gro文件)、拓扑文件(topol.top)、位置限制文件等。`-water spc`指定所用的水的模型；`-ignh`参数会让程序自动忽略输入pdb文件中的氢原子并根据相应的力场自动添加符合力场要求的氢原子；运行上面的命令之后，选择力场 `GROMOS96 53a6`；程序会自动给多肽链添加带正电的N端和带负电的C端；也可以通过`-ter`参数手动确定蛋白质的末端结构。

运行`pymol protein_processed.gro`检查坐标文件。



#### 2.2 配体拓扑结构处理

配体的处理需要用到PRODRG(http://prodrg2.dyndns.org/submit.html)网站，其可以生成配体的相关拓扑文件，使用的应该是GROMOS87力场。

利用该网站依次处理6个配体；首先在网站的输入框内粘贴配体的pdb文件内容的坐标部分，然后Chirality选项选择yes，Charges选项选择Full，EM选项选择No。Chirality意为是否保留分子的手性中心，选择保留；对于电荷而言，选择Full电荷，其使用43A1力场；真空模拟一般使用reduced电荷选项，其使用43B1力场，本研究非真空模拟，况且43B1力场精确性稍差；此步骤中需要保留配体的坐标，故而不进行能量最小化。之后运行就可以得到 prodrg生成的各种文件的压缩包，下载并根据配体编号进行命名，例如ligand_1.tgz 。

之后将压缩包解压，压缩包里的DRGGMX.ITP文件即是配体的itp文件，DRGAPH.GRO即是配体的分子坐标文件；分别将DRGGMX.ITP和DRGAPH.GRO重命名为相应的文件。

使用如下shell命令进行复制：

```shell
cp ligand_1/DRGGMX.ITP ligand_1.itp 
cp ligand_1/DRGAPH.GRO ligand_1.gro 
cp ligand_2/DRGGMX.ITP ligand_2.itp
cp ligand_2/DRGAPH.GRO ligand_2.gro
cp ligand_3/DRGAPH.GRO ligand_3.gro 
cp ligand_3/DRGGMX.ITP ligand_3.itp 
cp ligand_4/DRGAPH.GRO ligand_4.gro 
cp ligand_4/DRGGMX.ITP ligand_4.itp 
cp ligand_5/DRGAPH.GRO ligand_5.gro 
cp ligand_5/DRGGMX.ITP ligand_5.itp 
cp ligand_6/DRGAPH.GRO ligand_6.gro 
cp ligand_6/DRGGMX.ITP ligand_6.itp 
```

之后修改配体的itp文件和gro文件里面的配体命名，防止配体重名。itp文件 [moleculetype] 下的Name改为 \*ZIN，\* 为对应的配体编号，[atoms] 部分的resid也需要对应修改；对应的gro文件里的分子命名也相应改为 \*ZIN 。

可使用如下命令进行修改：

```shell
sed -i "s/ZIN/1ZIN/g" ligand_1.itp 
sed -i "s/ZIN/2ZIN/g" ligand_2.itp  
sed -i "s/ZIN/3ZIN/g" ligand_3.itp  
sed -i "s/ZIN/4ZIN/g" ligand_4.itp  
sed -i "s/ZIN/5ZIN/g" ligand_5.itp  
sed -i "s/ZIN/6ZIN/g" ligand_6.itp  

sed -i "s/1ZIN/1ZIN/g" ligand_1.gro  
sed -i "s/1ZIN/2ZIN/g" ligand_2.gro  
sed -i "s/1ZIN/3ZIN/g" ligand_3.gro  
sed -i "s/1ZIN/4ZIN/g" ligand_4.gro  
sed -i "s/1ZIN/5ZIN/g" ligand_5.gro  
sed -i "s/1ZIN/6ZIN/g" ligand_6.gro  
```

之后检查每个文件命名等是否正确。

### 3 复合物构建

复制protein_processed.gro文件并命名为complex.gro，将每个配体gro文件的坐标部分复制到complex.gro的蛋白质原子坐标末尾，最后更新complex.gro开头的原子总数。

之后需要修改topol.top，添加配体的itp文件以及修改系统组分内容，还需要将蛋白质位置限制所需的文件include进去。

修改之后的topol.top如下所示：

```txt
; Include forcefield parameters
#include "gromos53a6.ff/forcefield.itp"

; Include chain topologies
#include "topol_Protein_chain_A.itp"
#include "topol_Protein_chain_B.itp"
#include "topol_Protein_chain_C.itp"
#include "topol_Protein_chain_D.itp"
#include "topol_Protein_chain_E.itp"

; Include Position restraint file
#include "posre_Protein_chain_A.itp"
#include "posre_Protein_chain_B.itp"
#include "posre_Protein_chain_C.itp"
#include "posre_Protein_chain_D.itp"
#include "posre_Protein_chain_E.itp"

; Include Ligand topologies
#include "ligand_1.itp"
#include "ligand_2.itp"
#include "ligand_3.itp"
#include "ligand_4.itp"
#include "ligand_5.itp"
#include "ligand_6.itp"

; Include water topology
#include "gromos53a6.ff/spc.itp"

#ifdef POSRES_WATER
; Position restraint for each water oxygen
[ position_restraints ]
;  i funct       fcx        fcy        fcz
   1    1       1000       1000       1000
#endif

; Include topology for ions
#include "gromos53a6.ff/ions.itp"

[ system ]
; Name
Protein_Brazilin

[ molecules ]
; Compound        #mols
Protein_chain_A     1
Protein_chain_B     1
Protein_chain_C     1
Protein_chain_D     1
Protein_chain_E     1
1ZIN                1
2ZIN                1
3ZIN                1
4ZIN                1
5ZIN                1
6ZIN                1
```



### 4 周期性边界条件

```shell
# 定义盒子尺寸
gmx editconf -f complex.gro -o newbox.gro -bt cubic -box 8.0 8.0 8.0 
# 添加溶剂
gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro  
# 加入离子
gmx grompp -f em.mdp -c solv.gro -r solv.gro -p topol.top -o ions.tpr -maxwarn 2  
gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -np 5  
```

紧接着利用上面的命令定义周期性的边界条件，设置盒子类型为立方，盒子尺寸为边长8 nm ；之后往盒子里填充水分子。

后面两条命令可以往盒子里添加离子。第一条命令会因为力场(GROMOS96即将被废弃)和体系净电荷问题报错(蛋白质带有五个负电)；使用GROMOS力场对于本体系是合适且已有文献报道的，后文的力场报错也一律忽视；净电荷问题下一条命令即可解决。添加电荷时需要替换掉体系原来的粒子，这里选择替换溶剂(SOL)的分子；因为体系带5个负电，所以这里添加5个钠离子以平衡电荷。如果要添加阴离子，可以使用命令`-nn number`。

em.mdp内容如下：

```mdp
title		= Minimization	; Title of run

; Parameters describing what to do, when to stop and what to save
integrator	= steep		; Algorithm (steep = steepest descent minimization)
emtol		= 1000.0  	; Stop minimization when the maximum force < 10.0 kJ/mol
emstep      = 0.01      ; Energy step size
nsteps		= 50000	  	; Maximum number of (minimization) steps to perform
energygrps	= system	; Which energy group(s) to write to disk

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist		    = 1		    ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet
ns_type		    = grid		; Method to determine neighbor list (simple, grid)
rlist		    = 1.0		; Cut-off for making neighbor list (short range forces)
coulombtype	    = PME		; Treatment of long range electrostatic interactions
rcoulomb	    = 1.0		; long range electrostatic cut-off
rvdw		    = 1.0		; long range Van der Waals cut-off
pbc             = xyz 		; Periodic Boundary Conditions
```



### 5 能量最小化

在前文构建体系的过程中涉及到了一些体系原子的增删等操作，可能使得部分区域的原子之间应力过大，并可能影响最终的成品模拟，所以需要对体系进行弛豫，以减小局部不均衡的作用力。

能量最小化过程中可以进行能量检测，需要生成一个index.ndx文件以方便对蛋白质和配体之间的能量进行检测。

```shell
gmx make_ndx -f solv_ions.gro -o index.ndx  
```

运行命令之后将选择六个配体，将之组合为一个组1ZIN_2ZIN_3ZIN_4ZIN_5ZIN_6ZIN。

之后设置能量最小化的参数文件em_real.mdp如下，相关参数可以根据需要调整

```mdp
title		= Minimization	; Title of run

; -----------------------------------------------------------
; Parameters describing what to do, when to stop and what to save
integrator	= steep		; Algorithm (steep = steepest descent minimization)
emtol		= 500.0  	; Stop minimization when the maximum force < 5.0 kJ/mol
emstep      = 0.01      ; Energy step size
nsteps		= 50000	  	; Maximum number of (minimization) steps to perform
energygrps	= Protein 1ZIN_2ZIN_3ZIN_4ZIN_5ZIN_6ZIN	; Which energy group(s) to write to disk

; -----------------------------------------------------------
; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist		    = 1		    ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet
ns_type		    = grid		; Method to determine neighbor list (simple, grid)
rlist		    = 1.0		; Cut-off for making neighbor list (short range forces)
coulombtype	    = PME		; Treatment of long range electrostatic interactions
rcoulomb	    = 1.0		; long range electrostatic cut-off
rvdw		    = 1.0		; long range Van der Waals cut-off
pbc		        = xyz 		; Periodic Boundary Conditions
```

之后运行命令生成能量最小化所需的em.tpr：

```shell
gmx grompp -f em_real.mdp -c solv_ions.gro -r solv_ions.gro -p topol.top -n index.ndx -o em.tpr -maxwarn 2  
```

得到em.tpr文件之后通过以下命令运行能量最小化：

```shell
gmx mdrun -v -deffnm em  
```

在体系能量收敛之后，可以通过`gmx energy`抽提出能量最小化过程的能量变化并检验。



### 6 NVT预平衡

在主体模拟之前，还需要进行NVT预平衡和NPT预平衡，以使得体系稳定在合适的温度和压力下。

NVT预平衡的时长设定为100 ps，温度设定为310 K，这里的温度耦合组沿用传统设置，即设置为蛋白质和非蛋白质。在预平衡过程中需要对蛋白质施加位置限制，有益于平衡充分。在 mdp 文件中使用`define  = -DPOSRES`对蛋白质进行位置限制。另外在NVT预平衡时需要对体系粒子赋予初始速度，而在后面的模拟过程中则不需要再赋予速度而是直接沿用前面模拟的数据。

nvt.mdp文件如下：

```mdp
title       = Protein-ligand complex NVT equilibration 
define      = -DPOSRES  ; position restrain the protein

; Run parameters
; ----------------------------------------
integrator  = md        ; leap-frog integrator
nsteps      = 50000     ; 0.002 * 50000 = 100 ps
dt          = 0.002     ; 2 fs

; Output control
; ----------------------------------------
nstxout     = 500       ; save coordinates every 1.0 ps
nstvout     = 500       ; save velocities every 1.0 ps
nstenergy   = 500       ; save energies every 1.0 ps
nstlog      = 500       ; update log file every 1.0 ps
energygrps  = Protein 1ZIN_2ZIN_3ZIN_4ZIN_5ZIN_6ZIN

; Bond parameters
; ----------------------------------------
continuation    = no            ; first dynamics run
constraint_algorithm = lincs    ; holonomic constraints 
constraints     = all-bonds     ; all bonds (even heavy atom-H bonds) constrained
lincs_iter      = 1             ; accuracy of LINCS
lincs_order     = 4             ; also related to accuracy

; Neighborsearching
; ----------------------------------------
cutoff-scheme   = Verlet
ns_type         = grid      ; search neighboring grid cells
nstlist         = 10        ; 20 fs, largely irrelevant with Verlet
rcoulomb        = 1.4       ; short-range electrostatic cutoff (in nm)
rvdw            = 1.4       ; short-range van der Waals cutoff (in nm)

; Electrostatics
; ----------------------------------------
coulombtype     = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order       = 4         ; cubic interpolation
fourierspacing  = 0.16      ; grid spacing for FFT

; Temperature coupling
; ----------------------------------------
tcoupl      = V-rescale                     ; modified Berendsen thermostat
tc-grps     = Protein Non-Protein           ; two coupling groups - more accurate
tau_t       = 0.1   0.1                     ; time constant, in ps
ref_t       = 310   310                     ; reference temperature, one for each group, in K

; Pressure coupling
; ----------------------------------------
pcoupl      = no        ; no pressure coupling in NVT

; Periodic boundary conditions
; ----------------------------------------
pbc         = xyz       ; 3-D PBC

; Dispersion correction
; ----------------------------------------
DispCorr    = EnerPres  ; account for cut-off vdW scheme

; Velocity generation
; ----------------------------------------
gen_vel     = yes       ; assign velocities from Maxwell distribution
gen_temp    = 310       ; temperature for Maxwell distribution
gen_seed    = -1        ; generate a random seed
```

之后运行如下命令开始NVT预平衡：

```shell
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr -maxwarn 2  
gmx mdrun -deffnm nvt  
```

检查NVT预平衡的结果，注意温度是否达到稳定。



### 7 NPT预平衡

NPT预平衡的时候需要产生压力，本研究使用Berendsen 法进行控压，设定压力参考值为1 bar。NPT预平衡的时长也设置为100 ps。

npt.mdp文件如下：

```mdp
title       = Protein-ligand complex NPT equilibration 
define      = -DPOSRES  ; position restrain the protein and ligand

; Run parameters
;-----------------------------------------
integrator  = md        ; leap-frog integrator
nsteps      = 50000     ; 2 * 50000 = 100 ps
dt          = 0.002     ; 2 fs

; Output control
;-----------------------------------------
nstxout     = 500       ; save coordinates every 1.0 ps
nstvout     = 500       ; save velocities every 1.0 ps
nstenergy   = 500       ; save energies every 1.0 ps
nstlog      = 500       ; update log file every 1.0 ps
energygrps  = Protein 1ZIN_2ZIN_3ZIN_4ZIN_5ZIN_6ZIN

; Bond parameters
;-----------------------------------------
continuation    = yes           ; first dynamics run
constraint_algorithm = lincs    ; holonomic constraints 
constraints     = all-bonds     ; all bonds (even heavy atom-H bonds) constrained
lincs_iter      = 1             ; accuracy of LINCS
lincs_order     = 4             ; also related to accuracy

; Neighborsearching
;-----------------------------------------
cutoff-scheme   = Verlet
ns_type         = grid      ; search neighboring grid cells
nstlist         = 10        ; 20 fs, largely irrelevant with Verlet
rcoulomb        = 1.4       ; short-range electrostatic cutoff (in nm)
rvdw            = 1.4       ; short-range van der Waals cutoff (in nm)

; Electrostatics
;-----------------------------------------
coulombtype     = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order       = 4         ; cubic interpolation
fourierspacing  = 0.16      ; grid spacing for FFT

; Temperature coupling
;-----------------------------------------
tcoupl      = V-rescale                     ; modified Berendsen thermostat
tc-grps     = Protein Non-Protein           ; two coupling groups - more accurate
tau_t       = 0.1   0.1                     ; time constant, in ps
ref_t       = 310   310                     ; reference temperature, one for each group, in K

; Pressure coupling
;-----------------------------------------
pcoupl      = Berendsen   ; Parrinello-Rahman   ; pressure coupling is on for NPT
pcoupltype  = isotropic                     ; uniform scaling of box vectors
tau_p       = 2.0                           ; time constant, in ps
ref_p       = 1.0                           ; reference pressure, in bar
compressibility = 4.5e-5                    ; isothermal compressibility of water, bar^-1
refcoord_scaling    = com

; Periodic boundary conditions
;-----------------------------------------
pbc         = xyz       ; 3-D PBC

; Dispersion correction
;-----------------------------------------
DispCorr    = EnerPres  ; account for cut-off vdW scheme

; Velocity generation
;-----------------------------------------
gen_vel     = no        ; velocity generation off after NVT 
```

运行如下命令进行 NPT 预平衡：

```shell
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr -maxwarn 2  
gmx mdrun -deffnm npt  
```

检查NPT预平衡的结果。由于GROMACS控压时压力通常波动很大，故而压力的瞬时值不具有参考价值，可以关注压力的平均值。通常通过体系的密度来判断NPT预平衡是否充分。



### 8 预平衡检查

能量最小化、NVT预平衡和NPT预平衡检查的方式类似，可以通过监测体系参数变化以及配体蛋白的动态轨迹来进行。

体系参数监测的方法如下：

```shell 
gmx energy -f xxx.edr -o xxx_energy.xvg
```

上述语句可以从xxx.edr文件中抽提出相关的能量以及体系参数数据，xxx.edr文件是GROMACS程序产生的系统参数的二进制文件，里面包含了体系多种参数随时间变化的数据，如果在mdp控制文件中有设置能量组，则此文件中还包含能量组的相关势能数据。运行上述命令之后，选择需要输出到xvg文件的参数名，即可输出相应的参数变化数据并自动给出参数的平均值等统计数据。将xvg文件中的数据可视化之后即可直观了解参数随模拟时间的变化。对于较为合适的能量最小化过程而言，体系的potential需要在指定的最多步数以内收敛到指定值以下。对于NVT预平衡而言，最重要的是需要检查体系的温度是否达到设定的310 K并基本维持在310 K附近。对于NPT预平衡而言，则需要注意体系的密度或者盒子的体积是否达到稳定，也需要关注温度是否恒定；因为控压算法的关系，体系中的压力波动往往较大，故而压力的瞬时值会较之设定值有着较大的偏移，一般只需要保证体系的密度较为稳定即可。

另外还可以通过pymol监测蛋白配体的动态轨迹变化，查看有无异常：

```shell
gmx trjconv -f xxx.trr -o xxx.xtc
```

上述命令可以将trr格式的轨迹文件转变为pymol可读的xtc格式轨迹文件。之后运行pymol并输入如下指令：

```pymol
load xxx_theFormer.gro, xxx
load xxx.xtc, xxx
```

xxx_theFormer.gro表示本次模拟之前的gro文件名，如要载入NPT过程的轨迹，这里就需要输入nvt.gro，也即NVT过程结束的gro文件，逗号后面的xxx为标记，两次load的标记需要一致。如此即可在pymol中可视化观察模拟中蛋白配体的位置变化，通过转换显示模式，还可以显示出配体和蛋白之间氢键的动态形成与消失。由于在NVT 和NPT过程中对蛋白质施加了位置限制，所以需要在模拟结束之后检查蛋白质的位置以及结构有无较大的变化，还需要注意蛋白质周围是否被水分子合理填充等。

预平衡的过程不是必须的，但是对于主体模拟的采样和稳定有着很重要的作用。



### 9 成品模拟

成品模拟的时候不需要再进行蛋白质位置限制，手动注释掉topol.top中蛋白质位置限制的命令即可取消蛋白质位置限制。

成品模拟的时间设定为100 ns，采样间隔为10 ps；这里取消了能量监测，可以节省算力，后面得到轨迹之后可以重跑得到能量；控压的算法更改为Parrinello-Rahman 。

md.mdp设置如下：

```mdp
title       = Protein-ligand complex MD simulation 

; Run parameters
;-----------------------------------------
integrator  = md          ; leap-frog integrator
nsteps      = 50000000    ; 100ns
dt          = 0.002       ; 2 fs

; Output control
;-----------------------------------------
nstxout             = 0         ; suppress .trr output 
nstvout             = 0         ; suppress .trr output
nstenergy           = 5000      ; save energies every 10.0 ps
nstlog              = 5000      ; update log file every 10.0 ps
nstxout-compressed  = 5000      ; write .xtc trajectory every 10.0 ps
compressed-x-grps   = System
; energygrps          = Protein 1ZIN_2ZIN_3ZIN_4ZIN_5ZIN_6ZIN

; Bond parameters
;-----------------------------------------
continuation    = yes           ; first dynamics run
constraint_algorithm = lincs    ; holonomic constraints 
constraints     = all-bonds     ; all bonds (even heavy atom-H bonds) constrained
lincs_iter      = 1             ; accuracy of LINCS
lincs_order     = 4             ; also related to accuracy

; Neighborsearching
;-----------------------------------------
cutoff-scheme   = Verlet
ns_type         = grid      ; search neighboring grid cells
nstlist         = 10        ; 20 fs, largely irrelevant with Verlet
rcoulomb        = 1.4       ; short-range electrostatic cutoff (in nm)
rvdw            = 1.4       ; short-range van der Waals cutoff (in nm)

; Electrostatics
;-----------------------------------------
coulombtype     = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order       = 4         ; cubic interpolation
fourierspacing  = 0.16      ; grid spacing for FFT

; Temperature coupling
;-----------------------------------------
tcoupl      = V-rescale                     ; modified Berendsen thermostat
tc-grps     = Protein Non-Protein    ; two coupling groups - more accurate
tau_t       = 0.1   0.1                     ; time constant, in ps
ref_t       = 310   310                     ; reference temperature, one for each group, in K

; Pressure coupling 
;-----------------------------------------
pcoupl      = Parrinello-Rahman             ; pressure coupling is on for NPT
pcoupltype  = isotropic                     ; uniform scaling of box vectors
tau_p       = 2.0                           ; time constant, in ps
ref_p       = 1.0                           ; reference pressure, in bar
compressibility = 4.5e-5                    ; isothermal compressibility of water, bar^-1

; Periodic boundary conditions
;-----------------------------------------
pbc         = xyz       ; 3-D PBC

; Dispersion correction
;-----------------------------------------
DispCorr    = EnerPres  ; account for cut-off vdW scheme

; Velocity generation
;----------------------------------------
gen_vel     = no        ; assign velocities from Maxwell distribution
```

首先运行如下命令生成 tpr 文件：

```shell
gmx grompp -f md.mdp -c npt.gro -r npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_results_1.tpr -maxwarn 2  
```

运行主体模拟：

```shell
gmx mdrun -deffnm md_results_1 
```

如果所使用的GROMACS支持GPU加速，且前面没有设置能量监测组，那么可以在命令后面加上`-nb gpu`来使用GPU加速运算。


如果使用有PBS作业系统的集群或者服务器进行运算，那么可以通过`qsub`命令提交如下PBS命令文件来运行模拟：

```pbs
#PBS -N MD_ProLig
#PBS -l nodes=1:ppn=16
#PBS -o output.txt
#PBS -j oe
#PBS -q low

gmx mdrun -deffnm md_results_1
```

通过PBS命令定时检查运行状态，模拟结束之后下载数据进行分析。  



### Others

学习过程中多次参考Jerkwin的相关博客和其翻译的GROMACS tutorial，在此致谢。

本文仅是我分子动力学模拟的相关学习整理，以蛋白配体的某一模拟来说明MD的预处理过程，通用性尚未考虑。数据分析部分在下一文中进行说明。