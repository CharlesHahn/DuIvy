# generate the protein.gro
## 沿用带电荷的多肽末端
gmx pdb2gmx -f protein.pdb -o protein_processed.gro -water spc -ignh

## prodrg，to get the zip file
## 参数 ：  Chirality 选项选择 yes，Charges 选项选择 Full，EM选项选择 No

# ligand_cp.sh
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
# sed_ligand_itp  
sed -i "s/ZIN/1ZIN/g" ligand_1.itp 
sed -i "s/ZIN/2ZIN/g" ligand_2.itp  
sed -i "s/ZIN/3ZIN/g" ligand_3.itp  
sed -i "s/ZIN/4ZIN/g" ligand_4.itp  
sed -i "s/ZIN/5ZIN/g" ligand_5.itp  
sed -i "s/ZIN/6ZIN/g" ligand_6.itp  
# sed_ligand_gro 
sed -i "s/1ZIN/1ZIN/g" ligand_1.gro  
sed -i "s/1ZIN/2ZIN/g" ligand_2.gro  
sed -i "s/1ZIN/3ZIN/g" ligand_3.gro  
sed -i "s/1ZIN/4ZIN/g" ligand_4.gro  
sed -i "s/1ZIN/5ZIN/g" ligand_5.gro  
sed -i "s/1ZIN/6ZIN/g" ligand_6.gro  

cp protein_processed.gro complex.gro  

# cp the content of ligand.gro to complex.gro

# modify the topol.top : 
## 1. 添加 ligand itp 
## 2. ligand 名字
## 3. 蛋白质位置限制


# gmx 
# gmx editconf -f complex.gro -o newbox.gro -bt cubic -d 1.0  
gmx editconf -f complex.gro -o newbox.gro -bt cubic -box 8.0 8.0 8.0 
gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro  
gmx grompp -f em.mdp -c solv.gro -r solv.gro -p topol.top -o ions.tpr -maxwarn 2  
gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -np 5  
# gmx 对protein和ligands 进行能量最小化 
gmx make_ndx -f solv_ions.gro -o index.ndx  
gmx grompp -f em_real.mdp -c solv_ions.gro -r solv_ions.gro -p topol.top -n index.ndx -o em.tpr -maxwarn 2  
gmx mdrun -v -deffnm em  
# nvt
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr -maxwarn 2  
gmx mdrun -deffnm nvt  
# npt
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr -maxwarn 2  
gmx mdrun -deffnm npt  
# md 
gmx grompp -f md.mdp -c npt.gro -r npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_results_1.tpr -maxwarn 2  
gmx mdrun -deffnm md_results_1 -nb gpu  

# energy extract  
gmx trjconv -f md_results_1.xtc -s md_results_1.tpr -n -o lig.xtc  
# 下面这句不是必须的
# gmx trjconv -f md_results_1.gro -s md_results_1.tpr -n -o lig.gro    
gmx convert-tpr -s md_results_1.tpr -n  -o lig.tpr  
gmx mdrun -s lig.tpr -rerun lig.xtc -e lig.edr  


## PROTEIN LIGAND ENERGY 
gmx make_ndx -f npt.gro -o prolig.ndx
gmx trjconv -f md_results_1.xtc -s md_results_1.tpr -n prolig.ndx -o prolig.xtc
gmx convert-tpr -s md_results_1.tpr -n prolig.ndx -o prolig.tpr
gmx mdrun -s prolig.tpr -rerun prolig.xtc -e prolig.edr 

gmx make_ndx -f npt.gro -o pro.ndx
gmx trjconv -f md_results_1.xtc -s md_results_1.tpr -n pro.ndx -o pro.xtc
gmx convert-tpr -s md_results_1.tpr -n pro.ndx -o pro.tpr
gmx mdrun -s pro.tpr -rerun pro.xtc -e pro.edr 

gmx make_ndx -f npt.gro -o lig.ndx
gmx trjconv -f md_results_1.xtc -s md_results_1.tpr -n lig.ndx -o lig.xtc
gmx convert-tpr -s md_results_1.tpr -n lig.ndx -o lig.tpr
gmx mdrun -s lig.tpr -rerun lig.xtc -e lig.edr 



# 校验轨迹完整性
gmx make_ndx -f npt.gro -o prolig.ndx
gmx trjconv -s npt.tpr -f md_results_1.xtc -n prolig.ndx -o md_results_1_whole.xtc -pbc whole
# 手动创造一个只包含 蛋白质和配体 的 pdb 文件


########## 续跑
# 新建一个后继 run
gmx grompp -f md_2.mdp -c md_results_1.gro -r md_results_1.gro -t md_results _1.cpt -p topol.top -n index.ndx -o md_results_2.tpr -maxwarn 2
gmx mdrun -deffnm md_results_2

# 直接续跑
# 1. 更改原始md.mdp 的 nsteps, 生成新的 md_results_1.tpr
# 2. 下面的语句续跑
gmx mdrun -s md_results_1.tpr -cpi md_results_1.cpt -deffnm md_results_1 -append