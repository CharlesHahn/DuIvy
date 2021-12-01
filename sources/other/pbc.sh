# 本文档记录了该系统的周期校正方法
# 首先利用python读取蛋白质的每一个原子的原子坐标
# 计算得到所有原子的几何中心
# 然后变例所有原子找出与几何中心最近的原子，将之定义为center
# 本系统中的center为 2115

gmx make_ndx -f md.tpr -o prolig.ndx 
echo -e "\n[ center ]\n2115\n" >> prolig.ndx
gmx trjconv -s md.tpr -f md.xtc -o center.xtc -n prolig.ndx -pbc atom -center  
gmx trjconv -s md.tpr -f center.xtc -o mol.xtc -pbc mol -ur compact -n prolig.ndx   
gmx trjconv -s md.tpr -f mol.xtc -fit rot+trans -o fit.xtc -n prolig.ndx  
gmx trjconv -s md.tpr -f fit.xtc -o fit.pdb -dt 1000 -n prolig.ndx 
