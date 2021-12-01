# shell scripts for energy extract

echo " ================ extract the prolig energy =================== "
gmx make_ndx -f npt.gro -o prolig.ndx 
gmx trjconv -f md_results_1.xtc -s md_results_1.tpr -n prolig.ndx -o prolig.xtc
gmx convert-tpr -s md_results_1.tpr -n prolig.ndx -o prolig.tpr
gmx mdrun -s prolig.tpr -rerun prolig.xtc -e prolig.edr 
gmx energy -f prolig.edr -o prolig_SR.xvg

echo " ================ extract the pro energy =================== "
gmx make_ndx -f npt.gro -o pro.ndx
gmx trjconv -f md_results_1.xtc -s md_results_1.tpr -n pro.ndx -o pro.xtc
gmx convert-tpr -s md_results_1.tpr -n pro.ndx -o pro.tpr
gmx mdrun -s pro.tpr -rerun pro.xtc -e pro.edr 
gmx energy -f pro.edr -o pro_SR.xvg

echo " ================ extract the pro energy =================== "
gmx make_ndx -f npt.gro -o lig.ndx
gmx trjconv -f md_results_1.xtc -s md_results_1.tpr -n lig.ndx -o lig.xtc
gmx convert-tpr -s md_results_1.tpr -n lig.ndx -o lig.tpr
gmx mdrun -s lig.tpr -rerun lig.xtc -e lig.edr 
gmx energy -f lig.edr -o lig_SR.xvg

