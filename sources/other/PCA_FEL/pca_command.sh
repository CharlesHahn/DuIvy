########### note
# To run this script, you need several another scripts:
#    xtc2dcd.sh : convert xtc to dcd by mdconvert 
#    xpm2png.py
#    pepend_identify.py
#    pc_combine.py
#    PCA_analysis.R : to perform PCA 
###########

echo " ====== FEL analysis ====== "
gmx trjconv -s ../second_100ns/md.tpr -f ../second_100ns/fit.xtc -o pro_second.xtc -b 70000 
gmx trjconv -s ../third_100ns/md.tpr -f ../third_100ns/fit.xtc -o pro_third.xtc -b 70000
gmx trjcat -f pro_second.xtc pro_third.xtc -settime -o pro_combine.xtc 
gmx convert-tpr -s ../second_100ns/md.tpr -o pro.tpr 

gmx covar -s pro.tpr -f pro_combine.xtc -o eigenvalues.xvg -v eigenvectors.trr -xpma covapic.xpm 
# chose C-alpha for analysis
gmx anaeig -s pro.tpr -f pro_combine.xtc -v eigenvectors.trr -first 1 -last 1 -proj pc1.xvg  
gmx anaeig -s pro.tpr -f pro_combine.xtc -v eigenvectors.trr -first 2 -last 2 -proj pc2.xvg  
gmx anaeig -s pro.tpr -f pro_combine.xtc -v eigenvectors.trr -first 3 -last 3 -proj pc3.xvg  
python pc_combine.py pc1.xvg pc2.xvg pc12_sham.xvg
python pc_combine.py pc1.xvg pc3.xvg pc13_sham.xvg

gmx sham -tsham 310 -nlevels 100 -f pc12_sham.xvg -ls pc12_gibbs.xpm -g pc_12.log -lsh pc12_enthalpy.xpm -lss pc12_entropy.xpm
mv bindex.ndx pc12_bindex.ndx
mv prob.xpm pc12_prob.xpm
mv ener.xvg pc12_ener.xvg
python xpm2png.py -f pc12_gibbs.xpm -ip yes -o pc12_gibbs.png
gmx sham -tsham 310 -nlevels 100 -f pc13_sham.xvg -ls pc13_gibbs.xpm -g pc_13.log -lsh pc13_enthalpy.xpm -lss pc13_entropy.xpm
mv bindex.ndx pc13_bindex.ndx
mv prob.xpm pc13_prob.xpm
mv ener.xvg pc13_ener.xvg

echo " ====== PCA analysis ====== "
sleep 2s
gmx trjconv -s pro.tpr -f pro_combine.xtc -o pro_pca.pdb -e 0  
python pepend_identify.py -i pro_pca.pdb -o pro_pca_end.pdb  
echo " " 
cmd < xtc2dcd.sh
C:/Users/hhhhh/C_R/R-3.6.1/bin/Rscript.exe PCA_analysis.R  

echo " ====== DONE ! GOOD DAY ! ====== "
