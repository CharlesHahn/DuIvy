#### usage

- analysis_donor_count.py 
    - 分析蛋白配体AutoDock对接结果中成键情况
- analysis_deep.py
    - 分析蛋白配体AutoDock对接结果中成键情况，包括键长、成键位置等
- average_perAA.py
    - 利用analysis_deep.py的结果，进一步分析，求平均值等
- roc.py
    - 对vina的对接结果的打分进行统计并输出ROC曲线，判断分类器效果
- ligand_pdbFromDlg.py
    - 把配体的原子位置信息从对接结果文件dlg中分离出来
- ligand_pdbFromDlg_split.py
    - 功能同上，这个好像只在正则匹配的一些小地方有差异……
- ligand_pdbqt_to_config.py
    - 这是vina自动化对接的一部分，用于从pdbqt生成vina config文件
- generate_ligand_list.py
    - vina自动化对接的一部分，生成ligand_list，然后就可以利用raccoon自动化转换
- show_me_energy.sh
    - gromacs 蛋白配体自动计算相互之间势能的一个小脚本
- gmx_shell.log.sh
    - gromacs 蛋白配体模拟的一个命令记录文件

这些代码大部分我都记不清怎么用的了，需要的时候可能需要看代码了，代码虽然写得烂，好在代码不长，可以看懂；

很多代码可能针对性太强，需要进一步修改才能适应相应的体系。