# LightDock 介绍与高通量对接尝试

### Intro

lightdock是一个可以对接蛋白-蛋白、蛋白-DNA(RNA)、蛋白-多肽以及膜上蛋白-受体等二元体系的对接框架。

lightdock同时有web server版本，也可以自行通过pip或者conda安装。

论文里面对于lightdock的介绍是**a new multi-scale approach to protein–protein docking**。 加上framework这个词，意思就很明白了：“我们写了一些有意思的对接算法和功能，但是处理对接的流程和数据分析等工作我们不太感兴趣“，所以实际上就是给了你一堆模块，自己想怎么用就自己用去。

本文的主旨是想借lightdock对接粗略阐述高通量对接自动化的一种思路。lightdock对接的详细教程和信息，还请参照官方tutorial。



### 简单教程

LightDock官方有好几个教程，基本可以涵盖使用了。这里以一个简单的蛋白-蛋白对接教程举例，稍微介绍一下对接的基本流程。

#### 对接准备

lightdock3_setup.py这个模块是用于处理要对接的蛋白质的，包括蛋白质的去氢、去水、生成ANM弹性网络、残基限制等。

比如说可以执行：

```bash
lightdock3_setup.py 2UUY_rec.pdb 2UUY_lig.pdb --noxt --noh --now -anm
```

在生成的信息中，有一个比较重要的参数：

```bash
[lightdock3_setup] INFO: Number of calculated swarms is 182
```

这表明按照用户提供的设定，程序在受体蛋白表面找到了182个swarms。可以简单粗暴地把swarms理解成受体蛋白质表面的可能结合位点。在之后的对接中，每个位点都会进行一次独立的对接。程序会生成对接参数文件setup.json以及swarms数量的文件夹。

其中还有一个参数，glowworms，默认是200，代表对一个位点进行多少次位姿尝试，可以简单理解为蒙特卡洛类似的随机初始位置。

#### 运行对接

运行对接需要用到lightdock3.py模块：

```bash
lightdock3.py setup.json 100 -c 1 -l 0
```

这里的`-c 1`指的是一个运算核心的意思。`-l 0`就很有意思了，指的是第0个swarm，也即对第一个swarm运行对接；100指的是步数，也即针对这一个swarm的每一个glowworm都执行100步。

运行的结果文件会在对应的swarm文件夹中。

自然的，要完成一个蛋白质-蛋白质对接，我们需要对每一个swarm都执行类似的上面的指令，或者`-l 0 1 2 3 4`这样一次指定多个。貌似可以全部一下子指定完？那得看内存够不够了，所有的swarm会同时运行，内存直接拉满，然后out of memory崩溃。

#### 结果分析

每个swarm的结果就存储在一堆gso_x.out里面，x是步数，并且程序只会保存10的倍数的步数的数据。每个out里面有glowworm行数据，每一行代表在这一步时该glowworm的打分等数据。

所以假如说要拿到最佳打分的数据，可能就得自己去筛选了；不仅要不同步数、不同glowworm，还得不同swarm之间的打分进行比较。其实这个功能官方也有提供，在lgd_rank_swarm.py中。

如果需要拿到最佳打分的构象，得去到最佳打分的swarm文件夹里，对包含最佳打分数据的out文件夹执行类似如下的操作：

```bash
lgd_generate_conformations.py ../2UUY_rec.pdb ../2UUY_lig.pdb gso_100.out 200
```

然后从生成的glowworm个pdb中选出最佳打分的那个pdb（序数就是最佳打分的。

#### 其它模块

上述的最简单的对接流程中，涉及到的模块只是lightdock的模块中的一部分。还有很多有意思的功能，构象聚类、计算构象的打分、各种各样的打分函数等等，值得一看。

### 自动化

或许大家也会有需要高通量对接的时候，有时候有各种已有的高通量的程序，有时候则可能需要自己去做高通量的自动化。lightdock不涉及GUI，也不太涉及蛋白口袋等依赖蛋白质本身信息的数据，以及其较为简单的操作流程，可能比较适合作为自动化高通量的一个示例。

#### 功能需求

给定一些受体蛋白文件和配体蛋白文件，对于每一个受体蛋白，完成它与所有配体蛋白的对接，并记录每组受体配体对接的最佳打分和最佳对接构象。

#### 功能分解

基于上述的对接教程，很容易构思出单个受体和配体的对接流程，直接按照对接命令的顺序执行命令就完了。当然，需要对这一个受体配体对接的每一个swarm都执行对接。当这一个受体配体的对接结束之后，需要找到最佳的打分和最佳的构象。而针对每一个受体，我们需要使之与所有的配体都进行对接。

当然，真实的对接自动化过程，还需要考虑蛋白质文件格式的转换、去除多余原子多余帧、修复蛋白结构等问题，这通常需要较多的测试。比如在本例中，除了上述考虑之外，还需要注意受体和配体蛋白质的chain id不能有重复的，不然使用lightdock提取对接构象的时候，得到的蛋白质的链中会混杂受体配体的氨基酸。

我们粗略地将功能划分成三大块：

- 对于给定的配体和受体，执行对接
- 对于某一个配体受体对接，寻找最佳打分和最佳构象
- 负责遍历受体和配体，对每一个配体和受体进行对接等流程控制

#### 逐步实现

**执行对接**

直接用python执行cmd命令，假设对接的蛋白质受体和配体的名字分别是pname和lname。

先是对接前准备的部分：

```python
setup_cmd = f"lightdock3_setup.py {pname}.pdb {lname}.pdb"
setup_cmd += f" --noxt --noh --now -anm"
print(setup_cmd)
status, output = subprocess.getstatusoutput(setup_cmd)
```

然后是对接的部分：

```python
swarm_folders = [f for f in os.listdir() if f.startswith("swarm_")]
swarm_number = len(swarm_folders)
swarm_list = [str(i) for i in range(swarm_number)]
for index in range(0, swarm_number, batch_num):
    dock_cmd = f"lightdock3.py setup.json {dock_steps} -l " 
    dock_cmd += " ".join(swarm_list[index:index+batch_num])
    print(dock_cmd)
    status, output = subprocess.getstatusoutput(dock_cmd)
```

对接这部分的代码，循环对每一个swarm都执行了对接。其中的batch_num是为了能一次执行多个swarm的对接。

**最佳打分和最佳构象**

```python
swarm_folders = [f for f in os.listdir() if f.startswith("swarm_")]

data_list = []
for swarm in swarm_folders:
    out_files = [f for f in os.listdir(swarm) if f.startswith("gso_")]
    for out in out_files:
        with open(os.path.join(swarm, out), 'r') as fo:
            lines = [ line for line in fo.readlines() if not line.startswith("#")]
            for index, line in enumerate(lines): 
                data_list.append((swarm, out, index, float(line.split()[-1])))
data = sorted(data_list, key=lambda x:x[3], reverse=True)
# print(data[:best2rank])
with open("best_score.txt", 'w') as fo:
    fo.write("#swarm, outfile, index, score\n")
    for lis in data[:best2rank]:
        fo.write(f"{lis[0]}, {lis[1]}, {lis[2]}, {lis[3]}\n")

for rank, (swarm, out, index, _) in enumerate(data[:best2rank]):
    with open("gen_best_pose.bash", 'a') as fo:
        fo.write(f"cd {swarm}\n")
        fo.write(f"lgd_generate_conformations.py ../{pname}.pdb ")
        if dock_glowworms != None:
            fo.write(f"../{lname}.pdb {out} {dock_glowworms}\n")
        else:
            fo.write(f"../{lname}.pdb {out} 200\n")
        fo.write(f"cp lightdock_{index}.pdb ../best_rank_{rank}.pdb\n")
        fo.write("rm *.pdb\n")
        fo.write("cd ../ \n")

status, output = subprocess.getstatusoutput("bash gen_best_pose.bash")
```

上述代码中生成最佳构象的部分使用了bash脚本，也算是一个偷懒的方式哈哈哈哈。如果对全部的流程分解得当，全部的自动化都可以用一堆bash脚本来控制。

**流程控制**

假设蛋白受体和配体分别放置在文件夹protein和ligand里面，那我们首先需要针对不同的对接组合建立工作目录，然后把相应的受体配体文件放进去，之后调用对接执行和获取最佳构象和打分的代码就行了。

```python
protein_files = [f for f in os.listdir(protein_source)]
ligand_files = [f for f in os.listdir(ligand_source)]
protein_names = [f[:-4] for f in protein_files]
ligand_names = [f[:-4] for f in ligand_files]
main_path = os.getcwd()
# 遍历蛋白文件
for p, pname in enumerate(protein_names):
    if not os.path.exists(f"../dock/{pname}"):
        os.mkdir(f"../dock/{pname}")
    # 遍历配体文件
    for l, lname in enumerate(ligand_names):
        if not os.path.exists(f"../dock/{pname}/{lname}"):
            ## 创建某一个对接的工作目录
            os.mkdir(f"../dock/{pname}/{lname}")
        ## 把相应的受体配体文件copy进来就可以开始对接啦
        shutil.copyfile(f"{protein_source}{pname}.pdb", 
            f"../dock/{pname}/{lname}/{pname}.pdb")
        shutil.copyfile(f"{ligand_source}{lname}.pdb", 
            f"../dock/{pname}/{lname}/{lname}.pdb")
        os.chdir(f"../dock/{pname}/{lname}")
```



**最后**把这些功能组合起来就行了：

```python
# author: charlie
# date: 20230303

import os
import shutil
import subprocess

dock_swarms = 100     # set to None for default value of lightdock
dock_glowworms = None  # set to None for default value of lightdock
dock_steps = 100
batch_num = 100
ANM_flag = True

core_num = None
best2rank = 10
best2extract = 10

protein_source = "../protein/"
ligand_source = "../ligand/"

def fetch_score_pose(pname, lname, main_log):
    swarm_folders = [f for f in os.listdir() if f.startswith("swarm_")]
    # print(f"{len(swarm_folders)} swarm folders in current directory")

    data_list = []
    for swarm in swarm_folders:
        out_files = [f for f in os.listdir(swarm) if f.startswith("gso_")]
        for out in out_files:
            with open(os.path.join(swarm, out), 'r') as fo:
                lines = [ line for line in fo.readlines() if not line.startswith("#")]
                for index, line in enumerate(lines): 
                    data_list.append((swarm, out, index, float(line.split()[-1])))
    data = sorted(data_list, key=lambda x:x[3], reverse=True)
    # print(data[:best2rank])
    with open("best_score.txt", 'w') as fo:
        fo.write("#swarm, outfile, index, score\n")
        for lis in data[:best2rank]:
            fo.write(f"{lis[0]}, {lis[1]}, {lis[2]}, {lis[3]}\n")

    for rank, (swarm, out, index, _) in enumerate(data[:best2rank]):
        with open("gen_best_pose.bash", 'a') as fo:
            fo.write(f"cd {swarm}\n")
            fo.write(f"lgd_generate_conformations.py ../{pname}.pdb ")
            if dock_glowworms != None:
                fo.write(f"../{lname}.pdb {out} {dock_glowworms}\n")
            else:
                fo.write(f"../{lname}.pdb {out} 200\n")
            fo.write(f"cp lightdock_{index}.pdb ../best_rank_{rank}.pdb\n")
            fo.write("rm *.pdb\n")
            fo.write("cd ../ \n")

    status, output = subprocess.getstatusoutput("bash gen_best_pose.bash")
    with open("fetch_pose.log", 'a') as fo:
        fo.write(output)
    if status != 0:
        with open(main_log, 'a') as fo:
            fo.write(f">>> fetch_pose {pname} {lname}\n")

# 这个函数用于呈现与某一受体对接的最佳10个配体的相关情况
def collect_show(pname, ligand_names):
    best_scores = []
    for lname in ligand_names:
        if not os.path.exists(f"../dock/{pname}/{lname}"):
            continue
        if "best_score.txt" not in os.listdir(f"../dock/{pname}/{lname}"):
            print(f"ERROR>> no best_score.txt found in ../dock/{pname}/{lname}")
            continue
        with open(f"../dock/{pname}/{lname}/best_score.txt", 'r') as fo:
            lines = fo.readlines()
        if len(lines) < 2:
            print(f"ERROR>> no data in best_score.txt of ../dock/{pname}/{lname}")
            continue
        line = lines[1]
        score = float(line.split(",")[-1])
        best_scores.append((pname, lname, score))
    best_scores = sorted(best_scores, key=lambda x:x[2], reverse=True)

    if not os.path.exists(f"../dock/final2show"):
        os.mkdir(f"../dock/final2show")
    if not os.path.exists(f"../dock/final2show/{pname}"):
        os.mkdir(f"../dock/final2show/{pname}")
    for pname, lname, score in best_scores[:best2extract]:
        shutil.copyfile(f"../dock/{pname}/{lname}/best_rank_0.pdb", 
                        f"../dock/final2show/{pname}/{lname}.pdb")
    with open(f"../dock/final2show/{pname}/score_rank.txt", 'w') as fo:
        fo.write("#Acceptor name, Ligand name, Score\n")
        for pname, lname, score in best_scores:
            fo.write(f"{pname}, {lname}, {score}\n")
    print(f"Finished collectting scores and pdbs for {pname}")

def lightdock_run():
    protein_files = [f for f in os.listdir(protein_source)]
    ligand_files = [f for f in os.listdir(ligand_source)]
    protein_names = [f[:-4] for f in protein_files]
    ligand_names = [f[:-4] for f in ligand_files]
    main_path = os.getcwd()
    main_log = f"{main_path}/../dock/ERROR_log.log"
    for p, pname in enumerate(protein_names):
        if not os.path.exists(f"../dock/{pname}"):
            os.mkdir(f"../dock/{pname}")
        for l, lname in enumerate(ligand_names):
            if not os.path.exists(f"../dock/{pname}/{lname}"):
                os.mkdir(f"../dock/{pname}/{lname}")
            shutil.copyfile(f"{protein_source}{pname}.pdb", 
                f"../dock/{pname}/{lname}/{pname}.pdb")
            shutil.copyfile(f"{ligand_source}{lname}.pdb", 
                f"../dock/{pname}/{lname}/{lname}.pdb")
            os.chdir(f"../dock/{pname}/{lname}")
            
            # setup a certain dock
            setup_cmd = f"lightdock3_setup.py {pname}.pdb {lname}.pdb"
            setup_cmd += f" --noxt --noh --now "
            if ANM_flag == True:
                setup_cmd += "-anm "
            if dock_swarms != None:
                setup_cmd += f"-s {dock_swarms} "
            if dock_glowworms != None:
                setup_cmd += f"-g {dock_glowworms}"
            print(setup_cmd)
            status, output = subprocess.getstatusoutput(setup_cmd)
            with open("setup.log", 'w') as fo:
                fo.write(output)
            if status != 0:
                with open(main_log, 'a') as fo:
                    fo.write(f">>> setup {pname} {lname}\n")
                continue
            
            # run a certain dock
            swarm_folders = [f for f in os.listdir() if f.startswith("swarm_")]
            swarm_number = len(swarm_folders)
            swarm_list = [str(i) for i in range(swarm_number)]
            for index in range(0, swarm_number, batch_num):
                dock_cmd = f"lightdock3.py setup.json {dock_steps} -l " 
                dock_cmd += " ".join(swarm_list[index:index+batch_num])
                if core_num != None:
                    dock_cmd += f" -c {core_num}"
                print(dock_cmd)
                status, output = subprocess.getstatusoutput(dock_cmd)
                with open("dock.log", 'a') as fo:
                    fo.write(output)
                if status != 0:
                    with open(main_log, 'a') as fo:
                        fo.write(f">>> dock {pname} {lname}\n")
                    continue
            
            # compare score and fetch best pose
            fetch_score_pose(pname, lname, main_log)
            os.chdir(main_path)

        # collect pdb and score of best dock        
        collect_show(pname, ligand_names)

if __name__ == "__main__":
    lightdock_run()
```

相较于之前的代码，这里增加了一点点简单的路径管理以及一些日志管理。这样的一份代码虽然能够实现功能，但是可能还不够健壮。后续还需要增加**异常处理**、**并行运算**等等内容。



### Others

同样的功能，chatgpt生成的代码也还算不错，基本上大框架都妥妥当当。

以后的同学，学习编程会更难了哈哈哈。















