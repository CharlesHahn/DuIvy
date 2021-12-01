## PLIP

PLIP (Protein-Ligand Interaction Profiler) 是一个蛋白配体非共价相互作用的分析工具；可以分析蛋白配体复合物在原子水平的非共价相互作用，包括氢键，水桥，盐桥，卤键，疏水相互作用，π-堆叠，π-离子相互作用和金属复合物；检测机制主要是基于原子间的空间位置和几何关系。

PLIP使用简单方便，且开源，提供web应用和源代码。利用其web应用可以简单地检测分析蛋白配体之间的非共价相互作用，而源代码则可用于构建高通量的分析检测体系。

### web版

PLIP的[web应用](https://projects.biotec.tu-dresden.de/plip-web/plip/index "web应用") 主要提供了两种输入，一种是PDB中已有的蛋白ID，一种是蛋白配体复合物文件。将自己的蛋白配体复合物文件提交之后，PLIP的web应用就会对复合物自动进行分析和可视化。

![PLIP results 1](http://christinagong.pythonanywhere.com/webapp/static/PLIP_1.png)

上图是web应用分析结果的页面，主要有三个区域(蓝字标示)。如果提交的pdb文件有些许错误之处，最上面提供了修复的复合物文件(其分析的是修复的复合物文件)。最下面提供了分析结果的数据文件的下载，支持XML和RST两种格式。中间则是分析结果的可视化区域了，点开之后的结果举例如下：

![PLIP results 2](http://christinagong.pythonanywhere.com/webapp/static/PLIP_2.png)

最上面就是非共价相互作用的可视化了，PLIP用不同的颜色标示蛋白和配体以及各种不同的非共价相互作用。在可视化区域下面可以下载png格式的图片和PyMOL会话文件pse，pse格式的文件可以用PyMOL读取并按照自己的需要定义和绘制图形。最下面就是具体的非共价相互作用的信息，例子中只包含了疏水相互作用和氢键，可以看到，形成了相互作用的具体原子，键长键角等信息都是有的。

更多的PLIP使用方法，以及其原理，可以参考[tutorial](https://projects.biotec.tu-dresden.de/plip-web/plip/help "tutorial") 和其发表的[论文](https://academic.oup.com/nar/article/43/W1/W443/2467865 "论文")。

### PLIP源代码

[PLIP的Github主页](https://github.com/pharmai/plip "PLIP的Github主页") 包含了PLIP的Python源代码和安装方法。

如果只是少量的蛋白配体复合物需要进行非共价相互作用分析，那么web版是足够的；但如果分析的量比较大，几十几百几千，本地搭建PLIP的分析环境就很有必要了。

当然，如果对爬虫特别擅长且不想在自己电脑上折腾PLIP，也是可以通过爬虫的思路，利用PLIP的web服务搭建高通量分析体系的；但这样，非常不厚道哈哈哈。

在本地搭建PLIP分析体系不难，大致有三种方式。

如果安装有Docker的话，可以直接下载Docker镜像，非常的方便，不用自己处理依赖，缺点是镜像有点儿大。

还可以通过pip安装或者直接源码安装。目前由于PLIP更换了maintainer，所以最新的 v2.x.x 有较大的变动，not fully tested，pypi上也没有同步最新的PLIP，只有v1.4.2。

对于最新的PLIP v2.1.0-beta 来说，只需要本地电脑上安装**python 3.x** 、**OpenBabel 3** 、然后安装python的第三方库**openbabel**  就行了。

PLIP需要借助 Openbabel 3 软件进行一些格式转换，安装这个软件之后，记得添加到环境变量，能在cmd里面敲`obabel` 有相应的输出就成了。

```shell
$ obabel
No input file or format spec or possibly a misplaced option.
Most options must come after the input files. (-i -o -O -m can be anywhwere.)

Open Babel 3.0.0 -- Oct  7 2019 -- 20:18:16
Usage:
obabel [-i<input-type>] <infilename> [-o<output-type>] -O<outfilename> [Options]
Try  -H option for more information.
```

python的第三方库openbabel则是用来调用OpenBabel 3软件的，记得保证第三方库和软件的位数一样，还要和python3解释器的位数一样，最好都是64位的吧；不成的话就解释器什么位数，软件和第三方库就什么位数。一般来说，pip自动安装的第三方库的位数和python解释器的位数都是一样的，注意OpenBabel 3软件的位数就成了。

如果需要手动安装openbabel第三方库，可以到[pythonlibs](https://www.lfd.uci.edu/~gohlke/pythonlibs/ "pythonlibs")下载二进制文件手动安装。

在这些依赖都搞定之后，就可以开始从源码安装最新的PLIP了。首先Github的release界面下载最新的source code，然后解压，在解压后的目录里跑个`python setup.py install`就成了。这时候运行`pip list`应该就可以看到PLIP被安装好了。

测试的话，可以到解压后的目录中的 plip/test 目录。在test文件夹里有很多测试的代码，有个run_all_tests.sh的脚本，Linux的话，直接运行就成了，Windows的话可以复制run_all_tests.sh里的命令到命令行里运行。按照其论文中的表述，这些测试会对一些已知非共价相互作用的蛋白配体复合物进行测试，所得的结果与文献报道应该一致，不一致的话，测试不会通过。我之前针对 v1.4.5 进行过测试，发现分析结果和文献报道不符，也和web版不相符，但是web版也和文献报道不太相符，可能版本迭代过程中有些小问题；但好在，不符的地方很少；最新版本的测试代码经过新maintainer的重构，笔者还未测试过。

如果需要在命令行里使用的话，可能还需要 alias 一下：

```shell
alias plip='python ~/plip/plipcmd.py'
```

Linux用户自然明白怎么搞；Windows用户如果用cmder的话，也可以把这句命令加在cmder目录下的config/user_aliases.cmd文件中，如果没有cmder，也可以找找其他办法，或者直接每次运行的时候，`python C:\xxx\xxx\plipcmd.py` 也很成，虽然麻烦一点儿。

以上大致就是PLIP最新版的安装方法了。

对于某些已经安装了python 2.x或者 OpenBalbel 2.x的用户来说，可以安装PLIP较早的版本。PLIP最早是用 python 2写的，但是之后的 v1.4.x应该是兼容python 3的，但可能不兼容 OpenBabel 3(OpenBabel 3 改了python的接口)。具体的版本需要哪些依赖可以查看release信息。源码的安装方法还是同之前类似。

因为之前搭建分析环境的时候，PLIP不同时兼容python 3和OpenBabel 3，故此笔者重写了部分非核心代码，并提交了PR，有幸原maintainer merge了我的代码，因而那一个commit([761b9c6](https://github.com/CharlesHahn/plip/commit/761b9c6ebfa93ff6ca40ff01552bbb49b6ddb83b "761b9c6")) 是v1.4.4中同时兼容python 3和OpenBabel 3的。故而有需要的朋友也可以到[笔者的PLIP项目](https://github.com/CharlesHahn/plip/tree/master "笔者的PLIP项目")下载安装源代码 (注意同时兼容的是 master 分支下的代码)。

```shell
$ plip
usage: PLIP [-h] (-f INPUT [INPUT ...] | -i PDBID [PDBID ...])
            [-o OUTPATH | -O] [--rawstring] [-v] [-p] [-x] [-t] [-y]
            [--maxthreads MAXTHREADS] [--breakcomposite] [--altlocation]
            [--debug] [--nofix] [--nofixfile] [--nopdbcanmap] [--dnareceptor]
            [--name OUTPUTFILENAME] [--peptides PEPTIDES [PEPTIDES ...] |
            --intra INTRA] [--keepmod]
```



PLIP分析环境搭建好了之后，就可以通过python调用PLIP来进行高通量的分析啦！

最简单的调用方式莫过于直接用python执行命令行命令，利用python自带的subprocess模块执行getstatusoutput 函数即可。

```python
status, output = subprocess.getstatusoutput("plip -f complex.pdb -yvx")
```

上述代码可以生成分析结果的xml数据文件和pse文件。xml数据可以通过模块 xml.etree.ElementTree 等工具进行解析和处理。



### Others

若要对分子对接的结果进行高通量处理和分析：首先不论是AutoDock4还是vina等等，都可以通过python脚本解析拼接处理得到多个复合物pdb文件，然后利用python调用PLIP自动化对每一个复合物文件进行分析，再用python处理xml文件提取非共价相互作用的数据并处理分析，最后就可以较为方便地得到相互作用数量随氨基酸序列或者配体原子的变化、平均键长键角等信息。

若在搭建其它体系的时候有需要，PLIP这个工具也可以很方便地嵌入其中并协调工作。

另，如果体系中有多个配体，似乎无法导出pse文件，笔者看报错似乎与pymol中的center共功能有关，因为多个配体存在，故而无法决定将哪个配体放中间。

天雨，继续加油~

