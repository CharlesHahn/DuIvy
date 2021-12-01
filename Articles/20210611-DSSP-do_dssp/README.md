# DSSP计算蛋白质二级结构

## 0. preface

蛋白质的二级结构，最典型的应当就是α螺旋和β折叠了。通常一个蛋白质当中，α螺旋和β折叠以及loop都是很重要的组成部分。

![blueFlamingo_logo](C:\Users\hhhhh\Desktop\databank\公众号\20210611\blueFlamingo_logo.png)

如上图的典型的α螺旋，一圈包含3.6个氨基酸，每两个氨基酸之间的高度差通常是1.5埃。当然在实际的情况下，这些数值也不是固定的。α螺旋通过上面一圈的氨基酸与下面一圈氨基酸之间的氢键来稳定其结构，氢键的走向和轴向一致。

![WeChat Screenshot_20210612094418](C:\Users\hhhhh\Desktop\databank\公众号\20210611\WeChat Screenshot_20210612094418.png)

上图则显示了经典的antiparallel和parallel的β折叠。β折叠则是通过链与链之间的氢键来稳定的。

DSSP（Dictionary of Protein Secondary Structure）定义了一系列蛋白质的二级结构，包括3转角螺旋、α螺旋、5转角螺旋、转角、β折叠等等。在分子动力学模拟或者生物信息研究中，常常需要分析蛋白质的二级结构的变化。方法当然有很多，较为常见的还是利用DSSP软件对蛋白质pdb文件进行计算。

## DSSP的安装

### 在windows上安装DSSP

两种方法，一种是直接源码编译，一种是直接下载别人编译好的可执行文件。

windows上可以源码编译，有兴趣有能力的朋友可以探索探索！

预编译好的可执行程序的话，我这里有**32位的DSSP3.0.0**的版本，经过测试是可以搭配gmx2019的版本运行的，gmx2018或者2020我想应当也是没有问题的，64位的windows操作系统也是可以兼容的。感兴趣的朋友下方随意赞赏此文即可获得下载链接哦！

最新的DSSP4.0的windows预编译版本在https://github.com/PDB-REDO/dssp的issue里可以找到相关的下载链接，当然你也可以赞赏本文直接获得。**notice** : 下载得到的mkdssp.exe需要和`mmcif_pdbx_v50.dic`文件在同一路径下才能正常使用，且其是测试版本，可能不稳定。

在拿到DSSP的可执行程序之后，将程序添加到环境变量即可正常使用了（新建DSSP变量指向dssp.exe程序，也可以同时添加dssp.exe所在的路径到path变量），gmx应该也能正常调用。

![WeChat Image_20210612162523](C:\Users\hhhhh\Desktop\databank\公众号\20210611\WeChat Image_20210612162523.png)

如果dssp没有正确安装或者没有添加到环境变量，`do_dssp`就无法找到这个程序，会出现如下的fatal error ：

![WeChat Photo Editor_20210612145849](C:\Users\hhhhh\Desktop\databank\公众号\20210611\WeChat Photo Editor_20210612145849.jpg)

**NOTICE** ：DSSP4.0+的版本使用的命令较之前的版本有所改变。如果你使用的不是最新的GROMACS，那`do_dssp`程序可能无法正确调用DSSP4.0+，请使用较老的DSSP版本。

这是最新的DSSP4.0+版本：

![WeChat Photo Editor_20210612151534](C:\Users\hhhhh\Desktop\databank\公众号\20210611\WeChat Photo Editor_20210612151534.jpg)

这是较老的DSSP3.0.0的版本：

![WeChat Photo Editor_20210612151339](C:\Users\hhhhh\Desktop\databank\公众号\20210611\WeChat Photo Editor_20210612151339.jpg)

### 在Linux上安装DSSP

做生物信息或者分子模拟相关工作的朋友相信或多或少可能接触过Linux操作系统，在Linux上安装DSSP对于有Linux使用经验的朋友来说应当是不难的。可以通过源码编译DSSP也可以直接安装二进制程序。诸多做生信或者模拟的朋友使用的Linux可能都是centos或者红帽，又或者Ubuntu或者Debian之流。我可能比较另类，爱好ArchLinux。ArchLinux的AUR仓库中是有dssp3.1.4的。

![WeChat Screenshot_20210612110614](C:\Users\hhhhh\Desktop\databank\公众号\20210611\WeChat Screenshot_20210612110614.png)

故而可以通过下面的一行命令快速安装DSSP。

```shell
yay -S dssp
```

对于其它的发行版，可以查询下对应的源或者仓库里有没有别人预编译好的dssp可以直接安装，然后按照自己发行版对应的安装命令进行安装。

源码编译的话也并不复杂。从https://github.com/PDB-REDO/dssp.git克隆到源码之后，按照github上的介绍进行编译安装即可，这样得到的dssp是最新的（应该是4.0+的版本了，当然也可以编译3.0+的旧版本）：

```
./configure
make
sudo make install
```

在编译过程中可能会遇到一些和boost相关的问题，可以看看github的issue或者相关的帖子。这里提供一个比较详细的DSSP在Linux上的安装教程：https://blog.csdn.net/weixin_39417324/article/details/115136852，感谢这位<跑着去干饭>大佬哈哈哈。

不管是通过包管理器安装的还是源码安装的，最后得到的可能都是一个叫mkdssp的程序，重命名为dssp或者链接到dssp都行，之后添加到环境变量即可，这样就可以正常使用了，gmx的`do_dssp`命令也就可以识别了。

### 通过conda安装

通过conda可以一行命令安装。

dssp-3.0.0的安装命令：

```
conda install -c salilab dssp
```

dssp-2.2.1的安装命令：

```
conda install -c ostrokach dssp
```

上述两条命令在Manjaro上测试通过，在Linux上应该都还蛮可行的。但很可惜在win10上并没有找到对应的dssp程序，可能我个人配置原因，诸位还是可以一试。

## DSSP的使用

DSSP的用法应该是简单的

```
dssp -i input.pdb -o output.dssp
```

## gmx的do_dssp命令

`do_dssp`所需要的输入和相关的参数可以通过`gmx help do_dssp`命令查看。常用的输出通常是二级结构含量随时间变化的xpm文件和xvg文件了。

```
gmx do_dssp -s test.tpr -f test.xtc -o test.xpm -sc test_sc.xvg
```

test.xpm文件就是蛋白质二级结构图了，test_sc.xvg则是蛋白质二级结构百分比含量随时间变化的数据文件。

得到的xpm文件可以通过gmx的xpm2ps命令转换成ps或者eps文件再打开，诸多教程也是采用的这种方法。个人认为Jerkwin老师的xpm2all.bsh脚本更加好用，它可以把xpm文件转换成gnuplot的绘图脚本，如果安装了gnuplot软件，可以很方便地绘图。

## xpm文件的绘图

不管是Linux还是windows，gnuplot都可以很方便地安装。对于windows，可以在sourceforge上直接下载到安装包，地址如下https://sourceforge.net/projects/gnuplot/。xpm2all.bsh脚本可以在https://jerkwin.github.io/gmxtool/下载得到，这个界面上也有Jerkwin老师精心写的使用教程。

下面简单说明整个绘图的流程：

1. 将xpm文件转换成gpl脚本

```
bash xpm2all.bsh test.xpm -gpl xyz gmx
```

2. gnuplot载入test~.gpl进行绘图

```
echo load "test~.gpl" | gnuplot
```

```
$ gnuplot # 首先进入gnuplot的命令模式
gnuplot> load “test~.gpl"
gnuplot> exit
```

对于最新的xpm2all.bsh输出的gpl脚本，里面通常已经设定好了会将图像输出到图片，如果没有输出到图片，可以在gpl脚本前面加上如下两行，在进行gnuplot绘图。

```
set term png
set output "yourfilename.png"
```

xpm2all.bsh脚本的执行需要bash解释器，windows上也有诸多替代品，推荐cmder（包含bash、git、vim），简单方便。

如此我们就得到了模拟过程中蛋白质二级结构的图，但是默认出来的图可能并不太符合要求，可能还想微调一些细节等等。下面简单介绍一下这里的gpl脚本的内容。

```
set term png
set output "output.png"
unset colorbox
set pal defined(0 "#000000",1 "#00CC00",2 "#00B266",3 "#D4A800",4 "#FFF000",5 "#F28E8E",6 "#EC6262",7 "#E53535")

$data <<EOD
0 1 0 
...
... # 这里是特别多行数据
...
EOD

#set tmargin at screen 0.95
#set bmargin at screen 0.2
#set rmargin at screen 0.85
#set label " A-Helix " at screen 0.85,0.92 left textcolor rgb "#E53535"
#set label " 3-Helix " at screen 0.85,0.82 left textcolor rgb "#EC6262"
#set label " 5-Helix " at screen 0.85,0.72 left textcolor rgb "#F28E8E"
#set label " B-Sheet " at screen 0.85,0.62 left textcolor rgb "#FFF000"
#set label " B-Bridge" at screen 0.85,0.52 left textcolor rgb "#D4A800"
#set label " Turn    " at screen 0.85,0.42 left textcolor rgb "#00B266"
#set label " Bend    " at screen 0.85,0.32 left textcolor rgb "#00CC00"
#set label " Coil    " at screen 0.85,0.22 left textcolor rgb "#000000"
set term pngcairo enhanced truecolor font "HelveticaNeueLT Pro 85 Hv,85" fontscale 1 linewidth 20 pointscale 5 size 6000,3033
set tics out nomirror;
set key out reverse Left spacing 2 samplen 1/2
set xl"Time(ps)"; set yl"#res";
plot [-2.5:50002.5] [1:513] $data u 1:2:3 w imag notit, \
-1 w p ps 3 pt 5 lc rgb "#E53535" t"H A-Helix ", \
-1 w p ps 3 pt 5 lc rgb "#EC6262" t"G 3-Helix ", \
-1 w p ps 3 pt 5 lc rgb "#F28E8E" t"I 5-Helix ", \
-1 w p ps 3 pt 5 lc rgb "#FFF000" t"E B-Sheet ", \
-1 w p ps 3 pt 5 lc rgb "#D4A800" t"B B-Bridge", \
-1 w p ps 3 pt 5 lc rgb "#00B266" t"T Turn    ", \
-1 w p ps 3 pt 5 lc rgb "#00CC00" t"S Bend    ", \
-1 w p ps 3 pt 5 lc rgb "#000000" t"C Coil    "
```

以上就是这里的gpl文件的基本内容，包含了输出文件格式(png)，输出文件名字(output.png)，用于绘图的各类颜色和图例等等。‘#’开头的行都是注释。

改颜色：直接在文件开头和文件末尾改就行了，注意要对应着改。

改图片尺寸：`set term pngcairo enhanced truecolor`那一行的末尾两个数字。

改图片字体：`set term pngcairo enhanced truecolor`那一行的font后面双引号里面就是，逗号前面是字体，逗号后面是字号。

改横纵坐标：`set xl"Time(ps)"; set yl"#res";`一行，把引号里面的改改就成。

想控制绘图区域：将`set tmargin`，`set bmargin`，`set rmargin`三行取消注释，想要多大的绘图区域就自己调多大的绘图区域。

想要图例在下面：`set key out reverse`一行改成如下两行:

```
set key spacing 1 samplen 1/2
set key bottom outside horizontal center
```

想自己控制图例的位置：将`set label`的多行取消注释，自己调控制位置的那俩参数。

其它的部分我想应该都容易看懂，微调调就会比较不错啦。

![ss](D:\TAU\Tau\report\gmxtools\ss.png)

## 二级结构含量数据文件的绘图

参考Jerkwin老师https://jerkwin.github.io/2021/03/30/xpm2all更新-二级结构绘制_颜色方案/ 一文的末尾。

通常得到的二级结构含量文件如下所示：

![WeChat Image_20210612200052](C:\Users\hhhhh\Desktop\databank\公众号\20210611\WeChat Image_20210612200052.png)

在数据部分有诸多列，每一列表示不同的含义。最左边一列表示Time(ps)，紧接着的列分别表示Structure、Coil、B-Sheet、B-Bridge、Bend、Turn、3-Helix、Chain-Separator。当然不同的蛋白得到的这个文件可能不同，看看数据前面的说明就知道每一列表示什么意思了。

Jerkwin老师在他的博文最后给出了gnuplot绘制二级结构含量文件的脚本，但可能需要修改一下才能和自己的二级结构含量文件对应。对于这里的示例文件，相应的脚本如下：

```
set term png
set output "dsspsc.png"
set term pngcairo enhanced truecolor font "Arial,85" \
fontscale 1 linewidth 20 pointscale 5 size 6000,3500
set tics out nomirror;
set key out reverse Left spacing 2 samplen 1/5
set style fill solid 1.0 border
set xl"Time(ps)"; set yl"#Res%"

plot [0:40000] [-0.5:100] 'test_sc.xvg' \
   u 1:(($9+$8+$7+$6+$5+$4+$3)*100/739) s f w boxes lc rgb "#FF0000" t"Coil  ", \
'' u 1:(($9+$8+$7+$6+$5+$4   )*100/739) s f w boxes lc rgb "#2E8857" t"B-Sheet     ", \
'' u 1:(($9+$8+$7+$6+$5      )*100/739) s f w boxes lc rgb "#808080" t"B-Bridge ", \
'' u 1:(($9+$8+$7+$6         )*100/739) s f w boxes lc rgb "#0000FF" t"Bend     ", \
'' u 1:(($9+$8+$7            )*100/739) s f w boxes lc rgb "#3399FF" t"Turn     ", \
'' u 1:(($9+$8               )*100/739) s f w boxes lc rgb "#FF8C00" t"3-Helix  ", \
'' u 1:(($9                  )*100/739) s f w boxes lc rgb "#00FF00" t"Chain-Seperator  "
```

前面的内容上面已经讲过了，下面着重讲讲怎么对应更改，希望对不熟悉gnuplot的朋友有所帮助。

首先是`plot`一行，两个中括号里面的两对数字分别表示横纵坐标的起始和结束数值。之后的单引号里面是要绘制的二级结构含量文件的文件名，需要对应修改。

下面的几行中，最需要理解的是`$9`，`$2`等东西的含义。其实很简单，如果有朋友知道awk的域的概念的话。从上面的示例文件我们可以看到，在主要的数据部分，一共有9列数据。对于每一行数据而言，数字和数字之间通过空格或者tab键分开，这样的一个数字就叫做一个域，就可以通过`$`符号和一个数字来表示。比如`$9`这里就表示最后一列数据，`$1`就表示Structure列的数据。

那么在绘图的时候，为啥要把不同的域加起来呢？这里采用的是常用的堆积柱状图的方法。我们要在最上面表示`Coil`的话，就需要把所有数据加起来绘制在最底层，也就是`$9+$8+$7+$6+$5+$4+$3`，然后再在第二层绘制其它的部分，也就是下一层`$9+$8+$7+$6+$5+$4`，这样绘图区上方空出来的部分，也就是`$3`，也就是Coil这个域了。再画一层`$9+$8+$7+$6+$5`，这样绘图区上方又一次堆积出来的就是`$4`，也就是B-Sheet。一层一层的画下去，最后得到的就是我们需要的堆积柱状图了。

这里我没有绘制Structure这个域。在明白了这一堆dollar sign之后，还需要对应调整右边的二级结构的名字。比如第一行要绘制的是Coil结构域，右边就要对应写成`t"Coil"`，别搞错了。

还有一个地方需要注意，`*100/739`，739表示的是蛋白质的总氨基酸数目。这里也需要根据你的蛋白质进行修改。

![dsspsc](D:\TAU\210527\single\dsspsc.png)

可以自己调调每一个堆积的颜色，我这方面不擅长，也就这样凑合了。



## others

这两天下了好大的雨，把我的盆栽端出去淋了淋。科研不顺，总感觉自己在无中生有暗度陈仓orz，或许也需要淋一淋雨。

> 赞赏本文，送dssp-3.0.0.exe、mkdssp-4.0.exe，还送一本图文并茂非常nice的《protein structure and function》


