# 短肽三维结构生成

最近需要做这样的一件事情：已知一个短肽（小于八个氨基酸）的序列，如"AGKDAGKD"，需要生成它的三维坐标文件。

这个当然有很多办法可以做的。例如rdkit：可以把序列转成SMILES，还能转成3D结构，还能利用UFF力场做优化，还有个我很期待使用的功能，可以直接生成DL混合构型的短肽。当然，我还没测试过，这也不是此文的重点。

在探索这事儿的时候，意外发现了一篇13年的文章，介绍了一个从序列构建短肽三维结构的工具，**PeptideBuilder** （https://github.com/clauswilke/PeptideBuilder）。或许这刚好是我需要的，便简单看了看代码，发现挺简单的，加起来可能2.5k行代码，而且写得清楚明白。

简单介绍下它的功能：它就是一个python库，可以直接pip安装，然后用的时候直接喂给它序列，它就给你生成对应的三维结构，不能再简单了。

看一下它给的例子吧：

```python
from PeptideBuilder import Geometry
import PeptideBuilder

# create a peptide consisting of 6 glycines
geo = Geometry.geometry("G")
geo.phi = -60
geo.psi_im1 = -40
structure = PeptideBuilder.initialize_res(geo)
for i in range(5):
    PeptideBuilder.add_residue(structure, geo)
# add terminal oxygen (OXT) to the final glycine
PeptideBuilder.add_terminal_OXT(structure)

import Bio.PDB

out = Bio.PDB.PDBIO()
out.set_structure(structure)
out.save("example.pdb")
```

因为它代码的问题，所以需要先用你给定序列的头一个氨基酸初始化，也就是上面的5-8行。初始化之后就可以将后面的氨基酸挨个儿加上去，最后用BioPython包导出成pdb文件。

它这里对于每一个氨基酸，都指定了phi和psi两个角度。-60和-40这俩角度呢，对应的是alpha-helix的二级结构。如果你改成-140和130的话，对应的就是beta-sheet的二级结构。当然你也可以不指定这俩角度，那最后得到的其实就是一串直直的氨基酸。

对于一个指定的短肽，用类似下面的代码就可以生成了：

```python
geo = Geometry.geometry("D")
struture = PeptideBuilder.initialize_res(geo)
for aa in "CAHWLGE":
    geo = Geometry.geometry(aa)
    PeptideBuilder.add_residue(struture, geo)
PeptideBuilder.add_terminal_OXT(struture)

out = Bio.PDB.PDBIO()
out.set_structure(struture)
out.save("test.pdb")
```

这个库对20种天然氨基酸的geometry都做过相应的优化，得到的结果一般还是比较理想的。

不过，它只能构建L型氨基酸组成的多肽，官方不支持D型氨基酸等非天然氨基酸。

美则美矣，了则未了。

我想给它加上D型氨基酸的支持（甚至其它非天然氨基酸的支持），但是改代码、确定拓扑参数都不是短时间能搞定的，之后慢慢改吧。

当然，如果你现在就马上需要D型氨基酸短肽的话，对源代码做一个小小的修改，就可以立马用于生成D氨基酸短肽的三维结构生成。

我们知道，氨基酸的手性可以用CORN法则判定。对于L氨基酸而言，当氢原子朝向我们的时候，顺时针方向依次是CO（羰基），R侧链，N原子。而D型氨基酸则刚好是逆时针方向。故而实际上我们只需要更改氨基酸拓扑一个二面角的数据就可以了。这个库涉及到的20种氨基酸的拓扑都写在Geometry.py中，把这个文件中每一个氨基酸的`N_C_CA_CB_diangle`的角度改成相反数就可以了。当然，甘氨酸没有手性。

这样虽能得到D型氨基酸组成的短肽结构，但是因为我们只改了二面角数据，并没有对构象进行优化，在生成的结构中，某些不应该成键的原子会出现距离过近的情况。这在pymol或者vmd等分子可视化软件中，会导致错误的成键显示（生成的pdb中并没有connect信息），不过手性是没问题的，结构也是OK的。不满意？可以用ff19SB对短肽做个EM，就很合理了。

之后我再看看它的文章，然后争取尽早把我想要的功能加进去吧~

诸公也加油！