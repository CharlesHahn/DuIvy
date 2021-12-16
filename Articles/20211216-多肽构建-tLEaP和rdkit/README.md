## 多肽：从序列到三维结构

### tLEaP ： 多肽的构建

前儿分享了多肽构建工具PeptideBuilder，今天首先分享的是利用AmberTools里面的tLEaP工具构建多肽（包括天然多肽，D-多肽和DL-混合多肽）。

tLEaP里面有个`sequence`命令可以帮助我们产生想要的多肽。当然，只能产生类似于一条直线的多肽链儿，而不能像PeptideBuilder一样通过设置phi和psi两个角度来定义二级结构。不过，tLEaP可以通过flip命令生成D型多肽，这是PeptideBuilder目前没有的功能。

当然，使用tLEaP需要先安装AmberTools了。（之前liuyujie老师在他的博客里提供了一个预编译版，但是现在他的博客好像宕了。我有个备份，诸位有需要的话可以下载 http://charles8hahn.pythonanywhere.com/download/amber20_win64.7z 。下载的时候默念感谢liuyujie老师哈哈哈）

先总览一下生成多肽需要的tLEaP命令（`#`号后面是我加的注释）。

```bash
# 选用力场
source leaprc.protein.ff14sb
# 生成多肽序列KLVAFD，默认天然构型
mol = sequence { NLYS LEU VAL ALA PHE CASP }
# 检查一下
check mol
# 保存L构型的多肽
savepdb mol L-KLVAFD.pdb

# 利用flip构建全D型的多肽
# 首先选中所有氨基酸的alpha-C
select mol.1.CA
select mol.2.CA
select mol.3.CA
select mol.4.CA
select mol.5.CA
select mol.6.CA
# flip一下
flip mol
# 保存全D型的多肽
savepdb mol D-klvafd.pdb
```

下面简单解释下各条命令：

```bash
source leaprc.protein.ff14sb
```

载入力场，这里用的是蛋白的ff14sb力场。

```bash
mol = sequence { NLYS LEU VAL ALA PHE CASP }
```

通过`sequence`命令生成一个6肽链儿，并命名为mol。第一个氨基酸LYS前面加了N，表示这是N端的氨基酸。ASP前面加了C，表示这是C端的氨基酸。tLEaP会根据这个标识对链的两端进行合适的处理。新版的tLEaP好像已经能自动识别链的N端C端了。

```bash
savepdb mol L-KLVAFD.pdb
```

把mol，也就是天然多肽保存成pdb文件。

```bash
select mol.1.CA
```

选择mol对象的第一个氨基酸的alpha-C。

```bash
flip mol
```

对选中的原子进行翻转，也就是从L型换到D型。上面选择了所有氨基酸的alpha-C，所以这里flip之后得到的就是完全由D型氨基酸构成的多肽了。

如果你想要**DL混合型的多肽**，那可以先构建天然多肽序列，然后选择想要换成D型的氨基酸的alpha-C进行翻转即可。

SB系列的力场，对D型的氨基酸也是可以用的，不用自己再去搞D型氨基酸的参数了。如果你需要保存分析动力学模拟需要的参数和拓扑：

```bash
saveamberparm mol mol.prmtop mol.inpcrd
```



### rdkit : 多肽的构建

rdkit，非常著名的工具，当然也可以用来构建多肽（D型也可以）。

一遍上python代码一边解释吧。

先把完整的代码（生成天然六肽，D-六肽，DL-六肽）放上来：

```python
from rdkit import Chem
from rdkit.Chem import AllChem

# 生成天然六肽
L_seq = "KLVAFD"
L_mol = Chem.rdmolfiles.MolFromFASTA(L_seq, flavor=0)
# L_smi = Chem.MolToSmiles(L_mol)
# print(L_smi)
L_mol = Chem.AddHs(L_mol)
AllChem.EmbedMolecule(L_mol)
L_pdbblock = Chem.MolToPDBBlock(L_mol)
print(L_pdbblock, file=open("L-KLVAFD.pdb", "w"))

# 生成 D-六肽
D_seq = "klvafd"
D_mol = Chem.rdmolfiles.MolFromFASTA(D_seq, flavor=1)
D_mol = Chem.AddHs(D_mol)
AllChem.EmbedMolecule(D_mol)
D_pdbblock = Chem.MolToPDBBlock(D_mol)
print(D_pdbblock, file=open("D-KLVAFD.pdb", "w"))

# 生成 DL-六肽
D_seq = "KlvAfD"
D_mol = Chem.rdmolfiles.MolFromFASTA(D_seq, flavor=1)
D_mol = Chem.AddHs(D_mol)
AllChem.EmbedMolecule(D_mol)
D_pdbblock = Chem.MolToPDBBlock(D_mol)
print(D_pdbblock, file=open("DL-KlvAfD.pdb", "w"))
```

下面就一些重要的语句和参数进行解释：

```python
from rdkit import Chem
from rdkit.Chem import AllChem
```

导入需要的模块。

```python
L_seq = "KLVAFD"
L_mol = Chem.rdmolfiles.MolFromFASTA(L_seq, flavor=0)
```

这两句首先定义了多肽序列”KLVAFD"，然后利用MolFromFASTA函数将序列转换成分析对象。

**flavor**这个参数有着特殊的意义，详情请参照rdkit文档中相关的阐述（https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html#rdkit.Chem.rdmolfiles.MolFromFASTA）。仅对于蛋白或者多肽来说，这个参数可以取两个值（0 or 1），两个值的含义不同。如果flavor定义为0（默认情况），则不管定义的序列里面的字母是大写还是小写，通通认为是L型氨基酸（也即天然构型）进行读入。**如果flavor定义为1，则序列中的大小写会被区别对待，大写被认为是L型氨基酸，小写会被认为D型氨基酸。** 关于这个细节的描述，还可以参照https://github.com/rdkit/rdkit/issues/3088。

也正是有了flavor这个参数，我们才可以自由定义想要的氨基酸手性。

因而如下的这两行命令生成的多肽全是D型氨基酸：

```python
D_seq = "klvafd"
D_mol = Chem.rdmolfiles.MolFromFASTA(D_seq, flavor=1)
```

而如下两行命令生成的多肽是DL混合的，大写的是L型（K, A, D），小写的是D型（l, v, f）：

```python
D_seq = "KlvAfD"
D_mol = Chem.rdmolfiles.MolFromFASTA(D_seq, flavor=1)
```

`Chem.rdmolfiles.MolFromFASTA()`的功能就是从序列生成分子对象了，与它类似的函数还有`Chem.MolFromSequence()`，以及`Chem.MolFromFASTA()`，也都有flavor这个参数的。

```python
L_smi = Chem.MolToSmiles(L_mol)
print(L_smi)
```

这两句在上文中我注释掉了，主要作用就是把分子对象转化为SMILES表达式，在这里不起作用。（不过你有兴趣可以输出来看看，不同手性的氨基酸，其实SMILES表达式中就只有`@`符号数量的区别，参考https://zh.wikipedia.org/wiki/简化分子线性输入规范。）

```python
L_mol = Chem.AddHs(L_mol)
```

这一行是给分子加氢。（如果保存pdb之后再用其它工具加氢，也不一定能保证手性部分的氢就能加得对，索性在这里加了，比较靠谱。）

```python
AllChem.EmbedMolecule(L_mol)
```

这一行是生成分子的3D构象，或者可以认为是构象优化。先用距离几何初始化3D坐标，然后使用ETKDG算法进行优化。可以参考博文https://blog.csdn.net/dreadlesss/article/details/105613264 或者看官方文档。

```python
L_pdbblock = Chem.MolToPDBBlock(L_mol)
print(L_pdbblock, file=open("L-KLVAFD.pdb", "w"))
```

这两行的第一行就是利用`Chem.MolToPDBBlock()`把分子对象转换成PDB格式的文本，第二行就是存成一个文件。

到这里，我们就利用rdkit，从序列生成了多肽的PDB文件了。

如果你生成的是天然氨基酸构成的多肽，恭喜你，事儿成了！如果你生成的多肽包含D 型氨基酸，很抱歉，还没完！

rdkit对于不同手性的氨基酸标记了不同的名字。对于天然氨基酸，自然就是通用的命名，如Lys的名字就是LYS。但是对于D型氨基酸，rdkit给的命名是D和氨基酸名字的组合，如D型Lys的名字为DLY，D型ALA的名字为DAL。除此之外，PDB文件中本应是ATOM开头的部分行也被替换成HETATM开头。

同时，加的H，也被标记了不同的分子名字，为UNL。如果你需要H的分子名字跟随氨基酸名字，可以考虑不在rdkit里面加氢，导出成pdb之后用reduce等工具加氢。

所以，对于含D型氨基酸的pdb文件，可能还需要进一步处理。



### Others

前儿的PeptideBuilder，今儿的tLEaP和rdkit，也草草算是介绍了三种从序列构建多肽三维结构的工具了。

PeptideBuilder构建的多肽是不带氢原子的，官方也是没有支持D型氨基酸的。tLEaP和rdkit构建出来的多肽都是（可以）有氢原子的，并且支持D型氨基酸。

这两篇内容，都成文仓促。基本上都是下午开始，一遍测试一边写，晚上就发出来。难免疏漏错误，还望一定批评斧正。

祝诸公顺遂！



