# PLIP相互作用变化图

### Intro

前几天看到这么一张图，纵轴是氨基酸名字，横轴是拉伸动力学的拉伸距离，表现的也就是蛋白配体的相互作用随着拉伸距离的变化。

![WeChat Image_20230409220316](C:\Users\hhhhh\Desktop\20230412-PLIP绘图\WeChat Image_20230409220316.png)

更常见的，可能是相互作用随时间的变化的图，也即横坐标是时间。当然本质其实是一样的，都是相互作用随着某一变量的变化。

### 绘制

要绘制这样的一幅图，我们至少需要这样一些信息：相互作用类型、涉及的氨基酸及其编号，以及这个相互作用对应的x轴数据（可以是距离，也可以是时间等变量）。

#### 获取数据

蛋白配体的相互作用可以通过很多方式获得，这里我们通过PLIP来获取蛋白配体的相互作用。我们沿着拉伸方向截取若干帧复合物轨迹，导出成pdb，分别为conf0.pdb，conf1.pdb等，然后用PLIP分析得到相互作用的结果。

```python
plip -txyp -f conf0.pdb
```

上述命令会得到生成的相互作用图片、pymol的pse文件，以及两个记载着相互作用详细信息的文本文件：一个是方便人类阅读的txt文件，一个是方便程序读取的xml文件。

我们先按照顺序依次从xml文件中读取相互作用。

XML是一种可拓展的标记语言，可能大家更熟悉的是HTML，XML也是类似的由一个个节点嵌套来存储层次数据的一种数据表示格式。下面是我截取的PLIP结果中的一部分：

```html
<?xml version='1.0' encoding='ASCII'?>
<report>
  <plipversion>2.2.2</plipversion>
  <bindingsite id="1" has_interactions="True">
    <identifiers>
      <longname>UNK</longname>
    </identifiers>
    <lig_properties>
      <num_heavy_atoms>20</num_heavy_atoms>
    </lig_properties>
    <interacting_chains>
      <interacting_chain id="1">A</interacting_chain>
    </interacting_chains>
    <bs_residues>
      <bs_residue id="1" contact="False" min_dist="5.2" aa="ASP">67A</bs_residue>
    </bs_residues>
    <interactions>
      <hydrophobic_interactions>
        <hydrophobic_interaction id="1">
          <resnr>316</resnr>
          <restype>THR</restype>
          <reschain>A</reschain>
        </hydrophobic_interaction>
      </hydrophobic_interactions>
      <hydrogen_bonds>
        <hydrogen_bond id="1">
          <resnr>275</resnr>
          <restype>GLU</restype>
          <reschain>A</reschain>
        </hydrogen_bond>
      </hydrogen_bonds>
      <water_bridges/>
      <salt_bridges/>
      <pi_stacks/>
      <pi_cation_interactions/>
      <halogen_bonds/>
      <metal_complexes/>
    </interactions>
  </bindingsite>
  <filetype>PDB</filetype>
</report>
```

XML这类文件的解析方式有很多，相关的弯弯绕也不是此文能详尽的。因而只介绍python自带的xml解析工具xml.etree.ElementTree，其将xml文件很自然的解析成一棵树，我们只需要按照层级访问就可以了。

相关的资料可以参考：https://docs.python.org/zh-cn/3/library/xml.etree.elementtree.html 或者https://zhuanlan.zhihu.com/p/502584681 等。

对于上面的xml内容来说，我们可以有如下的解析代码：

```python
import xml.etree.ElementTree as ET

tree = ET.parse(xml)
root = tree.getroot()
for inter_type, interactions in enumerate(root[1][4]):
    for inter in interactions:
        resnr = inter.find("resnr").text
        restype = inter.find("restype").text
```

这里的`root`就是这个文本的根节点了，也就是`<report>`节点。这个根节点上面有两个子节点，分别是`plipversion`和`bindingsite`。`bindingsite`节点上又有许多节点：`identifiers`，`lig_properties`，`interacting_chains`，`bs_residues`，`interactions`。我们需要访问的是`interactions`节点，也就是`root[1][4]`了。`root[1][4]`的意思就是根节点的第2个子节点的第5个子节点，ElementTree是用类似于列表的方式存储的数据。

`interactions`节点中有很多子节点，每个子节点是一种PLIP中定义的相互作用类型，从上到下依次是：

- hydrophobic interactions
- hydrogen bonds
- water bridges
- salt bridges
- pi stacks
- pi cation interactions
- halogen bonds
- metal complexes

对于每种相互作用类型的节点来说，其中又包含多个具体的相互作用，如这里的疏水相互作用类型的父节点中包含一个疏水相互作用的子节点。

```xml
<hydrophobic_interactions>
  <hydrophobic_interaction id="1">
    <resnr>316</resnr>
    <restype>THR</restype>
    <reschain>A</reschain>
  </hydrophobic_interaction>
</hydrophobic_interactions>
```

上述代码中的两层循环的意义就是打开这两层相互作用的节点，进入到具体的一个相互作用的节点中去获取信息。`inter`变量就是一个具体的相互作用了。

PLIP中每种相互作用至少包含这么两个信息：残基编号`resnr`和残基名字`restype`，我们使用find函数找到这两个节点，并使用text属性得到其中的节点信息。当然也可以通过索引去访问，都是一样的。

如此，我们就得到了每个相互作用的残基信息。

相互作用类型的信息也是容易拿到的，注意到相互作用类型的节点在xml文件中的顺序是一样的，也即其索引就可以指代相互作用了。那我们可以直接存储这里的`inter_type`来表示相互作用类型的信息。当然，为了绘图时候呈现图例，我们之后需要将索引转化为对应的相互作用类型名称，不妨定义一个字典：

```python
type_dic = {
    0: "hydrophobic interactions",
    1: "hydrogen bonds",
    2: "water bridges",
    3: "salt bridges",
    4: "pi stacks",
    5: "pi cation",
    6: "halogen bonds",
    7: "metal complexes",
}
```

X轴的信息，其实可以自己定义，当然这里是跟我们输入的xml文件相关的，不妨就用xml文件的序号来作为x轴了，不很要紧的。

整理一下并保存数据：

```python
xmls = [f"conf{i}/report.xml" for i in range(41)]
data = []
type_dic = {
    0: "hydrophobic interactions",
    1: "hydrogen bonds",
    2: "water bridges",
    3: "salt bridges",
    4: "pi stacks",
    5: "pi cation",
    6: "halogen bonds",
    7: "metal complexes",
}

for x, xml in enumerate(xmls):
    tree = ET.parse(xml)
    root = tree.getroot()
    for inter_type, interactions in enumerate(root[1][4]):
        for inter in interactions:
            resnr = inter.find("resnr").text
            restype = inter.find("restype").text
            data.append((x, inter_type, f"{restype}{resnr}"))
```

这里我把x，相互作用类型以及残基信息用tuple存储在一起，data这个列表里面存储的就是多个相互作用啦。

#### 开始绘图

示例的图片中，同一Y轴上有相同的字段，如N51，但又通过不同的散点类型来表示不同的相互作用。所以实际上，Y轴呈现出来的残基信息与相互作用类型是有绑定的。也就是说，这样的一幅散点图，每个点的Y的数据实际上是残基信息与相互作用类型的合，然后再根据相互作用类型去区分点的类型和着色的。

因而，这里我们需要构造各个相互作用（散点）的Y数据：

```python
y_label = []
for d in data:
    if f"{d[2]}-{d[1]}" not in y_label:
        y_label.append(f"{d[2]}-{d[1]}")
```

对于每一个相互作用数据，我们将残基信息与相互作用类型组合了一下。

这里得到的数据就是”真实的“相互作用点的Y轴数据了，但是直接绘制字符串的Y不大行。我们可以考虑对用Y轴数据的索引来生成新的散点：

```python
x_list, y_list = [], []
for d in data:
	x_list.append(d[0])
    y_list.append(y_label.index(f"{d[2]}-{d[1]}"))
```

对于每一个相互作用的点，我们得到了x值与新的y值，这里的`y_label.index(f"{d[2]}-{d[1]}")`获得的就是这一个点的相互作用类型与残基信息在我们设置的`y_label`中的索引，也就是新的y值。之后绘制Y轴的时候，我们把`y_label`按索引贴上去就对应上了。

在绘制的时候，我们还需要根据相互作用类型的不同来绘制不同样式的点，因而我们直接根据相互作用类型来做一下筛选，同样的相互作用类型的散点直接同一批次绘制，就不用每次都自己去指定点的样式了。

绘制散点图这里采用matplotlib进行。

```python
for key in type_dic.keys():
    x_list, y_list = [], []
    for d in data:
        if d[1] == key:
            x_list.append(d[0])
            y_list.append(y_label.index(f"{d[2]}-{d[1]}"))
    plt.scatter(x_list, y_list, label=type_dic[key])
```

我们遍历每一种相互作用类型，如果数据data列表中记录的相互作用类型符合，则将这一个点的x值和y值存一下，之后统一用plt.scatter绘制散点，图例就为相互作用类型的文本。

Y轴数据的绘制不需要呈现出相互作用类型。我们之前组合的时候加了一个`-`，现在刚好根据这个对`y_label`做处理，生成新的Y轴数据：

```python
y_label2show = [item.split("-")[0] for item in y_label]
plt.yticks([i for i in range(len(y_label))], y_label2show)
```

我们把相互作用类型切掉了，然后按照索引把Y轴数据显示到了Y轴上。

最后收个尾：

```python
plt.xlabel("frame")
# 把图例显示在外边
plt.legend(bbox_to_anchor=(1,1))
plt.tight_layout()
plt.show()
```

这样虽然已经满足要求了，但是得到的Y轴的数据可能并不是按照我们期望的顺序进行的排列，例如我们可能希望Y轴数据按照残基编号从小到大显示，那我们还需要对`y_label`排一下序：

```python
y_label = sorted(y_label, key=lambda x:int(x[3:-2]))
```

我们要对残基序号排序，因而这里先把残基序号切出来，然后转成数字，最后用sorted函数进行排序。

贴一张成品图：

![Figure_1](C:\Users\hhhhh\Desktop\20230412-PLIP绘图\Figure_1.png)

附上所有的代码：

```python
## author: charlie
## date : 20230409

import xml.etree.ElementTree as ET
from matplotlib import pyplot as plt

plt.style.use("DIT.mplstyle")

def main():
    xmls = [f"conf{i}/report.xml" for i in range(41)]
    data = []
    type_dic = {
        0: "hydrophobic interactions",
        1: "hydrogen bonds",
        2: "water bridges",
        3: "salt bridges",
        4: "pi stacks",
        5: "pi cation",
        6: "halogen bonds",
        7: "metal complexes",
    }

    for frame, xml in enumerate(xmls):
        tree = ET.parse(xml)
        root = tree.getroot()
        for inter_type, interactions in enumerate(root[1][4]):
            for inter in interactions:
                resnr = inter.find("resnr").text
                restype = inter.find("restype").text
                data.append((frame, inter_type, f"{restype}{resnr}"))

    y_label = []
    for d in data:
        if f"{d[2]}-{d[1]}" not in y_label:
            y_label.append(f"{d[2]}-{d[1]}")
    y_label = sorted(y_label, key=lambda x:int(x.split("-")[0][3:]))
    for key in type_dic.keys():
        x_list, y_list = [], []
        for d in data:
            if d[1] == key:
                x_list.append(d[0])
                y_list.append(y_label.index(f"{d[2]}-{d[1]}"))
        plt.scatter(x_list, y_list, label=type_dic[key])
    y_label2show = [item.split("-")[0] for item in y_label]
    plt.yticks([i for i in range(len(y_label))], y_label2show)
    plt.xlabel("frame")
    plt.legend(bbox_to_anchor=(1,1))
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
```



#### 细节考虑

这个绘图的主要逻辑我们都已经完成了，还剩下一些绘图微调的工作，比如说点的样式、数据间距等。如果需要不同相互作用类型应用不同的点，那可以定义一个样式的列表，然后每次绘图的时候自动赋予就行了：

```python
markers = ["o", "v", "*", ">", "<", "1", "2", "^"]
# ...... some other lines
plt.scatter(x_list, y_list, label=type_dic[key], marker=markers[int(key)])
```

![3](C:\Users\hhhhh\Desktop\20230412-PLIP绘图\3.png)

有一些小**问题需要注意**：我们绘制这个散点图的时候，只考虑了残基信息和相互作用类型，没有考虑配体的相关信息，假如说同一残基与配体的不同原子发生多个了同一类型的相互作用，则在此图上，会表现为**多个彼此完全重叠的点**。另外这种散点的表示形式，散点之间的彼此重叠以及散点本身大小的因素，可能会导致看图者高估相互作用的数量，尤其是X轴数据比较密集的时候。或许热力图会更加合适？虽然也会存在类似的高估问题，但是可能呈现上会更舒服。自然，热力图的绘制可能不太一样了。

上图其实可以看到Y轴上有多个相同的标签，比如HIS227。一种相互作用类型一条横线其实不太简洁。常见的相互作用类型其实可能就两三种，比如说疏水、氢键和pi stacks。如果就三种的话，是不是可以考虑绘制在同一条Y值上，不同类型用不同的点的样式来代表。



### Others

加油







