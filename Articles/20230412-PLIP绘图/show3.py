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
                data.append((frame, inter_type, restype, resnr))

    # for d in data:
    #     print(d)
    y_label = []
    for d in data:
        if (d[1], d[2], d[3]) not in y_label:
            y_label.append((d[1], d[2], d[3]))
    y_label = sorted(y_label, key=lambda x:int(x[2]))
    y_label = [f"{y[1]}{y[2]}-{y[0]}" for y in y_label]
    markers = ["o", "v", "*", ">", "<", "1", "2", "^"]
    for key in type_dic.keys():
        x_list, y_list = [], []
        for d in data:
            if d[1] == key:
                x_list.append(d[0])
                y_list.append(y_label.index(f"{d[2]}{d[3]}-{d[1]}"))
        plt.scatter(x_list, y_list, label=type_dic[key], marker=markers[int(key)])
    y_label2show = [item.split("-")[0] for item in y_label]
    plt.yticks([i for i in range(len(y_label))], y_label2show)
    plt.xlabel("frame")
    plt.legend(bbox_to_anchor=(1,1))
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
