# 读入npt.gro 
with open("npt.gro", "r") as fo:
    lines = fo.readlines()
# 定义几何中心的变量
center_X = 0 
center_Y = 0
center_Z = 0 
# 遍历每一个蛋白的原子，计算几何中心
pro_start=3
pro_num = 5700
# 因为python计数从0开始，所以这里pro_start-1
for line in lines[pro_start-1:pro_start+pro_num-1]:
    # 按空格对每一行进行切分
    items = line[20:44].split()
    # 蛋白原子的坐标在第4、5、6列
    center_X += float(items[0])
    center_Y += float(items[1])
    center_Z += float(items[2])
# 计算几何中心的坐标
center_X = center_X/pro_num
center_Y = center_Y/pro_num
center_Z = center_Z/pro_num
# 定义变量存储寻找到的最近原子的信息
# atom_info记录寻找到的原子序号
atom_info = ""
# dist记录寻找到的原子与几何中心的距离
dist = 10
# 再遍历一遍原子，寻找离几何中心最近的原子
for line in lines[pro_start-1:pro_start+pro_num-1]:
    items = line[20:44].split()
    x = float(items[0])
    y = float(items[1])
    z = float(items[2])
    # 计算距离，小于dist就刷新atom_info和dist
    atom_center_dist = ((x-center_X)**2 + (y-center_Y)**2 + (z-center_Z)**2 )**0.5
    if atom_center_dist < dist:
        dist = atom_center_dist
        # 记录此原子信息
        atom_info = line
# 输出中心原子的信息
print(atom_info)
