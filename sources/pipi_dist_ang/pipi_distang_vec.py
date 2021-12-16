# author : charlie
# date : 20211129

import os
import math
import argparse


def calcDist(ring_1_frames, ring_2_frames):
    ## check the frames of two ring
    if len(ring_1_frames) != len(ring_2_frames):
        print("Error -> length of frames of coordinates isn't equal")
        exit()
    ## calculate the distance between two ring
    distance = []
    for i in range(len(ring_1_frames)):
        ring_1_coor = ring_1_frames[i]
        ring_2_coor = ring_2_frames[i]
        ring_1_atom_num = len(ring_1_coor)*1.0
        ring_2_atom_num = len(ring_2_coor)*1.0
        ## calculate the center of two ring
        ring_1_center = [ sum([coor[0] for coor in ring_1_coor]) / ring_1_atom_num, 
                sum([coor[1] for coor in ring_1_coor ]) / ring_1_atom_num,
                sum([coor[2] for coor in ring_1_coor ]) / ring_1_atom_num ]
        ring_2_center = [ sum([coor[0] for coor in ring_2_coor]) / ring_2_atom_num, 
                sum([coor[1] for coor in ring_2_coor ]) / ring_2_atom_num,
                sum([coor[2] for coor in ring_2_coor ]) / ring_2_atom_num ]
        ## calculate the distance
        dist = ((ring_1_center[0] - ring_2_center[0])**2 + 
                (ring_1_center[1] - ring_2_center[1])**2 + 
                (ring_1_center[2] - ring_2_center[2])**2 )**0.5
        distance.append(dist)

    return distance


def calcAng(ring_1_frames, ring_2_frames):
    ## check the frames of two ring
    if len(ring_1_frames) != len(ring_2_frames):
        print("Error -> length of frames of coordinates isn't equal")
        exit()
    ## calculate the angles
    angles = []
    for i in range(len(ring_1_frames)):
        ring_1_coor = ring_1_frames[i]
        ring_2_coor = ring_2_frames[i]
        ## p1 = [ a1, a2, a3 ], p2 = [ b1, b2, b3 ]
        ## p1xp2 = [ a2b3 - a3b2, a3b1 - a1b3, a1b2 - a2b1 ]
        ## calculate the plane normal of ring 1 by cross product
        r1_p1 = [ ring_1_coor[2][0] - ring_1_coor[0][0], 
                ring_1_coor[2][1] - ring_1_coor[0][1],
                ring_1_coor[2][2] - ring_1_coor[0][2]]
        r1_p2 = [ ring_1_coor[4][0] - ring_1_coor[0][0], 
                ring_1_coor[4][1] - ring_1_coor[0][1],
                ring_1_coor[4][2] - ring_1_coor[0][2]]
        r1_xProd = [ r1_p1[1]*r1_p2[2] - r1_p1[2]*r1_p2[1], 
                r1_p1[2]*r1_p2[0] - r1_p1[0]*r1_p2[2], 
                r1_p1[0]*r1_p2[1] - r1_p1[1]*r1_p2[0] ] 
        ## calculate the plane normal of ring 2 by cross product
        r2_p1 = [ ring_2_coor[2][0] - ring_2_coor[0][0], 
                ring_2_coor[2][1] - ring_2_coor[0][1],
                ring_2_coor[2][2] - ring_2_coor[0][2]]
        r2_p2 = [ ring_2_coor[4][0] - ring_2_coor[0][0], 
                ring_2_coor[4][1] - ring_2_coor[0][1],
                ring_2_coor[4][2] - ring_2_coor[0][2]]
        r2_xProd = [ r2_p1[1]*r2_p2[2] - r2_p1[2]*r2_p2[1], 
                r2_p1[2]*r2_p2[0] - r2_p1[0]*r2_p2[2], 
                r2_p1[0]*r2_p2[1] - r2_p1[1]*r2_p2[0] ] 
        ## calc the degree 
        dotProduct = r1_xProd[0] * r2_xProd[0] + r1_xProd[1] * r2_xProd[1]
        dotProduct += r1_xProd[2] * r2_xProd[2]
        r1_xProd_norm = (r1_xProd[0]**2 + r1_xProd[1]**2 + r1_xProd[2]**2)**0.5
        r2_xProd_norm = (r2_xProd[0]**2 + r2_xProd[1]**2 + r2_xProd[2]**2)**0.5
        cos_degree = dotProduct / (r1_xProd_norm * r2_xProd_norm)
        degree = math.acos(cos_degree) * 180 / math.pi 
        if degree > 90 :
            degree = 180 - degree
        angles.append(degree)

    return angles


def getCoor(gro_file, ring_1_id, ring_2_id):
    with open(gro_file, "r") as fo: 
        lines = fo.readlines()
    frames = []
    atom_lines = []
    for line in lines:
        ## not the atom line, judged by length of line
        if len(line) != 45 and len(line) != 69:
            # number line
            if len(line.strip().split()) == 1:
                if len(atom_lines) != 0:
                    frames.append(atom_lines)
                atom_lines = []
        elif len(line) == 45 or len(line) == 69:
            # suit for atom number < 10000
            if len(line[20:44].split()) == 3 and len(line[15:44].split()) == 4:
                atom_lines.append(line.rstrip())
    ## add the last frame
    if len(atom_lines) != 0:
        frames.append(atom_lines)
    ## read the atom coor
    ring_1_frames = []
    ring_2_frames = []
    for frame in frames:
        ring_1_coor = []
        ring_2_coor = []
        for line in frame:
            if int(line[15:20]) in ring_1_id:
                ring_1_coor.append([
                    float(line[20:28]), float(line[28:36]), float(line[36:44]) ])
            if int(line[15:20]) in ring_2_id:
                ring_2_coor.append([
                    float(line[20:28]), float(line[28:36]), float(line[36:44]) ])
        if len(ring_1_coor) != len(ring_1_id) or len(ring_2_coor) != len(ring_2_id):
            print("Error -> shit happens when reading coordinates of ring")
            exit()
        ring_1_frames.append(ring_1_coor)
        ring_2_frames.append(ring_2_coor)
    ## new a time sequence
    time = [ i for i in range(len(frames)) ]            

    return time, ring_1_frames, ring_2_frames


def dealNdx(ndx_file, select, vg=False):
    with open(ndx_file, 'r') as fo:
        content = fo.read()
    ## read in each group
    ndx_dic = {}
    content = content.strip().strip("[").replace("\n", " ")
    ndx_groups = content.split("[")
    for group in ndx_groups:
        items = group.split("]")
        if items[0].strip() in ndx_dic.keys():
            print("Error -> two groups with the same name")
            exit()
        ndx_dic[items[0].strip()] = [ int(n.strip()) 
                for n in items[1].split() if n != "" ]
    if select == None:
        ## print to get input from user
        print("Info -> reading your index file:")
        for name, num_lis in ndx_dic.items():
            print("    {:30}   {:8} atoms".format(name, len(num_lis)))
            # print(" ".join([ str(i) for i in num_lis]))
        if vg == False:
            ring_1_name = input("Type the name of first ring -> ")
            print("Info -> you have chosed " + ring_1_name + " as first ring.")
            ring_2_name = input("Type the name of second ring -> ")
            print("Info -> you have chosed " + ring_2_name + " as second ring.")
        elif vg == True:
            ring_1_name = input("Type the name of ring group -> ")
            print("Info -> you have chosed " + ring_1_name + " as first ring.")
            ring_2_name = input("Type the name of vector group-> ")
            print("Info -> you have chosed " + ring_2_name + " as vector group.")
    elif len(select) == 2:
        if vg == False:
            ring_1_name = select[0]
            print("Info -> you have chosed " + ring_1_name + " as first ring.")
            ring_2_name = select[1]
            print("Info -> you have chosed " + ring_2_name + " as second ring.")
        elif vg == True:
            ring_1_name = select[0]
            print("Info -> you have chosed " + ring_1_name + " as first ring.")
            ring_2_name = select[1]
            print("Info -> you have chosed " + ring_2_name + " as vector group.")
    else:
        print("Error -> Wrong parameter number of -select")
        exit()
        
    ## atom id of two ring groups
    ring_1_id = ndx_dic[ring_1_name]
    ring_2_id = ndx_dic[ring_2_name]

    ## return the atom ids of these two rings
    return ring_1_id, ring_2_id 


def dealNdx_single(ndx_file, select):
    with open(ndx_file, 'r') as fo:
        content = fo.read()
    ## read in each group
    ndx_dic = {}
    content = content.strip().strip("[").replace("\n", " ")
    ndx_groups = content.split("[")
    for group in ndx_groups:
        items = group.split("]")
        if items[0].strip() in ndx_dic.keys():
            print("Error -> two groups with the same name")
            exit()
        ndx_dic[items[0].strip()] = [ int(n.strip()) 
                for n in items[1].split() if n != "" ]
    ## print to get input from user
    print("Info -> reading your index file:")
    for name, num_lis in ndx_dic.items():
        print("    {:30}   {:8} atoms".format(name, len(num_lis)))
    # print(" ".join([ str(i) for i in num_lis]))
    if select == None:
        ring_1_name = input("Type the name of first ring -> ")
        print("Info -> you have chosed " + ring_1_name + " as first ring.")
    elif len(select) == 1:
        ring_1_name = select[0]
        print("Info -> you have chosed " + ring_1_name + " as first ring.")
    else:
        print("Error -> you may only need one parameter for select")
        exit()

    ## atom id of two ring groups
    ring_1_id = ndx_dic[ring_1_name]

    ## return the atom ids of these two rings
    return ring_1_id


def calcAng_RingVec(ring_frames, vec_frames):
    ## check the frames 
    if len(ring_frames) != len(vec_frames):
        print("Error -> length of frames of ring and vector isn't equal")
        exit()
    ## calculate the angles
    angles = []
    for i in range(len(ring_frames)):
        ring_1_coor = ring_frames[i]
        vector = vec_frames[i]
        ## p1 = [ a1, a2, a3 ], p2 = [ b1, b2, b3 ]
        ## p1xp2 = [ a2b3 - a3b2, a3b1 - a1b3, a1b2 - a2b1 ]
        ## calculate the plane normal of ring 1 by cross product
        r1_p1 = [ ring_1_coor[2][0] - ring_1_coor[0][0], 
                ring_1_coor[2][1] - ring_1_coor[0][1],
                ring_1_coor[2][2] - ring_1_coor[0][2]]
        r1_p2 = [ ring_1_coor[4][0] - ring_1_coor[0][0], 
                ring_1_coor[4][1] - ring_1_coor[0][1],
                ring_1_coor[4][2] - ring_1_coor[0][2]]
        r1_xProd = [ r1_p1[1]*r1_p2[2] - r1_p1[2]*r1_p2[1], 
                r1_p1[2]*r1_p2[0] - r1_p1[0]*r1_p2[2], 
                r1_p1[0]*r1_p2[1] - r1_p1[1]*r1_p2[0] ] 
        ## calc the degree 
        dotProduct = r1_xProd[0] * vector[0] + r1_xProd[1] * vector[1]
        dotProduct += r1_xProd[2] * vector[2]
        r1_xProd_norm = (r1_xProd[0]**2 + r1_xProd[1]**2 + r1_xProd[2]**2)**0.5
        vector_norm = (vector[0]**2 + vector[1]**2 + vector[2]**2)**0.5
        cos_degree = dotProduct / (r1_xProd_norm * vector_norm)
        degree = math.acos(cos_degree) * 180 / math.pi 
        if degree > 90 :
            degree = 180 - degree
        angles.append(degree)
    return angles


def calcVec(vg_frames):
    vg_vec_frames = []
    for atoms in vg_frames:
        vec = [ 0, 0, 0 ]
        for i in range(1,len(atoms)):
            vec[0] += atoms[i][0] - atoms[i-1][0]
            vec[1] += atoms[i][1] - atoms[i-1][1]
            vec[2] += atoms[i][2] - atoms[i-1][2]
        vec = [ v/(len(atoms)-1) for v in vec ]
        vg_vec_frames.append(vec)
    return vg_vec_frames


def dealTwoRings(ndx_file, gro_file, time_b, time_dt, output_file, select):
    ## get the atom id
    ring_1_id, ring_2_id = dealNdx(ndx_file, select, False)
    if len(ring_1_id) < 5 or len(ring_1_id) > 7 or \
            len(ring_2_id) < 5 or len(ring_2_id) > 7 : 
        print("Error -> index of your ring is more than 7 or less than 5")
        print("Error -> only support 5, 6 or 7 membered ring which is in a plane ")
        print("Error -> please check your index file")
        print("Error -> your index : ", ring_1_id, ring_2_id)
        exit()
    ## get the coordinates of two rings
    time, ring_1_frames, ring_2_frames = getCoor(gro_file, ring_1_id, ring_2_id)
    ## modify the time sequence
    time = [ t*time_dt + time_b for t in time ]
    # print(len(time), len(ring_1_frames), len(ring_2_frames))
    ## calculate the distance of two rings
    distance = calcDist(ring_1_frames, ring_2_frames)
    ## calculate the angles of the normals of two rings
    angles = calcAng(ring_1_frames, ring_2_frames)
    # print(len(time), len(distance), len(angles))
    ## check data and output
    print("Info -> there is ", len(time), " frames in your gro file")
    if len(time) != len(distance) or len(time) != len(angles):
        print("Error -> length of time, dist, ang are not equal")
        print(len(time), len(distance), len(angles))
    out_content = "{:<10} {:>10} {:>10} \n".format("time", "distance", "angle")
    out_content += "\n".join(["{:<10.0f} {:>10.3f} {:>10.3f} ".format(
        time[i], distance[i], angles[i] ) for i in range(len(distance))])
    with open(output_file, 'w') as fo:
        fo.write(out_content)

    ## calc the angle distribution
    ang0_30, ang30_60, ang60_90 = 0, 0, 0
    for ang in angles:
        if ang >= 0 and ang < 30:
            ang0_30 += 1
        elif ang >= 30 and ang < 60:
            ang30_60 += 1
        elif ang > 60 :
            ang60_90 += 1
    print("Info =>  0 <= angle < 30 : {}/{} = {:>6.2%}".format(
        ang0_30, len(time), ang0_30*1.0/len(time)))
    print("Info => 30 <= angle < 60 : {}/{} = {:>6.2%}".format(
        ang30_60, len(time), ang30_60*1.0/len(time)))
    print("Info => 60 <= angle < 90 : {}/{} = {:>6.2%}".format(
        ang60_90, len(time), ang60_90*1.0/len(time)))
    ## calc the average distance
    print("Info => average distance : {:>10.4f} nm".format(
        sum(distance)/len(time)))


def dealRingVG(ndx_file, gro_file, time_b, time_dt, output_file, select):
    ring_id, vg_id = dealNdx(ndx_file, select, True)
    if len(ring_id) < 5 or len(ring_id) > 7 :
        print("Error -> index of your ring is more than 7 or less than 5")
        print("Error -> only support 5, 6 or 7 membered ring which is in a plane ")
        print("Error -> please check your index file")
        print("Error -> your index : ", ring_id)
        exit()
    if len(vg_id) < 2 :
        print("Error -> less than 2 atom index in vector group you input")
        print("Error -> 2 or more atom index are needed")
        print("Error -> check your vector index :", vg_id)
        exit()
    ## get the coordinates of ring and vg
    time, ring_frames, vg_frames = getCoor(gro_file, ring_id, vg_id)
    ## modify the time sequence
    time = [ t*time_dt + time_b for t in time ]
    # calculate the vector vs time
    vg_vec_frames = calcVec(vg_frames)
    ## calculate the angles 
    angles = calcAng_RingVec(ring_frames, vg_vec_frames)
    # save results
    print("Info -> there is ", len(time), " frames in your gro file")
    if len(time) != len(angles):
        print("Error -> length of time, angles are not equal")
        print(len(time), len(angles))
    out_content = "{:<10} {:>10} \n".format("time", "angle")
    out_content += "\n".join(["{:<10.0f} {:>10.3f} ".format(
        time[i], angles[i] ) for i in range(len(angles))])
    with open(output_file, 'w') as fo:
        fo.write(out_content)
    ## calc the angle distribution
    ang0_30, ang30_60, ang60_90 = 0, 0, 0
    for ang in angles:
        if ang >= 0 and ang < 30:
            ang0_30 += 1
        elif ang >= 30 and ang < 60:
            ang30_60 += 1
        elif ang > 60 :
            ang60_90 += 1
    print("Info =>  0 <= angle < 30 : {}/{} = {:>6.2%}".format(
        ang0_30, len(time), ang0_30*1.0/len(time)))
    print("Info => 30 <= angle < 60 : {}/{} = {:>6.2%}".format(
        ang30_60, len(time), ang30_60*1.0/len(time)))
    print("Info => 60 <= angle < 90 : {}/{} = {:>6.2%}".format(
        ang60_90, len(time), ang60_90*1.0/len(time)))


def dealRingVec(ndx_file, gro_file, time_b, time_dt, output_file,
                vec, select):
    ring_id = dealNdx_single(ndx_file, select)
    if len(ring_id) < 5 or len(ring_id) > 7 :
        print("Error -> index of your ring is more than 7 or less than 5")
        print("Error -> only support 5, 6 or 7 membered ring which is in a plane ")
        print("Error -> please check your index file")
        print("Error -> your index : ", ring_id)
        exit()
    time, ring_frames, _ = getCoor(gro_file, ring_id, ring_id)
    ## modify the time sequence
    time = [ t*time_dt + time_b for t in time ]
    # deal with vec 
    vec = [ float(i) for i in vec ]
    if 0 == vec[0] == vec[1] == vec[2]:
        print("Error -> You can't input an all zero vector")
        exit()
    vec_frames = [ vec for i in range(len(time)) ]
    # calculate the angles
    angles = calcAng_RingVec(ring_frames, vec_frames)
    # save results
    print("Info -> there is ", len(time), " frames in your gro file")
    if len(time) != len(angles):
        print("Error -> length of time, angles are not equal")
        print(len(time), len(angles))
    out_content = "{:<10} {:>10} \n".format("time", "angle")
    out_content += "\n".join(["{:<10.0f} {:>10.3f} ".format(
        time[i], angles[i] ) for i in range(len(angles))])
    with open(output_file, 'w') as fo:
        fo.write(out_content)
    ## calc the angle distribution
    ang0_30, ang30_60, ang60_90 = 0, 0, 0
    for ang in angles:
        if ang >= 0 and ang < 30:
            ang0_30 += 1
        elif ang >= 30 and ang < 60:
            ang30_60 += 1
        elif ang > 60 :
            ang60_90 += 1
    print("Info =>  0 <= angle < 30 : {}/{} = {:>6.2%}".format(
        ang0_30, len(time), ang0_30*1.0/len(time)))
    print("Info => 30 <= angle < 60 : {}/{} = {:>6.2%}".format(
        ang30_60, len(time), ang30_60*1.0/len(time)))
    print("Info => 60 <= angle < 90 : {}/{} = {:>6.2%}".format(
        ang60_90, len(time), ang60_90*1.0/len(time)))


def main():
    ## parse the input argvs
    print("Info -> To calculate the distance and angles", end = "")
    print(" of two 5 or 6 membered rings")
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", help = "index file contains two groups of rings")
    parser.add_argument("-f", help = "gro file")
    parser.add_argument("-b", default=0, type=int, 
            help = "set the start time, default=0")
    parser.add_argument("-dt", default=1, type=int, 
            help = "set the time interval, default=1")
    parser.add_argument("-o", default="output.xvg", 
            help = "the results data, default output.xvg")
    parser.add_argument("-vg", action="store_true",
            help = "whether to get vector by index group")
    parser.add_argument("-vec", nargs=3, 
            help = "get vector by your input, eg. -vec 6 6 6")
    parser.add_argument("-select", nargs="*",
            help = "select the groups, eg. -select ring1 ring2")
    args = parser.parse_args()
    ndx_file = args.n
    gro_file = args.f
    time_b = args.b
    time_dt = args.dt
    output_file = args.o
    vg = args.vg
    vec = args.vec
    select = args.select

    if ndx_file not in os.listdir():
        print("Error -> no ", ndx_file, " in current directory")
        exit()
    if gro_file not in os.listdir():
        print("Error -> no ", gro_file, " in current directory")
        exit()

    if vec == None and vg == False:
        dealTwoRings(ndx_file, gro_file, time_b, time_dt, output_file, select)
    elif vec == None and vg == True:
        dealRingVG(ndx_file, gro_file, time_b, time_dt, output_file, select)
    elif vec != None and vg == False:
        dealRingVec(ndx_file, gro_file, time_b, time_dt, output_file,
                    vec, select)
    elif vec != None and vg == True:
        print("Error -> You can't set -vg and -vec at the same time")
        exit()

    print("Done -> good day ! ")


if __name__ == "__main__":
    main()

