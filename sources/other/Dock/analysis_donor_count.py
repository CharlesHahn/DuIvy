# author : charlie

import os
import time


def write_aminos_log(line):
    with open("bond_per_amino.aminolog", 'a') as fo:
        fo.write(line+'\n')


def write_ligands_log(line):
    with open("bond_per_ligand.ligandlog", 'a') as fo:
        fo.write(line+'\n')


def load_data(filename):
    with open(filename, 'r') as fo:
        file_content = fo.read()
    data = []
    file_lines = file_content.strip().strip('\n').split('\n')
    for line in file_lines:
        data.append(line.split(','))
    return data


def average_per_amino(data, file_name):
    x_names = ['17LEU', '18VAL', '19PHE', '20PHE', '21ALA', '22GLU', '23ASP', '24VAL', '25GLY', '26SER',
               '27ASN', '28LYS', '29GLY', '30ALA', '31ILE', '32ILE', '33GLY', '34LEU', '35MET', '36VAL',
               '37GLY', '38GLY', '39VAL', '40VAL', '41ILE', '42ALA']
    aminos = {}
    for line in data[1:]:
        key = line[1] + line[2]
        keys_list = aminos.keys()
        if key not in keys_list:
            aminos[key] = [1] + [float(i) for i in line[6:]]

        elif key in keys_list:
            value = aminos[key]
            final_value = []
            final_value.append(value[0] + 1)
            if len(value[1:]) == len(line[6:]):
                for i in range(1, len(value)):
                    final_value.append(value[i] + float(line[i+5]))
                aminos[key] = final_value
            else:
                print(
                    "length of value[0:] if not equal to length of line[6:] ")
                write_aminos_log(
                    "length of value[0:] if not equal to length of line[6:] ")
        else:
            print("WRONG in average_per_amino function ! ")
            write_aminos_log("WRONG in average_per_amino function ! ")

    # print("amino -> [ count , "+ ", ".join(data[0][4:]) + " ]")
    write_aminos_log("amino, count , " + ", ".join(data[0][6:]) )

    for key, value in aminos.items():
        aminos[key] = [aminos[key][0]] + \
            [float(i)/int(aminos[key][0]) for i in aminos[key][1:]]
        # print(key, end = ' -> ')
        # print(["{:.4f}".format(i) for i in aminos[key]])

        # write_aminos_log(file_name +', ' + key + ', ' + ',  '.join(['{:.4f}'.format(i) for i in aminos[key]]))
    for x in x_names:
        if x in aminos.keys(): 
            write_aminos_log( file_name +', ' + x + ', ' + ',  '.join(['{:.4f}'.format(i) for i in aminos[x]]) )
        elif x not in aminos.keys():
            write_aminos_log(file_name + ', ' + x + ', 0')



def average_per_ligand(data, filename):
    ligands = {}
    for line in data[1:]:
        key = line[4]
        keys_list = ligands.keys()
        if key not in keys_list:
            ligands[key] = [1] + [float(i) for i in line[6:]] + [1 if 'donor' in line[5] else 0 ]

        elif key in keys_list:
            value = ligands[key]
            final_value = []
            final_value.append(value[0] + 1)
            if len(value[1:-1]) == len(line[6:]):
                for i in range(1, len(value)-1):
                    final_value.append(value[i] + float(line[i+5]))
                ligands[key] = final_value + [ligands[key][-1] +1 if 'donor' in line[5] else ligands[key][-1]]
            else:
                print(
                    "length of value[0:] if not equal to length of line[6:] ")
                write_ligands_log(
                    "length of value[0:] if not equal to length of line[6:] ")
        else:
            print("WRONG in average_per_ligand function ! ")
            write_ligands_log("WRONG in average_per_ligand function ! ")

    # print("ligand -> [ count , "+ ", ".join(data[0][4:]) + " ]")
    write_ligands_log("ligand, count, " + ", ".join(data[0][6:]) + ", donor")

    for key, value in ligands.items():
        ligands[key] = [ligands[key][0]] + \
            [float(i)/int(ligands[key][0]) for i in ligands[key][1:-1]] + [ligands[key][-1]]
        # print(key, end = ' -> ')
        # print(["{:.4f}".format(i) for i in ligands[key]])
        write_ligands_log(
            filename + ', ' + key + ', ' + ',  '.join(['{:.4f}'.format(i) for i in ligands[key]]))



def main():
    files = os.listdir()
    working_files = []
    for file in files:
        if 'originlog' in file:
            working_files.append(file)
    for working_file in working_files:
        data = load_data(working_file)
        average_per_amino(data, working_file.split('.')[0] )
        average_per_ligand(data, working_file.split('.')[0] )





if __name__ == '__main__':
    main()
