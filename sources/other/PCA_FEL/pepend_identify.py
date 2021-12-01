# author : charlie

import sys
import getopt 
import time
import os


def readpdb(inputfile):
    if inputfile not in os.listdir():
        print('Error => {} is not in this directory !'.format(inputfile))
        sys.exit(0)
    with open(inputfile, 'r', encoding='utf-8') as fo:
        content = fo.read()
    return content


def writepdb(outputfile, content2write, write_flag):
    if write_flag == 0 and outputfile in os.listdir():
        print('Error => {} is already exist in this directory !'.format(
            outputfile))
        sys.exit(0)
    with open(outputfile, 'a', encoding='utf-8') as fo:
        fo.write(content2write)


def distinguish(content, inputfile, outputfile):
    model_count = 0
    write_flag = 0
    models = content.strip().strip('\n').strip('ENDMDL').split('ENDMDL')
    len_models = len(models)
    print('Detect -> {} MODELs in {}'.format(len_models, inputfile))
    for model in models:
        content2write = ''
        identifier = ''
        for i in ['A','B','C','D','E','F','G','H','I','J']:
            for j in range(1140):
                identifier += i
        # print(len(identifier))
        lines = model.strip('\n').split('\n')
        # print("number of lines -> {}".format(len(lines)))
        for line in lines:
            if 'ATOM' in line:
                atom_num = int(line[5:11].strip())
                if atom_num <= 11400:
                     # remind u: atom_num is start from 1 instead of naught
                     line2write = line[:21]+identifier[atom_num-1]+line[22:]+'\n'
                else:
                    line2write = line + '\n'
            else:
                line2write = line + '\n'
            content2write += line2write
        # add ENDMDL keyword
        content2write += 'ENDMDL\n'
        # judge whether the writing should end
        writepdb(outputfile, content2write, write_flag)
        # set write_flag = 1 means the successive writing shall continue
        write_flag = 1
        model_count += 1
        print("Progress -> |{:-<50}|".format(
            int((model_count/len_models)*50) * '>' ))

def main():
    time_start = time.perf_counter()
    print("Help => python pepend_identify.py -i inputfile -o outputfile")
    inputfile, outputfile = '','outputfile_from_pepend_identify.pdb'
    try:
        opts,args = getopt.getopt(sys.argv[1:],'i:o:')
    except:
        print("Error => wrong arguments, input again ! ")
        sys.exit(0)
    else:
        for opt, arg in opts:
            if opt == '-i':
                inputfile = arg
            elif opt == '-o':
                outputfile = arg 
    if inputfile == '':
        print("Error => no input file !")
        sys.exit(0)
    if outputfile == 'outputfile_from_pepend_identify.pdb':
        print("info -> not specify output filename, use default")
    content = readpdb(inputfile)
    distinguish(content, inputfile, outputfile)
    time_end = time.perf_counter()
    print("Done in {:.4f} s -> {} chain identifier successfully added !".format(
        time_end - time_start, inputfile))


if __name__ == "__main__":
    main()
