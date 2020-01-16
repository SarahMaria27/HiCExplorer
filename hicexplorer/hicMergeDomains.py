#!/usr/bin/env python
import argparse


def parse_arguments(args=None):
    """
    get command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        conflict_handler='resolve',
        description="""""")

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--domain1', '-d1',
                                help='The domains.bed file of the first matrix is required',
                                required=True)

    parserRequired.add_argument('--domain2', '-d2',
                                help='The domains.bed file of the second: matrix is required',
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--value', '-v',
                           help='Determine a value by how much the boundaries of two TADs must at least differ to view them as two separate areas.',
                           type=int, default=0)
    return parser

def create_list_of_file(file):
    with open(file) as f:
        newList = [line.rstrip() for line in f]
    splittedList = []
    for line in newList:
        x = line.split("\t")
        splittedList.append(x[0:3])
    return splittedList


def merge_list(d1, d2, pValue):
    pos1, pos2 = 0, 0
    merged_list = []
    while (pos1 < len(d1) and pos2 < len(d2)):
        while (d1[pos1][0] == d2[pos2][0]):
            if (int(d1[pos1][1]) <= int(d2[pos2][1])):
                if ((abs(int(d1[pos1][1])-int(d2[pos2][1]))>pValue) or ((abs(int(d1[pos1][2])-int(d2[pos2][2]))>pValue))):
                    merged_list.append(d1[pos1])
                if (pos1+1 != len(d1)): pos1 += 1
                else: break
            elif (int(d1[pos1][1]) > int(d2[pos2][1])):
                if ((abs(int(d1[pos1][1])-int(d2[pos2][1]))>pValue) or ((abs(int(d1[pos1][2])-int(d2[pos2][2]))>pValue))):
                    merged_list.append(d2[pos2])
                if (pos2+1 != len(d2)): pos2 += 1
                else: break
        if (d1[pos1-1][0] == d2[pos2][0]):
            while ((pos2 < len(d2)) and (d1[pos1-1][0] == d2[pos2][0])):
                merged_list.append(d2[pos2])
                pos2 += 1
        if (d1[pos1][0] == d2[pos2-1][0]):
            while ((pos1 < len(d1))and (d1[pos1][0] == d2[pos2-1][0])):
                merged_list.append(d1[pos1])
                pos1 += 1
    return merged_list
            

def add_id(domainList):
    id_number = 1
    for tad in domainList:
        tad.append("ID_" + str(id_number))
        id_number += 1
    return domainList
        

def write_in_file(l):
    filename = "mergedDomains.bed"
    myfile = open(filename, 'w')
    i = 0
    while (i < len(l)):
        myfile.write(l[i][0]+'\t'+l[i][1]+'\t'+l[i][2]+'\t'+l[i][3]+'\n')
        i += 1
    myfile.close()


def main(args=None):

    args = parse_arguments().parse_args(args)
    pValue = args.value
    domainList1= create_list_of_file(args.domain1)
    domainList2 = create_list_of_file(args.domain2)
    merged_list = merge_list(domainList1, domainList2, pValue)
    merged_list_with_id = add_id(merged_list)
    write_in_file(merged_list_with_id)

main()


