#!/usr/bin/env python
import argparse
from scipy.cluster.hierarchy import dendrogram, linkage
import numpy as np
from matplotlib import pyplot as plt


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

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--domainList', '-l',
                                help='The domains.bed file of the second: matrix is required',
                           nargs='+')

    parserOpt.add_argument('--value', '-v',
                           help='Determine a value by how much the boundaries of two TADs must at least differ to view them as two separate areas.',
                           type=int, default=0)

    parserOpt.add_argument('--percent', '-p',
                           help='For the relationship determination, a percentage is required from which area coverage the TADs are related to each other.'
			   'For example, a relationship should be entered from 5 percent area coverage -p 0.05',
                           type=float, default=0.5)

    return parser

def create_list_of_file(file):
    with open(file) as f:
        newList = [line.rstrip() for line in f]
    splittedList = []
    for line in newList:
        x = line.split("\t")
        splittedList.append(x)
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
        tad[3] = "ID_" + str(id_number)
        id_number += 1
    return domainList


def create_relationsship_list(domainList, pPercent):
    relationList = []
    tad1, tad2 = 0, 1
    while (tad1 < len(domainList)):
        while (tad2 < len(domainList)):
            if (int(domainList[tad1][2]) < int(domainList[tad2][1])) or (domainList[tad1][0] != domainList[tad2][0]):
                break
            halfTad1 = (float(domainList[tad1][2])-float(domainList[tad1][1]))*pPercent
            halfTad2 = (float(domainList[tad2][2])-float(domainList[tad2][1]))*pPercent
            # Überprüfung auf Überlappende TAD´s
            if (((float(domainList[tad1][2])-halfTad1) > float(domainList[tad2][1])) and ((float(domainList[tad1][2])+halfTad1) <= float(domainList[tad2][2]))):
                if ((float(domainList[tad1][2])-float(domainList[tad1][1])) > (float(domainList[tad2][2])-int(domainList[tad2][1]))):
                    relationList.append([domainList[tad1][0], domainList[tad1][3], domainList[tad2][3]])
                else:
                    relationList.append([domainList[tad2][0], domainList[tad2][3], domainList[tad1][3]])
            # Überprüfung auf Innenliegende TAD´s
            elif (((float(domainList[tad1][1])) <= float(domainList[tad2][1])) and (float(domainList[tad1][2]) >= int(domainList[tad2][2]))):
                relationList.append([domainList[tad1][0], domainList[tad1][3], domainList[tad2][3]])
            tad2 += 1
        tad2 = tad1 + 2
        tad1 += 1
    return relationList

        

def write_in_file(l, name):
    filename = name
    myfile = open(filename, 'w')
    i = 0
    while (i < len(l)):
        element = 0
        string = ""
        while (element < len(l[i])-1):
            string += l[i][element] + '\t'
            element += 1
        string += l[i][element] + '\n'
        myfile.write(string)
        i += 1
    myfile.close()


def main(args=None):

    args = parse_arguments().parse_args(args)
    pValue = args.value
    mergedList= create_list_of_file(args.domain1)
    listOfDomains = []
    for domain in args.domainList:
        listOfDomains.append(create_list_of_file(domain))    
    for domain in listOfDomains:
        mergedList = merge_list(mergedList, domain, pValue)
    mergedListWithId = add_id(mergedList)
    write_in_file(mergedListWithId, "mergedDomains.bed")
    relationList = create_relationsship_list(mergedListWithId, args.percent)
    write_in_file(relationList, "relationList.bed")



