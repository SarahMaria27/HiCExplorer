#!/usr/bin/env python
import argparse
from scipy.cluster.hierarchy import dendrogram, linkage
import numpy as np
from matplotlib import pyplot as plt
from graphviz import Digraph


def parse_arguments(args=None):
    """
    get command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        conflict_handler='resolve',
        description="""""")

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--domain1', '-d',
                                help='The domains.bed file of the first matrix is required',
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--domainList', '-l',
                                help='The domains.bed file of the second: matrix is required',
                           nargs='+')

    parserOpt.add_argument('--ctcfFile', '-c',
                           help='In order to be able to better assess the relationship between TADs, the associated CTCF file can be included.',
                           type=int, default=0)

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


def merge_list(d1, d2, pValue= 5000):
    pos1, pos2 = 0, 0
    merged_list = []
    while (pos1 < len(d1) and pos2 < len(d2)):
        while (d1[pos1][0] == d2[pos2][0]):
            if (int(d1[pos1][1]) <= int(d2[pos2][1])):
                if ((abs(int(d1[pos1][1])-int(d2[pos2][1]))>pValue) or ((abs(int(d1[pos1][2])-int(d2[pos2][2]))>pValue))):
                    merged_list.append(d1[pos1])
                if (pos1+1 != len(d1)):
                    pos1 += 1
                else: 
                    while ((pos2 < len(d2)) and (d1[pos1][0] == d2[pos2][0])):
                        merged_list.append(d2[pos2])
                        pos2 += 1
                    break
            elif (int(d1[pos1][1]) > int(d2[pos2][1])):
                if ((abs(int(d1[pos1][1])-int(d2[pos2][1]))>pValue) or ((abs(int(d1[pos1][2])-int(d2[pos2][2]))>pValue))):
                    merged_list.append(d2[pos2])
                if (pos2+1 != len(d2)):
                    pos2 += 1         
                else: 
                    while ((pos1 < len(d1))and (d1[pos1][0] == d2[pos2][0])):
                        merged_list.append(d1[pos1])
                        pos1 += 1
                    break
        print("Merged", d1[pos1-1][0])
        if (pos1 < len(d1) and pos2 < len(d2) and d1[pos1][0] != d2[pos2][0]):
            if (d1[pos1-1][0] == d2[pos2][0]):
                pos2 += 1
            else:
                pos1 += 1
    return merged_list
            

def add_id(domainList):
    id_number = 1
    for tad in domainList:
        tad[3] = "ID_" + str(id_number)
        id_number += 1
    return domainList


def create_relationsship_list(domainList, pProzent = 0.5):
    relationList = []
    tad1, tad2 = 0, 1
    while (tad1 < len(domainList)):
        while (tad2 < len(domainList)):
            if (int(domainList[tad1][2]) < int(domainList[tad2][1])) or (domainList[tad1][0] != domainList[tad2][0]):
                break
            halfTad1 = (float(domainList[tad1][2])-float(domainList[tad1][1]))*pProzent
            halfTad2 = (float(domainList[tad2][2])-float(domainList[tad2][1]))*pProzent
            # Überprüfung auf Überlappende TAD´s
            if (((float(domainList[tad1][2])-halfTad1) > float(domainList[tad2][1])) and ((float(domainList[tad1][2])+halfTad1) <= float(domainList[tad2][2]))):
                if ((float(domainList[tad1][2])-float(domainList[tad1][1])) > (float(domainList[tad2][2])-int(domainList[tad2][1]))):
                    add_relation_to_list(relationList, domainList[tad1][0], domainList[tad1][3], domainList[tad2][3])
                else:
                    add_relation_to_list(relationList, domainList[tad2][0], domainList[tad2][3], domainList[tad1][3])   
            # Überprüfung auf Innenliegende TAD´s
            elif (((float(domainList[tad1][1])) <= float(domainList[tad2][1])) and (float(domainList[tad1][2]) >= int(domainList[tad2][2]))):
                add_relation_to_list(relationList, domainList[tad1][0], domainList[tad1][3], domainList[tad2][3])
            tad2 += 1
        tad2 = tad1 + 2
        tad1 += 1
    return relationList

def add_relation_to_list(rList, chromosom, parent, child):
    if (len(rList) == 0):
        rList.append([chromosom, parent, [child]])
        return rList
    pos = len(rList)-1
    while(True):
        if (rList[pos][1] == parent):
            rList[pos][2].append(child)
            return rList
        elif (int(rList[pos][1][3:]) < int(parent[3:])):
            rList.append([chromosom, parent, [child]])
            return rList
        else:
            pos -= 1
        

def write_in_file(l, name, relation = False):
    filename = name
    myfile = open(filename, 'w')
    i = 0
    if (relation):
        while (i < len(l)):
            element = 0
            while (element < len(l[i][2])):
                string = l[i][0] + '\t' + l[i][1] + '\t' + l[i][2][element] + '\n'
                myfile.write(string)
                element += 1
            i += 1
        myfile.close()
    else:
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

def create_tree(rList, dList):
    name = rList[0][0] + '_relations'
    g = Digraph(filename=name)
    sList = create_small_list(dList)
    chrom = rList[0][0]
    pos = 0
    while (pos < len(rList)):
        if (chrom == rList[pos][0]):
            g.edge(rList[pos][0], rList[pos][1])
            for child in rList[pos][2]:
                g.edge(rList[pos][1], child)
        else:
            """for single in sList:
                if (rList[pos-1][0] == single[0]):
                    g.edge(single[0], single[1])
                    sList.remove([single[0], single[1]])"""
            print("Saved relation tree of " + chrom)
            g.render()
            chrom = rList[pos][0]
            name = chrom + '_relations'
            g = Digraph(filename=name)
        pos += 1
    print("Saved relation tree of " + chrom)
    g.render()
    """for single in sList:
        if (rList[pos-1][0] == single[0]):
            g.edge(single[0], single[1])
            sList.remove([single[0], single[1]])"""

def create_small_list(dList):
    sList = []
    for tad in dList:
        sList.append([tad[0], tad[3]])
    return sList

def read_ctcf(file):
    with open(file) as f:
        newList = [line.rstrip() for line in f]
    splittedList = [[]]
    actualChr = newList[0].split("\t")[0]
    for line in newList:
        x = line.split("\t")
        if (x[0] == actualChr):
            splittedList[len(splittedList)-1].append(x[0:3])
        else:
            actualChr = x[0]
            splittedList.append([x[0:3]])
    return splittedList

def check_boundaries_ctcf(bList, cList):
    posB = 0
    while(True):
        actualChrom = bList[pos[b]][0]
        posC = 0
        while (cList[posC][0] != actualChrom):
            posC += 1

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
    write_in_file(relationList, "relationList.bed", True)
    create_tree(relationList, mergedListWithId)
    #print(relationList)



