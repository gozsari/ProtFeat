#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Authors: Chris (chris@nohup.cc) & Young (young@nohup.cc)

def usage():
    print("possum.py usage:")
    print("python possum.py <options> <source files> ")
    print("-i,--input: input a file in fasta format.")
    print("-o,--ouput: output a file of the generated feature.")
    print("-t,--type: specify a feature encoding algorithm.")
    print("-p,--pssmdir: specify the directory of pssm files.")
    print("-h,--help: show the help information.")

import fileinput
import sys, getopt
from os import listdir
from os.path import isfile, join
import re
import numpy as np
from possum_ft import *

opts, args = getopt.getopt(sys.argv[1:], 'i:o:t:p:a:b:h', ['input=','output=','type=','pssmdir=','argument=','veriable=','help'])
inputFile=""
outputFile=""
algoType=""
pssmdir=""
argument=""
veriable=""

for opt, arg in opts:
    if opt in ('-i','--input'):
        inputFile = arg
    elif opt in ('-o','--output'):
        outputFile = arg
    elif opt in ('-t','--type'):
        algoType = arg
    elif opt in ('-p','--pssmdir'):
        pssmdir = arg
    elif opt in ('-a','--argument'):
        argument = int(arg)
    elif opt in ('-b','--veriable'):
        veriable = int(arg)
    elif opt in ('-h', '--help'):
        usage()
        sys.exit(2)
    else:
        usage()
        sys.exit(2)
check_head = re.compile(r'\>')

smplist = []
smpcnt = 0
seq = ""

for line, strin in enumerate(fileinput.input(inputFile)):
    if not check_head.match(strin):
        seq += strin.strip()
        #smplist.append(strin.strip())
        #smpcnt += 1
    elif seq!="":
        smplist.append(seq)
        smpcnt += 1
        seq= ""
smplist.append(seq)
smpcnt += 1

#print "smplist="
#print smplist
onlyfiles = [ f for f in listdir(pssmdir) if isfile(join(pssmdir,f)) ]
#print "onlyfiles="
#print onlyfiles
fastaDict = {}

for fi in onlyfiles:
    cntnt = ''
    pssmContentMatrix=readToMatrix(fileinput.input(pssmdir+'/'+fi))
    #print(pssmContentMatrix)
    pssmContentMatrix=np.array(pssmContentMatrix)
    #print(pssmContentMatrix)
    sequence=pssmContentMatrix[:,0]
    seqLength=len(sequence)
    for i in range(seqLength):
        cntnt+=sequence[i]
    if cntnt in fastaDict:
        #print strin
        continue
    fastaDict[cntnt] = fi

finalist = []
for smp in smplist: 
    #print "smp="+smp
    #print "fastaDict[smp]="+fastaDict[smp]
    finalist.append(pssmdir+'/'+fastaDict[smp])

file_out = open(outputFile,'w')

for fi in finalist:
    #print "fi=" + fi
    #pass in original matrix
    input_matrix=fileinput.input(fi)
    #output a feature vector
    #input_matrix=readToMatrix(input_matrix)
    #print(input_matrix)
    feature_vector=calculateDescriptors(input_matrix, algoType, argument, veriable)
    np.savetxt(file_out, feature_vector, delimiter=",")
