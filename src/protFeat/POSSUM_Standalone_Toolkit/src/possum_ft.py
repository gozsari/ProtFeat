#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Authors: Chris (chris@nohup.cc) & Young (young@nohup.cc)

import sys
import numpy as np
import math
import re
import fileinput
from featureGenerator import *

def readToMatrix(input_matrix):
    #print "start to read PSSM matrix"
    PSSM = []
    p = re.compile(r'-*[0-9]+')
    for line, strin in enumerate(input_matrix):
        if line > 2:
            str_vec = []
            overall_vec = strin.split()
            #print len(overall_vec)
            if len(overall_vec) == 0:
                break
            str_vec.extend(overall_vec[1])
            if(len(overall_vec) < 44):
                print("There is a mistake in the pssm file")
                print("Try to correct it")
                for cur_str in overall_vec[2:]:
                    str_vec.extend(p.findall(cur_str))
                    if(len(str_vec) >= 21):
                        if(len(str_vec)) >21:
                            #print len(str_vec)
                            #print str_vec
                            #print overall_vec
                            #print "Exit with an error"
                            exit(1)
                        break;
                print("Done")
            else:
                str_vec = strin.split()[1:42]
            if len(str_vec) == 0:
                break
            #str_vec_positive=map(int, str_vec[1:])
            PSSM.append(str_vec)
    fileinput.close()
    #print "finish to read PSSM matrix"
    PSSM = np.array(PSSM)
    return PSSM

def calculateDescriptors(input_matrix, algoType, argument, veriable):
    input_matrix=readToMatrix(input_matrix)

    if algoType=="aac_pssm":
        return aac_pssm(input_matrix)
    elif algoType=="d_fpssm":
        return d_fpssm(input_matrix)
    elif algoType=="smoothed_pssm":
        return smoothed_pssm(input_matrix, argument, veriable)
    elif algoType=="ab_pssm":
        return ab_pssm(input_matrix)
    elif algoType=="pssm_composition":
        return pssm_composition(input_matrix)
    elif algoType=="rpm_pssm":
        return rpm_pssm(input_matrix)
    elif algoType=="s_fpssm":
        return s_fpssm(input_matrix)
    elif algoType=="dpc_pssm":
        return dpc_pssm(input_matrix)
    elif algoType=="k_separated_bigrams_pssm":
        return k_separated_bigrams_pssm(input_matrix, argument)
    elif algoType=="tri_gram_pssm":
        return tri_gram_pssm(input_matrix)
    elif algoType=="eedp":
        return eedp(input_matrix)
    elif algoType=="tpc":
        return tpc(input_matrix)
    elif algoType=="edp":
        return edp(input_matrix)
    elif algoType=="rpssm":
        return rpssm(input_matrix)
    elif algoType=="pse_pssm":
        return pse_pssm(input_matrix, argument)
    elif algoType=="dp_pssm":
        return dp_pssm(input_matrix, argument)
    elif algoType=="pssm_ac":
        return pssm_ac(input_matrix, argument)
    elif algoType=="pssm_cc":
        return pssm_cc(input_matrix, argument)
    elif algoType=="aadp_pssm":
        return aadp_pssm(input_matrix)
    elif algoType=="aatp":
        return aatp(input_matrix)
    elif algoType=="medp":
        return medp(input_matrix)
    
