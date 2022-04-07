#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Authors: Chris (chris@nohup.cc) & Young (young@nohup.cc)

import sys
import numpy as np
import math
import re
import fileinput
from matrixTransformer import *

'''
Row transformation based feature encoding algorithms
'''
def aac_pssm(input_matrix):
    #print "start aac_pssm function"
    SWITCH = 0
    COUNT = 20
    seq_cn=float(np.shape(input_matrix)[0])
    aac_pssm_matrix=handleRows(input_matrix,SWITCH,COUNT)
    aac_pssm_matrix=np.array(aac_pssm_matrix)
    aac_pssm_vector=average(aac_pssm_matrix,seq_cn)
    #print "end aac_pssm function"
    return aac_pssm_vector

def d_fpssm(input_matrix):
    #print "start d_fpssm function"
    SWITCH = 1
    COUNT = 20
    d_fpssm_matrix=handleRows(input_matrix,SWITCH,COUNT)
    seqLen = float(np.shape(input_matrix)[0])
    seq_cn = 1.0
    element_max = np.amax(d_fpssm_matrix[0], axis=0)
    element_min = np.amin(d_fpssm_matrix[0], axis=0)
    d_fpssm_vector=average(d_fpssm_matrix, seq_cn)
    d_fpssm_vector=np.add(d_fpssm_vector, -element_min)
    d_fpssm_vector=np.divide(d_fpssm_vector, element_max)
    d_fpssm_vector=np.divide(d_fpssm_vector, seqLen)
    #print "end d_fpssm function"
    return d_fpssm_vector

def smoothed_pssm(input_matrix, SMOOTH, SLIDE):
    #print "start smoothed_pssm function"
    SMOOTH = 7
    SLIDE = 50
    smoothed_pssm_matrix = window(input_matrix, SMOOTH, SLIDE)
    #print "finish smoothed_pssm function"
    return smoothed_pssm_matrix

def ab_pssm(input_matrix):
    #print "start ab_pssm function"
    seq_cn=np.shape(input_matrix)[0]
    BLOCK=int(seq_cn/20)
    #print BLOCKdp_pssm
    matrix_final=[]
    for i in range(19):
        tmp=input_matrix[i*BLOCK:(i+1)*BLOCK]
        #print tmp
        matrix_final.append(aac_pssm(tmp)[0])
    tmp=input_matrix[19*BLOCK:]
    #print tmp
    matrix_final.append(aac_pssm(tmp)[0])
    ab_pssm_matrix=average(matrix_final,1.0)
    #print "finish ab_pssm function"
    return ab_pssm_matrix

def pssm_composition(input_matrix):
    #print "start pssm_composition function"
    SWITCH = 0
    COUNT = 400
    seq_cn=float(np.shape(input_matrix)[0])
    pssm_composition_matrix=handleRows(input_matrix, SWITCH, COUNT)
    pssm_composition_vector=average(pssm_composition_matrix,seq_cn)
    #print "end pssm_composition function"
    return pssm_composition_vector

def rpm_pssm(input_matrix):
    #print "start rpm_pssm function"
    SWITCH = 1
    COUNT = 400
    seq_cn=float(np.shape(input_matrix)[0])
    rpm_pssm_matrix=handleRows(input_matrix,SWITCH,COUNT)
    rpm_pssm_vector=average(rpm_pssm_matrix,seq_cn)
    #print "end rpm_pssm function"
    return rpm_pssm_vector

def s_fpssm(input_matrix):
    #print "start s_fpssm function"
    SWITCH = 2
    COUNT = 400
    seq_cn = 1
    s_fpssm_matrix=handleRows(input_matrix, SWITCH, COUNT)
    s_fpssm_matrix = np.array(s_fpssm_matrix)
    s_fpssm_matrix_shape = np.shape(s_fpssm_matrix)
    matrix_average = [(np.reshape(s_fpssm_matrix, (s_fpssm_matrix_shape[0] * s_fpssm_matrix_shape[1], )))]
    #print "end s_fpssm function"
    return matrix_average

'''
Column transformation based feature encoding algorithms
'''
def dpc_pssm(input_matrix):
    #print "start dpc_pssm function"
    PART = 0
    STEP = 1
    ID = 0
    KEY = 0
    matrix_final = preHandleColumns(input_matrix, STEP, PART, ID)
    seq_cn = float(np.shape(input_matrix)[0])
    dpc_pssm_vector = average(matrix_final, seq_cn-STEP)
    #print "end dpc_pssm function"
    return dpc_pssm_vector

def k_separated_bigrams_pssm(input_matrix,STEP):
    #print "start k_separated_bigrams_pssm function"
    PART=1
    ID=0
    KEY=0
    STEP = 1
    matrix_final = preHandleColumns(input_matrix, STEP, PART, ID)
    seq_cn = float(np.shape(input_matrix)[0])
    k_separated_bigrams_pssm_vector=average(matrix_final,10000.0)
    #print "end k_separated_bigrams_pssm function"
    return k_separated_bigrams_pssm_vector

def tri_gram_pssm(input_matrix):
    #print "start tri_gram_pssm function"
    tri_gram_pssm_matrix = handleTriColumns(input_matrix)
    #print "end tri_gram_pssm function"
    return tri_gram_pssm_matrix

def eedp(input_matrix):
    #print "start eedp function"
    STEP = 2
    PART = 0
    ID = 1
    KEY = 0
    seq_cn = float(np.shape(input_matrix)[0])
    matrix_final = preHandleColumns(input_matrix, STEP, PART, ID)
    eedp_vector = average(matrix_final, seq_cn-STEP)
    #print "end eedp function"
    return eedp_vector

def tpc(input_matrix):
    #print "start tpc function"
    PART=0
    STEP=1
    ID=0
    KEY=1
    matrix_final=preHandleColumns(input_matrix, STEP, PART, ID)
    matrix_tmp=[0.0] * 20
    matrix_tmp=np.array(matrix_tmp)
    for i in range(20):
        matrix_tmp=list(map(sum, list(zip(matrix_final[i], matrix_tmp))))
    for i in range(20):
        for j in range(20):
            matrix_final[i][j]=matrix_final[i][j]/matrix_tmp[j]
    tpc_vector = average(matrix_final, 1.0)
    #print "end tpc function"
    return tpc_vector

'''
Mixture of row and column transformation based feature encoding algorithms
'''
def edp(input_matrix):
    #print "start edp function"
    STEP=2
    PART=0
    ID=1
    edp_matrix=[ [0.0] * 20 ] * 1
    edp_matrix=np.array(edp_matrix)
    seq_cn=float(np.shape(input_matrix)[0])
    output_matrix=preHandleColumns(input_matrix,STEP,PART,ID)
    output_matrix=np.array(output_matrix)
    for i in range(20):
        edp_matrix[0]=list(map(sum, list(zip(output_matrix[i], edp_matrix[0]))))
    edp_matrix=np.divide(edp_matrix, (seq_cn-STEP)*20.0)
    #print "end edp function"
    return edp_matrix

def rpssm(input_matrix):
    #print "start rpssm function"
    rpssm_matrix=handleMixed3(input_matrix)
    #print "finish rpssm function"
    return rpssm_matrix

def pse_pssm(input_matrix, ALPHA):
    #print "start pse_pssm function"
    ALPHA=1
    pse_pssm_matrix=handleMixed(input_matrix,ALPHA)
    #print "end pse_pssm function"
    return pse_pssm_matrix

def dp_pssm(input_matrix, ALPHA):
    #print "start dp_pssm function"
    ALPHA=5
    dp_pssm_matrix=handleMixed2(input_matrix,ALPHA)
    #print "end dp_pssm function"
    return dp_pssm_matrix

def pssm_ac(input_matrix, GROUP):
    #print "start pssm_ac function"
    ID=0
    GROUP = 10
    pssm_ac_matrix=correlation(input_matrix, ID, GROUP)
    #print "end pssm_ac function"
    return pssm_ac_matrix

def pssm_cc(input_matrix, GROUP):
    #print "start pssm_cc function"
    ID=1
    GROUP = 10
    pssm_cc_matrix=correlation(input_matrix, ID, GROUP)

    #print "finish pssm_cc function"
    return pssm_cc_matrix
'''
Feature encoding algorithms that simply combine the above single feature
'''
def aadp_pssm(input_matrix):
    aac_pssm_matrix=aac_pssm(input_matrix)
    dpc_pssm_matrix=dpc_pssm(input_matrix)
    aac_pssm_matrix=np.array(aac_pssm_matrix)
    dpc_pssm_matrix=np.array(dpc_pssm_matrix)
    aadp_pssm_matrix=np.hstack((aac_pssm_matrix, dpc_pssm_matrix))
    #print np.shape(aadp_pssm_matrix)
    return aadp_pssm_matrix

def aatp(input_matrix):
    aac_pssm_matrix=aac_pssm(input_matrix)
    tpc_matrix=tpc(input_matrix)
    aac_pssm_matrix=np.array(aac_pssm_matrix)
    tpc_matrix=np.array(tpc_matrix)
    aatp_matrix=np.hstack((aac_pssm_matrix,tpc_matrix))
    #print np.shape(aatp_matrix)
    return aatp_matrix

def medp(input_matrix):
    edp_matrix=edp(input_matrix)
    eedp_matrix=eedp(input_matrix)
    edp_matrix=np.array(edp_matrix)
    eedp_matrix=np.array(eedp_matrix)
    #print edp_matrix
    #print eedp_matrix
    medp_matrix=np.hstack((edp_matrix, eedp_matrix))
    #print medp_matrix
    #print np.shape(medp_matrix)
    return medp_matrix
