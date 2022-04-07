#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Authors: Chris (chris@nohup.cc) & Young (young@nohup.cc)

import sys
import numpy as np
import math
import re
import fileinput


def average(matrixSum, seqLen):
    # average the summary of rows
    matrix_array = np.array(matrixSum)
    matrix_array = np.divide(matrix_array, seqLen)
    matrix_array_shp = np.shape(matrix_array)
    matrix_average = [(np.reshape(matrix_array, (matrix_array_shp[0] * matrix_array_shp[1], )))]
    return matrix_average

def normalizePSSM(PSSM):
    PSSM=PSSM[:,1:21]
    PSSM=PSSM.astype(float)
    PSSM=np.array(PSSM)
    seq_cn=np.shape(PSSM)[0]
    PSSM_norm=[ [0.0] * 20 ] * seq_cn
    PSSM_norm=np.array(PSSM_norm)
    mean_matrix=np.mean(PSSM, axis=1)
    std_matrix=np.std(PSSM, axis=1)

    for i in range(seq_cn):
        for j in range(20):
            if std_matrix[i]==0.0:
                PSSM_norm[i][j] = PSSM[i][j]-mean_matrix[i]
            else:
                PSSM_norm[i][j]=(PSSM[i][j]-mean_matrix[i])/std_matrix[i]
    return PSSM_norm

def window(PSSM, w_smth, w_slide):
    # 0-19 represents amino acid 'ARNDCQEGHILKMFPSTWYV'
    w_smth=int(w_smth)
    w_slide=int(w_slide)
    Amino_vec = "ARNDCQEGHILKMFPSTWYV"

    PSSM=PSSM[:,1:21]
    PSSM=PSSM.astype(float)

    seq_cn = np.shape(PSSM)[0]

    #original PSSM
    PSSM_smth = np.array([[0.0]*20]*seq_cn)
    PSSM_orig = np.array(PSSM)
    #print PSSM_orig
    #section for PSSM_smth features
    PSSM_smth_full=pssm_smth(PSSM_orig,PSSM_smth,w_smth,seq_cn)
    #print PSSM_smth_full
    PSSM_smth_final=[[0.0]*20]*w_slide
    PSSM_smth_final=np.array(PSSM_smth_final)
    for i in range(w_slide):
        PSSM_smth_final[i]=PSSM_smth_full[i]
    matrix_final=average(PSSM_smth_final,1.0)
    return matrix_final

def pssm_smth(PSSM_orig, PSSM_smth, w_smth, l):
    for i in range(l):
        if i <(w_smth-1)/2:
            for j in range(i+int((w_smth-1)/2)+1):
                PSSM_smth[i]+=PSSM_orig[j]
        elif i>=(l-(w_smth-1)/2):
            for j in range(i-int((w_smth-1)/2),l):
                PSSM_smth[i]+=PSSM_orig[j]
        else:
            for j in range(i-int((w_smth-1)/2),i+int((w_smth-1)/2)+1):
                PSSM_smth[i]+=PSSM_orig[j]
    return PSSM_smth

def handleRows(PSSM, SWITCH, COUNT):
    '''
    if SWITCH=0, we filter no element.
    if SWITCH=1, we filter all the negative elements.
    if SWITCH=2, we filter all the negative and positive elements greater than expected.
    '''
    '''
    if COUNT=20, we generate a 20-dimension vector.
    if COUNT=400, we generate a 400-dimension vector.
    '''
    # 0-19 represents amino acid 'ARNDCQEGHILKMFPSTWYV'
    Amino_vec = "ARNDCQEGHILKMFPSTWYV"
    c = COUNT/20
    matrix_final = [ [0.0] * 20 ] * int(c)
    matrix_final=np.array(matrix_final)
 
    seq_cn = 0

    PSSM_shape=np.shape(PSSM)
    for i in range(PSSM_shape[0]):
        seq_cn += 1
        str_vec=PSSM[i]
        str_vec_positive=list(map(int, str_vec[1:21]))
        str_vec_positive=np.array(str_vec_positive)
        if SWITCH==1:
            str_vec_positive[str_vec_positive<0]=0
        elif SWITCH==2:
            str_vec_positive[str_vec_positive<0]=0
            str_vec_positive[str_vec_positive>7]=0
        #print "str_vec_positive="
        #print str_vec_positive
        if COUNT==20:
            
            matrix_final[0]=list(map(sum, list(zip(str_vec_positive, matrix_final[0]))))
        elif COUNT==400:
            
            matrix_final[Amino_vec.index(str_vec[0])] = list(map(sum, list(zip(str_vec_positive, matrix_final[Amino_vec.index(str_vec[0])]))))

        #print "matrix_final="
        #print matrix_final

    return matrix_final
    
def preHandleColumns(PSSM,STEP,PART,ID):
    '''
    if STEP=k, we calculate the relation betweem one residue and the kth residue afterward.
    '''
    '''
    if PART=0, we calculate the left part of PSSM.
    if PART=1, we calculate the right part of PSSM.
    '''
    '''
    if ID=0, we product the residue-pair.
    if ID=1, we minus the residue-pair.
    '''
    '''
    if KEY=1, we divide each element by the sum of elements in its column.
    if KEY=0, we don't perform the above process.
    '''
    if PART==0:
        #print "PART=",PART
        PSSM=PSSM[:,1:21]
    elif PART==1:
        #print "PART=",PART
        PSSM=PSSM[:, 21:]
    PSSM=PSSM.astype(float)
    matrix_final = [ [0.0] * 20 ] * 20
    matrix_final=np.array(matrix_final)
    seq_cn=np.shape(PSSM)[0]
    #print "seq_cn=",seq_cn

    if ID==0:
        #print "ID=",ID
        for i in range(20):
            for j in range(20):
                for k in range(seq_cn - STEP):
                    matrix_final[i][j]+=(PSSM[k][i]*PSSM[k+STEP][j])

    elif ID==1:
        #print "ID=",ID
        for i in range(20):
            for j in range(20):
                for k in range(seq_cn - STEP):
                    matrix_final[i][j] += ((PSSM[k][i]-PSSM[k+STEP][j]) * (PSSM[k][i]-PSSM[k+STEP][j])/4.0)
    #print matrix_final
    return matrix_final

def handleTriColumns(PSSM):
    matrix_final=[ [ [0.0] * 20] * 20] * 20
    matrix_final=np.array(matrix_final)
    PSSM=PSSM[:, 21:]
    PSSM=PSSM.astype(float)
    PSSM=np.array(PSSM)
    #print "PSSM="
    #print PSSM
    seq_cn=np.shape(PSSM)[0]
    for m in range(20):
        for n in range(20):
            for r in range(20):
                for i in range(seq_cn-2):
                    matrix_final[m][n][r]+=(PSSM[i][m]*PSSM[i+1][n]*PSSM[i+2][r])
    #print "matrix_final="
    #print matrix_final
    matrix_final=np.divide(matrix_final,1000000.0)
    matrix_final_shape=np.shape(matrix_final)
    #print matrix_final
    matrix_result = [(np.reshape(matrix_final, (matrix_final_shape[0] * matrix_final_shape[1] * matrix_final_shape[2], )))]
    return matrix_result

def handleMixed(PSSM,ALPHA):
    row1=[0.0] * 20
    row2=[0.0] * 20

    matrix_final = [ [0.0] * 40 ] * 1
    row1=np.array(row1)
    row2=np.array(row2)
    matrix_final=np.array(matrix_final)

    PSSM_norm=normalizePSSM(PSSM)
    seq_cn=np.shape(PSSM)[0]
    for i in range(seq_cn):
        #print PSSM_norm[i]
        row1=list(map(sum, list(zip(row1, PSSM_norm[i]))))
    #print row1
    row1=np.divide(row1,seq_cn)

    for j in range(20):
        for i in range(seq_cn-ALPHA):
            row2[j]+=(PSSM_norm[i][j]-PSSM_norm[i+ALPHA][j])*(PSSM_norm[i][j]-PSSM_norm[i+ALPHA][j])
    #print row2
    row2=np.divide(row2,seq_cn-ALPHA)

    row=np.hstack((row1,row2))
    matrix_final[0]=row
    #print np.shape(matrix_final)
    return matrix_final

def handleMixed2(PSSM,ALPHA):
    row1=[0.0] * 40
    row2=[ [0.0] * (2*ALPHA) ]*20
    matrix_final = [ [0.0] * (40+40*ALPHA) ] * 1

    row1=np.array(row1)
    row2=np.array(row2)
    matrix_final=np.array(matrix_final)

    PSSM_norm=normalizePSSM(PSSM)
    seq_cn=np.shape(PSSM)[0]
    for j in range(20):
        positive_count_1=0
        negative_count_1=0
        for i in range(seq_cn):
            if PSSM_norm[i][j]>=0:

                positive_count_1+=1
                row1[2*j]+=PSSM_norm[i][j]
            elif PSSM_norm[i][j]<0:

                negative_count_1+=1
                row1[2*j+1]+=PSSM_norm[i][j]
        #print "positive_count_1="
        #print positive_count_1
        #print "negative_count_1="
        #print negative_count_1
        row1[2*j]=row1[2*j]/positive_count_1
        row1[2*j+1]=row1[2*j+1]/negative_count_1

    for j in range(20):
        for alpha in range(1,ALPHA+1):
            positive_count_2=0
            negative_count_2=0
            for i in range(seq_cn-alpha):
                if (PSSM_norm[i][j]-PSSM_norm[i+alpha][j])>=0:
                    positive_count_2+=1
                    row2[j][2*alpha-2]+=(PSSM_norm[i][j]-PSSM_norm[i+alpha][j])*(PSSM_norm[i][j]-PSSM_norm[i+alpha][j])
                elif (PSSM_norm[i][j]-PSSM_norm[i+alpha][j])<0:
                    negative_count_2+=1
                    row2[j][2*alpha-1]+=(PSSM_norm[i][j]-PSSM_norm[i+alpha][j])*(PSSM_norm[i][j]-PSSM_norm[i+alpha][j])
            row2[j][2*alpha-2]=row2[j][2*alpha-2]/positive_count_2
            row2[j][2*alpha-1]=row2[j][2*alpha-1]/negative_count_2

    row2=average(row2,1.0)
    row=np.hstack((row1,row2[0]))
    matrix_final[0]=row
    #print np.shape(matrix_final)
    return matrix_final

def handleMixed3(PSSM):
    row=[ [0.0] * 110 ]*1
    row1=[0.0] * 100
    row2=[0.0] * 10
    row=np.array(row)
    row1=np.array(row1)
    row2=np.array(row2, dtype='float')

    seq_cn=np.shape(PSSM)[0]
    Amino_vec = "ARNDCQEGHILKMFPSTWYV"
    RPSSM=[ [0.0]*10 ]*seq_cn
    RPSSM=np.array(RPSSM)

    PSSM=PSSM[:,1:21]
    PSSM=PSSM.astype(float)
    PSSM=np.array(PSSM)
    #print "PSSM="
    #print PSSM
    RPSSM[:,0]=np.divide(list(map(sum,list(zip(PSSM[:,13],PSSM[:,17],PSSM[:,18])))), 3.0)
    RPSSM[:,1]=np.divide(list(map(sum,list(zip(PSSM[:,10],PSSM[:,12])))), 2.0)
    RPSSM[:,2]=np.divide(list(map(sum,list(zip(PSSM[:,9],PSSM[:,19])))), 2.0)
    RPSSM[:,3]=np.divide(list(map(sum,list(zip(PSSM[:,0],PSSM[:,15],PSSM[:,16])))), 3.0)
    RPSSM[:,4]=np.divide(list(map(sum,list(zip(PSSM[:,2],PSSM[:,8])))), 2.0)
    RPSSM[:,5]=np.divide(list(map(sum,list(zip(PSSM[:,5],PSSM[:,6],PSSM[:,3])))), 3.0)
    RPSSM[:,6]=np.divide(list(map(sum,list(zip(PSSM[:,1],PSSM[:,11])))), 2.0)
    RPSSM[:,7]=PSSM[:,4]
    RPSSM[:,8]=PSSM[:,7]
    RPSSM[:,9]=PSSM[:,14]
    # print "RPSSM="
    # print RPSSM
    mean_matrix=np.mean(RPSSM, axis=0)
    # print "mean_matrix="
    # print mean_matrix
    for j in range(10):
        for i in range(seq_cn):
            row2[j]+=(RPSSM[i][j]-mean_matrix[j])*(RPSSM[i][j]-mean_matrix[j])
    #print "row2="
    #print row2
    row2=np.divide(row2,seq_cn)
    matrix_final = [ [0.0] * 10 ] * 10
    matrix_final=np.array(matrix_final)
    for i in range(10):
        for j in range(10):
            for k in range(seq_cn-1):
                matrix_final[i][j] += ((RPSSM[k][i]-RPSSM[k+1][j]) * (RPSSM[k][i]-RPSSM[k+1][j])/2.0)
    #print "matrix_final="
    #print matrix_final
    row1=average(matrix_final,seq_cn-1)[0]
    #print "row1="
    #print row1
    row[0]=np.hstack((row1,row2))
    return row

def correlation(PSSM,ID,GROUP):
    # 0-19 represents amino acid 'ARNDCQEGHILKMFPSTWYV'
    Amino_vec = "ARNDCQEGHILKMFPSTWYV"
    #GROUP=10
    PSSM=PSSM[:,1:21]
    PSSM=PSSM.astype(float)
    #print PSSM
    #section for PSSM_AC features
    seq_cn = np.shape(PSSM)[0]
    g=GROUP
    l=seq_cn
    if ID==0:
        matrix_final=pssm_ac_cal(PSSM,g,l)
    elif ID==1:
        matrix_final=pssm_cc_cal(PSSM,g,l)

    matrix_final =  average(matrix_final, l)
    return matrix_final

def pssm_ac_cal(PSSM, g, l):

    PSSM_AC = np.array([ [0.0] * 20 ] * g)

    for pg in range(g):
        l_g = l - pg - 1
        for pj in range(20):
            sum_jl = 0.0
            for i in range(l):
                sum_jl += PSSM[i][pj]
            sum_jl /= l

            pssm_acjg = 0.0
            for i in range(l_g):
                pssm_acjg +=  (PSSM[i][pj]-sum_jl) * (PSSM[i+pg+1][pj]-sum_jl)
            pssm_acjg /= l_g
            PSSM_AC[pg][pj] = pssm_acjg
    return PSSM_AC

def pssm_cc_cal(PSSM, g, l):
    PSSM_CC = np.array([ [0.0] * 380 ] * g)
    for pg in range(g):
        l_g = l - pg - 1
        for pj_1 in range(20):
            sum_jl_1 = 0.0
            for i in range(l):
                sum_jl_1 += PSSM[i][pj_1]
            sum_jl_1 /= l

            for pj_2 in range(20):
                if pj_2!=pj_1:
                    sum_jl_2 = 0.0
                    for i in range(l):
                        sum_jl_2 += PSSM[i][pj_2]
                    sum_jl_2 /= l
                    pssm_acjg = 0.0
                    for i in range(l_g):
                        pssm_acjg += (PSSM[i][pj_1]-sum_jl_1) * (PSSM[i+pg+1][pj_2]-sum_jl_2)
                    pssm_acjg /= l_g
                    if(pj_1<pj_2):
                        PSSM_CC[pg][19*pj_1+(pj_2-1)] = pssm_acjg
                    else:
                        PSSM_CC[pg][19*pj_1+(pj_2)] = pssm_acjg

    return PSSM_CC
