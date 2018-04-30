#--------------biblio & modules-------------------------------------
import sys
from sys import argv
from Bio import AlignIO 
from Bio import SeqIO
from Bio.Seq import Seq
#--------------extracting file-----------------------------

#-----------modules-----------------------------------------

round_ = 4  # The round parameter! 


def nuc_identification_module(line): #looking for nuc and nuc% in line
    array = list(line)
    result = {i: array.count(i) for i in array};
    return result


def duplet_compute_module(duplet_list, A, B):
    A_keys = []
    B_keys = []
    for key in A.keys():
        A_keys = A_keys + [key]
    for key in B.keys():
        B_keys = B_keys + [key]

    p1 = round(A[A_keys[0]]/(A[A_keys[0]] + A[A_keys[1]]), round_)
    p2 = round(A[A_keys[1]]/(A[A_keys[0]] + A[A_keys[1]]), round_)
    q1 = round(B[B_keys[0]]/(B[B_keys[0]] + B[B_keys[1]]), round_)
    q2 = round(B[B_keys[1]]/(B[B_keys[0]] + B[B_keys[1]]), round_)

    p1q1 = str(A_keys[0]) + str(B_keys[0])
    p1q2 = str(A_keys[0]) + str(B_keys[1])
    p2q1 = str(A_keys[1]) + str(B_keys[0])
    p2q2 = str(A_keys[1]) + str(B_keys[1])

    counter_p1q1_pr = round(float( duplet_list.count(p1q1) / len(duplet_list)), round_)
    counter_p1q2_pr = round(float( duplet_list.count(p1q2) / len(duplet_list)), round_)
    counter_p2q1_pr = round(float( duplet_list.count(p2q1) / len(duplet_list)), round_)
    counter_p2q2_pr = round(float( duplet_list.count(p2q2) / len(duplet_list)), round_)
    return counter_p1q1_pr, counter_p1q2_pr, counter_p2q1_pr, counter_p2q2_pr, p1, p2, q1, q2



#------------------------------------------------------------------------
 
out_2_slices_file = open(str(argv[1]), "r") 

LD_out_file = open("LD_" + str(argv[2]) + ".txt", "w")
LD_out_file.write("Len" + "\t" + "LD" + "\n")

line = out_2_slices_file.readline().strip().split()
calc = 1               


while len(line) != 0:
    if line[0] == ">": 
        l = line[4]
        l = int(l)
        
        line1 = out_2_slices_file.readline().strip()
        line2 = out_2_slices_file.readline().strip()
        A = nuc_identification_module(line1)
        B = nuc_identification_module(line2)

        duplet_list = []
        for i in range(len(line1)):
            duplet_list = duplet_list + [line1[i] + line2[i]]
        AB = duplet_compute_module(duplet_list, A, B) 
 
        LD = AB[0]*AB[3] - AB[1]*AB[2]
        if LD >= 0:
            LD = round(LD/(min(AB[4]*AB[7], AB[5]*AB[6])), round_)
        if LD < 0:
            LD = round(LD/(max( -AB[4]*AB[6], -AB[5]*AB[7])), round_)   
  
    LD_out_file.write(str(l) + "\t" + str(LD) + "\n")

    line = out_2_slices_file.readline().strip().split()

out_2_slices_file.close()    

