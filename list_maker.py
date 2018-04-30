#--------------biblio & modules-------------------------------------
import sys
from sys import argv
from Bio import AlignIO 
from Bio import SeqIO
from Bio.Seq import Seq
#--------------extracting file-----------------------------

if len(sys.argv) == 1:
    print("Invalid arguments")
    sys.exit (1)

in_file = open(argv[1], "r")  

directory = str(argv[2]) + "_" + str(argv[3])  #Creating a directory

out_slice_file = open("out_slice_file_" + directory + ".txt", "w")
out_slice_file.write("№" + "\t" + "Slice" + "\n")

list_file = open("list_file.txt", "w")


#-----------modules-----------------------------------------

round_ = 4  # The round parameter! 
h = float(argv[3]) #The mask parameter!


def slice_module(slice_number, align_slice):   #It takes a slice bar os sequences and write in to the out_file
    align_slice = align_slice.upper()
    A = align_slice.count("A")
    T = align_slice.count("T")
    G = align_slice.count("G")
    C = align_slice.count("C")
    counter_A_pr = round(float( align_slice.count("A") / len(align_slice)), round_)
    counter_T_pr = round(float( align_slice.count("T") / len(align_slice)), round_)
    counter_G_pr = round(float( align_slice.count("G") / len(align_slice)), round_)
    counter_C_pr = round(float( align_slice.count("C") / len(align_slice)), round_)
    counter_N_pr = round(float( align_slice.count("N") / len(align_slice)), round_)
    counter__pr = round(float( align_slice.count("-") / len(align_slice)), round_)
    if ((counter_A_pr or counter_T_pr or counter_G_pr or counter_C_pr or counter_N_pr or counter__pr) == 1) or (counter_N_pr or counter__pr != 0):
        return 0  
    else:
        if A+T == len(align_slice) or A+G == len(align_slice) or A+C == len(align_slice) or T+G == len(align_slice) or T+C == len(align_slice) or G+C == len(align_slice):
            if min(counter_A_pr, counter_T_pr) > h or min(counter_A_pr, counter_G_pr) > h or min(counter_A_pr, counter_C_pr) > h or min(counter_T_pr, counter_G_pr) > h or min(counter_T_pr, counter_C_pr) > h or min(counter_G_pr, counter_C_pr) > h:
                out_slice_file.write(str(slice_number) + "\t" + str(align_slice) + "\n")
                list_file.write(str(slice_number) + "\n")
        return 1

#--------------import alignment----------------------------------

alignment = AlignIO.read(in_file, "fasta") #Importing all alignments from in file
in_file.close()


align_slice = ''

genome_size = len(alignment[0, :])
for slice_number in range(genome_size):  # Taking a nucleotide from each bar/slice in alignment and from each sequence. For nuc in range of length of seq
    for seq in range(len(alignment)):  # for seq from all seqs
        align_slice = align_slice + str(alignment[seq, slice_number]) #Making slices
    a = slice_module(slice_number, align_slice) #Collecting it and putting in the slice_module
    align_slice = ''

out_slice_file.close()    

#----------------------------------------------------------------

slices_file = open("out_slice_file_" + directory + ".txt", "r")
line_id = slices_file.readline() # Take out the heat of the table
slice_id_dic = {}
t_calc = 0
																																																																																																																																																																													
while len(line_id) > 0:
    line_id = slices_file.readline().strip().split()  # Reading a slice in the file
    if len(line_id) != 0:
        slice_id_dic[line_id[0]] = line_id[1]  # Making a dict with seqs and keys
        t_calc += 1

out_slice_file.close()
snps = (len(slice_id_dic))
print("			- Done.")
																																																																											

#---------------------------------------------------------------------------------


c = 0
for key in list(slice_id_dic):  # For keys in dict
    one = slice_id_dic[key]  #retrieving a slice 
    out_2_slices_file = open("out_2_slices_file_" + directory + "_" + key + ".txt", "w")
    for q in list(slice_id_dic):
        if key != q:
            summ = min(abs(int(q) - int(key)), genome_size - abs(int(key) - int(q))) #Calculating the distance
            out_2_slices_file.write(">" + "\t" + str(key) + "\t" + "-" + "\t" + str(q) + "\t" + str(summ) + "\n") #Writing into out_2_slices_file.txt - the temp file - all non snps slices
            out_2_slices_file.write(one + "\n")
            out_2_slices_file.write(slice_id_dic[q] + "\n")
    c += 1 
    out_2_slices_file.close() 

#------------------------------------------------------------------------

