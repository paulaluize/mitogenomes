#This script is part of supplementary documents of "Impact of Introns and Homing Endonucleases on Structural Mitogenome Shaping in Hypocreales"
#submitted to Frontiers in Microbiology, section Fungi and Their Interactions
#Manuscript ID: 531057
#Authors:  Paula Fonseca, Fernanda Badotti, Ruth De-Paula, Daniel Araújo, Dener Eduardo Bortolini, Luiz-Eduardo Del-Bem, Vasco Ariston De Carvalho Azevedo, 
#Bertram Brenig, Eric Roberto Guimarães Rocha Aguiar, Aristóteles Góes-Neto

#This script uses the output file of Mfannot to list and save uORFs of a target species
#As described in https://github.com/BFL-lab/Mfannot, Mfannot is a program for annotation of mitochondrial and plastid genomes 

#******************************************************************************#
#                          Run the code in Python 3+                           #
#******************************************************************************#

# -*- Coding: UTF-8 -*-
#coding: utf-8

import sys
import os.path
from os import path


#Function to check if files are OK
def checkMfannotFile():
    #Check if all the necessary files names are passed as arguments
    if (len(sys.argv)!=2):
        print("\n--------------------------------------------------------------------------------------------\n")
        print ("\nUsage:\npython Mfannot2uORFs.py [file_path_name]")
        print("\n--------------------------------------------------------------------------------------------\n")
        sys.exit(0)

    mfannot_file_name=sys.argv[1]

    #Check if path/files exists
    if (path.exists(mfannot_file_name)==False):
        print("\n--------------------------------------------------------------------------------------------\n")
        print("\nOne or more files not found! Check the path and file names.\n")
        print("\n--------------------------------------------------------------------------------------------\n")
        exit(0)

    #Open input file
    input_file=open(mfannot_file_name,'r') 

    #Check if input file is a Mfannot
    if (input_file.readline().find("mfannot")==-1):
        print("\nThe file is empty or is not a mfannot output file\n")
        input_file.close()
        exit(0)

    #Get ID specie from Mfannot file name
    #The strip will remove '.\' that appear on console in Windows 10 before path\filename 
    if (os.name=="nt"):
        mfannot_file_name=mfannot_file_name.strip(".\\")
    output_file_name=mfannot_file_name[0:mfannot_file_name.find(".")]
    print(mfannot_file_name)
    print(output_file_name)

    #Open uORFs output file
    output_file=open(output_file_name+".uORFs",'w')

    return input_file, output_file

#This function look if a string contain a uORfs_name and store it in uORfs_name_vector
def findOrfs(str_name, uORfs_name_vector):
    if (str_name.find("orf")!=-1):
        uORfs_name_vector.append(str_name) 

def getuORFSNamesMfannot(uORfs_name_vector, input_file):
    #find_str_gene is a boolean variable that signalize where the block of gene names start and end
    #The Mfannot file lists all gene names in a tab format starting at line 4
    find_str_gene=0 
    #Loop to read the Mfannot input file
    for line in input_file:
        #If the string "List of genes added", then the boolean find_str_gene receives 1, signalizing that we are reading the block with gene names
        if (line.find("List of genes added")!=-1):
            find_str_gene=1 
        #If find_str_gene is true, then we read gene names
        if (find_str_gene==1):
            #Check if it is at end of gene names block
            if (line.find("end mfannot")!=-1): 
                find_str_gene=0 #para terminar de ler até o end do mfannot.
            #The gene names are structured in 3 columns of regular spaced sizes, so we read each one and stores in uORfs_name_vector
            else: 
                findOrfs(line[8:29].rstrip(' '),uORfs_name_vector)
                findOrfs(line[29:50].rstrip(' '),uORfs_name_vector)
                findOrfs(line[50:70].rstrip(' '),uORfs_name_vector) 


def getuORFsStartEndSeq(uORfs_name_vector, input_file,output_file):
    for uORfs_name in uORfs_name_vector:
        #This sets the position of reading the input_file at the start
        input_file.seek(0)
        #This boolean tell us when a sequence of a specific uORF begins
        bool_start_seq=0
        #orf_seq store the orf's sequence
        orf_seq=""
        orf_start_position=""
        orf_end_position=""
        detailed_orf_name=""
        #num_index store the index +1 after the number position in sequences lines
        num_index=-1
        #Loop to read the input_file
        for line in input_file:
            #When we find the start line with the uORFS_name, bool_start_seq receives 1 (true)
            if (line.find("-" + uORfs_name)!=-1 and line.find(" ==> start")!=-1):
                detailed_orf_name=line[1:line.find(" ==> start")].strip()
                bool_start_seq=1
            else:
                if (bool_start_seq==1):
                    #num_index store the index +1 after the number position in sequences lines 
                    num_index=line.find("  ",2)
                    #If orf_start_position is empty and bool_start_seq==1, then we get the start position of the orf
                    if (orf_start_position==""):
                        orf_start_position=line[:num_index].strip()
                    #Check if its the end of the sequence of uORfs_name, then break case true
                    if (line.find("-" + uORfs_name)!=-1 and line.find(" ==> end")!=-1):
                        print(">"+detailed_orf_name)
                        print("+"+orf_start_position)
                        print("-"+str(orf_end_position))
                        print("@"+orf_seq)
                        output_file.write(">"+detailed_orf_name+"\n")
                        output_file.write("+"+orf_start_position+"\n")
                        output_file.write("-"+str(orf_end_position)+"\n")
                        output_file.write("@"+orf_seq+"\n\n")
                        #Here we reset variables and let the loop go to the end, as is possible to have another copy
                        #forward in the file
                        orf_seq=""
                        orf_start_position=""
                        orf_end_position=""
                        detailed_orf_name=""
                        bool_start_seq=0
                    elif(line.find(";")==-1):
                        orf_seq=orf_seq + line[num_index:].strip()
                        #Calculate the orf_end_position
                        orf_end_position=int(line[:num_index].strip())+len(line[num_index:].strip())-1


def main():
    input_file,output_file=checkMfannotFile()

    #uORfs_name_vector is a array that stores the name of the uORFs listed in Mfannot file
    uORfs_name_vector=[]

    getuORFSNamesMfannot(uORfs_name_vector, input_file)

    getuORFsStartEndSeq(uORfs_name_vector, input_file, output_file)

    print("\n\n____________________________________________________________")
    print("\nResults saved in: "+output_file.name)
    print("____________________________________________________________\n\n\n")

    output_file.close() 
    input_file.close() 

    

if __name__ == '__main__':
    main()