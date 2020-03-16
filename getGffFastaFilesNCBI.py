#This script is part of supplementary documents of "Impact of Introns and Homing Endonucleases on Structural Mitogenome Shaping in Hypocreales"
#submitted to Frontiers in Microbiology, section Fungi and Their Interactions
#Manuscript ID: 531057
#Authors:  Paula Fonseca, Fernanda Badotti, Ruth De-Paula, Daniel Araújo, Dener Eduardo Bortolini, Luiz-Eduardo Del-Bem, Vasco Ariston De Carvalho Azevedo, 
#Bertram Brenig, Eric Roberto Guimarães Rocha Aguiar, Aristóteles Góes-Neto

#This script download the gff and fasta files using the curl command
#Use a txt file to get data from one or more species
#As result it will create a folder with fasta and gff files for each specie in the txt file

#******************************************************************************#
#                          Run the code in Python 3+                           #
#******************************************************************************#

# -*- Coding: UTF-8 -*-
#coding: utf-8

import sys
import subprocess
import os.path
from os import path
from pathlib import Path

#This function reads the target(s) id(s) specie(s) from a 'txt' file
def readIDs(fileName_txt):
    ids_file=open(fileName_txt,'r')
    query_IDs = []
    for line in ids_file:
        query_IDs.append(line.rstrip())
    ids_file.close()
    return query_IDs

def checkInputFiles():

    #Check if all the necessary files names are passed as arguments
    if (len(sys.argv)!=2 or sys.argv[1].find(".txt")==-1):
        print("\n--------------------------------------------------------------------------------------------\n")
        print ("\nUsage:\npython getGffFastaFilesNCBI.py [file_path_name.txt]")
        print("\n--------------------------------------------------------------------------------------------\n")
        sys.exit(0)

    fileName_txt=sys.argv[1]

    #Check if path/files exists
    if (path.exists(fileName_txt)==False):
        print("\nTXT file not found! Check the path and file name.\n")
        exit(0)

    return fileName_txt

def main():
    
    fileName_txt=checkInputFiles()

    query_IDs=readIDs(fileName_txt)

    #This variable store ids that returns a empty result
    error_ids=""

    for item in query_IDs:
        print("\n\nQuerying ID: "+item+"\n\n")

        #Retrieve gff file. The --insecure options is to avoid revocation function error 
        gff_document=subprocess.check_output(["curl", "--insecure", "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id="+item])
        gff_document_dec=gff_document.decode("utf-8")

        #Retrieve fasta file. The --insecure options is to avoid revocation function error 
        fasta_document=subprocess.check_output(["curl", "--insecure", "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=fasta&id="+item])
        fasta_document_dec=fasta_document.decode("utf-8")
        
        if (fasta_document_dec.find("Failed to understand id")==-1 and gff_document_dec.find("Failed to understand id")==-1):
            #check if folder exists
            folder = Path(item)
            if not folder.exists():
                #create the folder
                os.mkdir(item)
            #Open and write gff file
            gff_output_file=open(item+"/"+item+".gff",'w')
            gff_output_file.write(gff_document_dec)
            gff_output_file.close()
            #Open and write fasta file
            fasta_output_file=open(item+"/"+item+".fasta",'w')
            fasta_output_file.write(fasta_document_dec)
            fasta_output_file.close()
        else:
            error_ids=error_ids+item+"\n"
        print("\n--------------------------------------------------------------------------------------------\n")
        
    #If some ID returned empty, list them
    if (len(error_ids)>0):
        print("\n--------------------------------------------------------------------------------------------\n")
        print("\nThe following IDs returned a empty result:\n"+error_ids+"\nCheck these IDs and try again\n")
        print("\n--------------------------------------------------------------------------------------------\n")

if __name__ == '__main__':
    main()


    

        