#This script is part of supplementary documents of "Impact of Introns and Homing Endonucleases on Structural Mitogenome Shaping in Hypocreales"
#submitted to Frontiers in Microbiology, section Fungi and Their Interactions
#Manuscript ID: 531057
#Authors:  Paula Fonseca, Fernanda Badotti, Ruth De-Paula, Daniel Araújo, Dener Eduardo Bortolini, Luiz-Eduardo Del-Bem, Vasco Ariston De Carvalho Azevedo, 
#Bertram Brenig, Eric Roberto Guimarães Rocha Aguiar, Aristóteles Góes-Neto

#This script uses a gff and fasta files of a target species to get the sequences of genes of interest (GOI)
#The gff file provide the start and end positions of each GOI.
#The files used are as follow: 
# -gff
# -fasta
# Both can be downloaded from NCBI using getGffFastFilesNCBI.py script
#The output file is in fasta format, with the names (preceded by '>') and the sequence of  the genes of interest (GOI)

#******************************************************************************#
#                          Run the code in Python 3+                           #
#******************************************************************************#



import sys
import os.path
from os import path


def checkInputFiles():
    #Check if all the necessary files names are passed as arguments
    if (len(sys.argv)!=3 or sys.argv[1].find(".gff")==-1 or sys.argv[2].find(".fasta")==-1):
        print("\n--------------------------------------------------------------------------------------------\n")
        print ("\nUsage:\npython getGeneSeqGff.py [file_path_name.gff] [file_path_name.fasta]\n\n")
        print("\n--------------------------------------------------------------------------------------------\n")
        sys.exit(0)

    gff_file_name=sys.argv[1]
    fasta_file_name=sys.argv[2]


    #Check if path/files exists
    if (path.exists(gff_file_name)==False or path.exists(fasta_file_name)==False):
        print("\n--------------------------------------------------------------------------------------------\n")
        print("\nOne or more files not found! Check the path and file names.\n")
        print("\n--------------------------------------------------------------------------------------------\n")
        exit(0)

    return gff_file_name, fasta_file_name

#Reads whole sequence from input fasta file
def readFasta(fasta_file):
    #Position 0 of whole_genome will not be used
    whole_genome=" "
    for line in fasta_file:
        if (line.find(">")==-1):
            whole_genome=whole_genome+line.strip()
    fasta_file.close()
    return whole_genome

#This function check if the gene_name is present in the whitelist array, case true, it saves the name and sequence of the gene 
def search_genes(gene_name,gene_sequence, output_file, GOI):
    for item in GOI:
        if (gene_name.startswith(item)):
            output_file.write(">"+gene_name+"\n"+gene_sequence+"\n\n")
            break

#Read the gff file to extract data.
#The gff file contains 1 gene per row with several values ordered by 'tab'. Its straight forward to get the name and positions of a single gene
#and retrieve th sequence from the whole_genome
def readGffSelGenes(GOI, gff_file, output_file, whole_genome):
    for line in gff_file:
        values=line.split("\t")
        if (len(values)==9):
            #Start position is at index 3
            start_gene_position=values[3]
            #End position is at index 4
            end_gene_position=values[4]
            #Name is at index 8
            gene_name=values[8][values[8].find("Name=")+5:].strip()
            #Get gene sequence
            gene_sequence=whole_genome[int(start_gene_position):int(end_gene_position)+1]
            #Verify if it is on gene_whitelist to save in output file
            search_genes(gene_name,gene_sequence,output_file,GOI)

def main():

    GOI={"rrnL","rps3","nad2","nad3","atp9","cox2","nad4l","nad5","cob","cox1","nad1","nad4","atp8","atp6","rrnS","cox3","nad6"}

    gff_file_name, fasta_file_name=checkInputFiles()

    gff_file=open(gff_file_name,'r')
    fasta_file=open(fasta_file_name,'r')

    whole_genome=readFasta(fasta_file)

    #Get ID specie from gff file name
    #The strip will remove '.\' that appear on console in Windows 10 before path\filename 
    if (os.name=="nt"):
        gff_file_name=gff_file_name.strip(".\\")

    output_file_name=gff_file_name[0:gff_file_name.find(".")]+"_GOI.fasta"
    #Open output file with '_GOI.fasta' extension
    output_file=open(output_file_name,'w')

    readGffSelGenes(GOI,gff_file,output_file, whole_genome)

    gff_file.close()
    output_file.close()
    print("\n\n____________________________________________________________")
    print("\nResults saved in: "+output_file_name)
    print("____________________________________________________________\n\n\n")

if __name__ == '__main__':
    main()