
#This script is part of supplementary documents of "Impact of Introns and Homing Endonucleases on Structural Mitogenome Shaping in Hypocreales"
#submitted to Frontiers in Microbiology, section Fungi and Their Interactions
#Manuscript ID: 531057
#Authors:  Paula Fonseca, Fernanda Badotti, Ruth De-Paula, Daniel Araújo, Dener Eduardo Bortolini, Luiz-Eduardo Del-Bem, Vasco Ariston De Carvalho Azevedo, 
#Bertram Brenig, Eric Roberto Guimarães Rocha Aguiar, Aristóteles Góes-Neto


#This script uses the NCBI API to make a query and retrieve data from GenBank with all gene annotation of target species
#Using a text file as input, the script can retrieve data from multiples species at time
#The result is a '.cds' file, with name, ID, size of genome, start, end positions and the name of all genes

#******************************************************************************#
#                          Run the code in Python 3+                           #
#******************************************************************************#


# -*- Coding: UTF-8 -*-
#coding: utf-8
import sys
import urllib.request
import re
import os.path
from os import path

#This function reads the target(s) id(s) specie(s) from a 'txt' file
def readIDs(fileName_txt):
    ids_file=open(fileName_txt,'r')
    query_IDs = []
    for line in ids_file:
        #the '[accn]' substring is added as requirement by NCBI API
        query_IDs.append(line.rstrip()+"[accn]")
    ids_file.close()
    return query_IDs

#This function will try to retrieve the data of a target specie using NCBI API. The result is a xml string with all data.
def getXMLNCBI(str_ID):
    urlBase = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

    #create URL for esearch
    url = urlBase+"esearch.fcgi?db=nuccore&term="+str_ID+"&usehistory=y"
    f = urllib.request.urlopen(url)
    #Read xml from url
    xml = f.read()
    strXml=xml.decode("utf-8")

    #If WebEnv e QueryKey exists in this firstxml, fetch the URL query 
    objRe = re.search('<WebEnv>(\S+)<\/WebEnv>',strXml)
    webEnv = objRe.group()
    webEnv = webEnv[8:len(webEnv)-9]
    objRe = re.search('<QueryKey>(\d+)<\/QueryKey>',strXml)
    queryKey = objRe.group()
    queryKey = queryKey[10:len(queryKey)-11]
    #URL efetch
    url = urlBase+"efetch.fcgi?db=nuccore&query_key="+queryKey+"&WebEnv="+webEnv+"&rettype=gb&retmode=xml"
    f = urllib.request.urlopen(url)
    xml_esearch_bin = f.read()
    xml_esearch = xml_esearch_bin.decode("utf-8")
    return xml_esearch

#This function generates xml and cds files with the data retrieved from NCBI. The CDS file will contain all the annotated genes within the genome of a specie
def generateXMLCDS(item, xml_esearch, key_words):

    genome_ID=""
    genome_size=0
    genome_size_CDS=0
    #To treat gene overlapping, the array cds_vector will store the nucleotides positons that are part of a coding region
    #If a nucleotide is part of CDS, then its position on cds_vector will be '1', while the positions with '0' will represent the nucleotides of non coding region
    cds_vector=[]
    #The index_key_words is a counter that controls the exploration of the ordered indendation and alignment in the xml format, retrieving the gene data when found
    index_key_words=0

    #Open output files .xml and .cds
    output_xml_file=open(item[:len(item)-6]+".xml",'w')
    output_cds_file=open(item[:len(item)-6]+".cds","w")
    #List that stores start and end positions of a gene
    rangeCDS=[]
    for line in xml_esearch:
        output_xml_file.write(line+"\n")
        #The command below search for the key_word indicated by index_key_words in the present line of the xml
        id_str_found=line.find(key_words[index_key_words])
        if(id_str_found!=-1):
            #if it is the first index_key_word, get the genome ID and set to next key_word 
            if (index_key_words==0):
                genome_ID=line[id_str_found+len(key_words[index_key_words]):line.find("<",id_str_found+1)]
                index_key_words=index_key_words+1
            #If it is the second one, get genome total size and set to next key_word 
            elif (index_key_words==1):
                genome_size=int(line[id_str_found+len(key_words[index_key_words]):line.find("<",id_str_found+1)])
                #Instantiate cds_vector with size + 1 of the whole genome. The 0 position will not be used
                cds_vector = [0]*(genome_size+1)
                genome_size_CDS=0
                index_key_words=index_key_words+1
            #In the third key_word, get specie name header, write in screen and cds file and set to next key_word
            elif (index_key_words==2):
                output_cds_file.write(line[id_str_found+len(key_words[index_key_words]):line.find("<",id_str_found+1)]+"\n")
                output_cds_file.write("Genome ID: "+genome_ID+"\n")
                output_cds_file.write("Genome size: "+str(genome_size)+"\n")
                output_cds_file.write("Genes:\n")
                print(line[id_str_found+len(key_words[index_key_words]):line.find("<",id_str_found+1)])
                print("Genome ID: "+genome_ID)
                print("Genome size: "+str(genome_size))
                index_key_words=index_key_words+1
            #Here we get the information if the next data is part of a gene and go further into the xml with another key_word
            elif (index_key_words==3):
                index_key_words=index_key_words+1
            #In the fifth one, we get the string with the start and end positions of the gene, calling a function that print the values on screen and cds file, 
            #besides marking the nuclotides positions of the gene on cds_vector
            elif (index_key_words==4):
                #Treat a Join if necessary
                if (line.find("join")==-1):
                    rangeCDS.append(re.sub('[^0-9.]','',line))
                else:
                    aux_rangeCDS=re.sub('[^0-9.,]','',line)
                    rangeCDS.append(aux_rangeCDS[:aux_rangeCDS.find(",")])
                    rangeCDS.append(aux_rangeCDS[aux_rangeCDS.find(",")+1:])
                index_key_words=index_key_words+1
            #We go further into the xml indendation
            elif (index_key_words==5):
                index_key_words=index_key_words+1
            #And get the gene name, printing on screen and cds file
            elif (index_key_words==6):
                gene_name=line[line.find("<GBQualifier_value>")+19:line.find("</GBQualifier_value>",18)]
                for cds_range in rangeCDS:
                    write_start_end_gene(cds_range,output_cds_file,cds_vector)
                    output_cds_file.write("#"+gene_name+"\n")
                    print(" ("+str(gene_name)+")")
                #Then we retrocede 3 positions in index_key_value to look for another genes
                index_key_words=index_key_words-3
                #And clear rangeCDS
                rangeCDS=[]

    #Get the number of nucleotides in coding regions
    genome_size_CDS=sum(cds_vector)
    print("Sum of nucleotides in the coding regions (CDS) of the genome ID= "+genome_ID+": "+str(genome_size_CDS)+" of "+str(genome_size)+" nucleotides ("+str(round(genome_size_CDS*100/genome_size,2))+"%)")
    output_cds_file.write("Sum of nucleotides in the coding regions (CDS) of the genome: "+str(genome_size_CDS)+" of "+str(genome_size)+" nucleotides ("+str(round(genome_size_CDS*100/genome_size,2))+"%)")   
    output_cds_file.close()
    output_xml_file.close()

#This function extracts from a string the start and end position of a gene, print the data on screen and on cds file and register the position of all nucleotides on cds_vector
#by changing the '0' value to 1. The change only occurs once for a nucleotide.
def write_start_end_gene(range_gene, output_cds_file,cds_vector):
    indexRange=range_gene.find("..")
    print("\t\t"+range_gene[:indexRange]+"\t"+range_gene[indexRange+2:], end = '')
    output_cds_file.write(range_gene[:indexRange]+";"+range_gene[indexRange+2:])
    for i in range(int(range_gene[:indexRange]),int(range_gene[indexRange+2:])+1):
        if (cds_vector[i]==0):
            cds_vector[i]=cds_vector[i]+1

#Function to check if files are OK
def checkInputFiles():

    #Check if all the necessary files names are passed as arguments
    if (len(sys.argv)!=2 or sys.argv[1].find(".txt")==-1):
        print ("\nUsage:\npython getCDSGenBank.py [file_path_name.txt]")
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

    #The key words are used to read the '.xml' format return by NCBI API and extract the genes data from it
    key_words=["<GBSeq_locus>", "<GBSeq_length>", "<GBSeq_definition>", "<GBFeature_key>gene</GBFeature_key>", "<GBFeature_location>","<GBQualifier_name>","<GBQualifier_value>"]
    #This variable store ids that returns a empty result
    error_ids=""

    for item in query_IDs:
        print("\n\nQuerying ID: "+item[:len(item)-6]+"\n\n")

        #Get xml with GenBAnk data from NCBI
        xml_esearch=getXMLNCBI(item)
        if (xml_esearch.find("<ERROR>Empty result - nothing to do</ERROR>")==-1):
            xml_esearch=xml_esearch.splitlines()
            generateXMLCDS(item,xml_esearch,key_words)
        else:
            error_ids=error_ids+item[:len(item)-6]+"\n"
        print("\n--------------------------------------------------------------------------------------------\n")
        
    #If some ID returned empty, list them
    if (len(error_ids)>0):
        print("\n--------------------------------------------------------------------------------------------\n")
        print("\nThe following IDs returned a empty result:\n"+error_ids+"\nCheck these IDs and try again\n")
        print("\n--------------------------------------------------------------------------------------------\n")


if __name__ == '__main__':
    main()
        