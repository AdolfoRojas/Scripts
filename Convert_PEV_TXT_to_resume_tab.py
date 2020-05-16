#!/usr/bin/env python3
##################################################################
import sys
sys.path.append("/usr/local/lib/python2.7/dist-packages")
import re 
##################################################################

#file1 = open(input("PEV-TXT File : "), 'r')

workfile = "5cN8XnXzquyzuf4B.txt"
out_file = "Resume_PEV_" + str(workfile) + ".tsv"
file1 = open(workfile, 'r')
file2 = open(out_file, 'w+')
file2.write("Upload_variant\tLocation\tAllele\tConsequence\tIMPACT\tSYMBOL\tFeature_type\tNFeatures(Ntranscripts)\tFeatures\tBIOTYPE\tExisting_variation\n")
Lines = file1.readlines()[1:]
print(str(len(Lines)) + " Efectos en el Archivo")
Variants = list()
for line in Lines:
    l1 = line.strip().split("\t")
    Variant = l1[0]
    Variants.append(Variant)
print(len(Variants))

Variants = list(dict.fromkeys(Variants))
print(str(len(Variants)) + " Variantes involucradas")



for element in Variants:
    Location = list()
    Allele = list()
    Consequence = list()
    IMPACT = list()
    SYMBOL = list()
    Feature_type = list()
    Feature = list()
    BIOTYPE = list()
    Existing_variation = list()
    for line in Lines:               
        if re.match(element, line):
            l2 = line.strip().split("\t")
            Location.append(l2[1])
            Allele.append(l2[2])
            Consequence.append(l2[3])
            IMPACT.append(l2[4])
            SYMBOL.append(l2[5])
            Feature_type.append(l2[7])
            Feature.append(l2[8])
            BIOTYPE.append(l2[9])
            Existing_variation.append(l2[19])
    Location = list(dict.fromkeys(Location))
    Allele = list(dict.fromkeys(Allele))
    Consequence = list(dict.fromkeys(Consequence))
    IMPACT = list(dict.fromkeys(IMPACT))
    SYMBOL = list(dict.fromkeys(SYMBOL))
    count_transcripts = Feature_type.count("Transcript")
    Feature_type = list(dict.fromkeys(Feature_type))    
    Feature = list(dict.fromkeys(Feature))
    BIOTYPE = list(dict.fromkeys(BIOTYPE))
    Existing_variation = list(dict.fromkeys(Existing_variation))
    out_line = (str(element) + "\t" + str(Location) + "\t" + str(Allele) + "\t"+ str(Consequence) + "\t"+ str(IMPACT) + "\t"+ str(SYMBOL) + "\t"+ str(Feature_type) + "\t" + str(len(Feature)) +"(" + str(count_transcripts)+")" + "\t" + str(Feature) + "\t" + str(BIOTYPE) + "\t" + str(Existing_variation) + "\n")
    out_line= out_line.replace("'", "")
    out_line= out_line.replace("]", "")
    out_line= out_line.replace("[", "")
    out_line= out_line.replace(",", ";")
    out_line= out_line.replace("; -", "")
    out_line= out_line.replace("-", "")
    file2.write(out_line)