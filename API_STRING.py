#!/usr/bin/env python3
##################################################################
import sys
sys.path.append("/usr/local/lib/python2.7/dist-packages")
from Bio import Entrez
import re 
import requests ## python -m pip install requests
##################################################################

#####rs_ID = input("Inserte rs: ")
def Gene_anotation():
    file1 = open(input("Lista de SNP: "), 'r')
    file2 = open('Genes_afectados.tsv', 'w+')
    file2.write("Rs\tTipo de Polimorfismo\tN°_de _genes\tGenes afectados\n")
    Lines = file1.readlines()
    Lines = [w.replace('\n', '') for w in Lines]
    print(str(len(Lines)) + " Variantes geneticas a buscar")
    count = 0
    for line in Lines:        
        count+=1
        print(str(count/len(Lines)*100)+"% de Variantes buscadas ("+ str(count) +" de " + str(len(Lines)) + ")")
        Entrez.email = "adolfo.rh@postqyf.uchile.cl" 
        handle = Entrez.efetch(db="snp", id=line, retmode="text")
        Variant_Type = "N.A."
        NAME = list()
        GENE_ID = list()    
        for element in handle.readline().strip().split("<"):
            if re.match("SNP_CLASS>\w+", element):
                Variant_Type = element.replace("SNP_CLASS>", "")
            if re.match("NAME>\w+", element):
                NAME1 = element.replace("NAME>", "")
                NAME.append(NAME1)
            if re.match("GENE_ID>\w+", element):
                GENE_ID1 = element.replace("GENE_ID>", "")
                GENE_ID.append(GENE_ID1)
        if len(GENE_ID) > 0:
            file2.write("\t".join([line, Variant_Type, str(len(GENE_ID)), str(NAME)]) + '\n')                
        elif len(GENE_ID) == 0:
            file2.write(line + "\t" + Variant_Type + "\t" + "¿Intergenico?, no presenta un gen asociado en dbSNP" + '\n')
        handle.close() 
#Gene_anotation()

##############################  STRING  #########################################

def STRING():
    
    file3 = open('Genes_afectados.tsv', 'r')    
    file4 = open('Interacciones_encontradas.tsv', 'w')    
    file4.write("Rs\tN°_Interactores\tInteracciones_sobre_0.7\n") ### Nombre de los encabezados de las columnas
    Lines2 = file3.readlines()[1:]
    print(str(len(Lines2)) + " Variantes geneticas a buscar")
    count2 =0
    for line in Lines2:
        l1 = line.strip().split("\t")
        count2 +=1
        print(str(count2/len(Lines2)*100)+"% de potenciales interacciones buscadas ("+ str(count2) +" de " + str(len(Lines2)) + ")")       
        if len(l1) <2: #len(l1) <4
            continue 
        #Number_Genes = l1[2]  #Number_Genes = l1[2]      
        NAME = l1[1] #NAME = l1[3]
        Rs = l1[0]  #Rs = l1[0]           
        string_api_url = "https://string-db.org/api"
        output_format = "tsv"
        method = "interaction_partners" 

        request_url = "/".join([string_api_url, output_format, method]) ## Construct URL
        NAME = NAME.replace("[", "")
        NAME = NAME.replace("]", "")
        NAME = NAME.replace("'", "")        
        my_genes = NAME.split(", ")
        params = {}        
        params = {
        "identifiers" : "%0d".join(my_genes), # your protein
        "species" : 9606, # species NCBI identifier 
        #"required_score" : 700
        }              
        response = requests.post(request_url, data=params) ## Call STRING
    
        Interactors= list()
        count3 = 0
        interaction = list()
        partner_name = "N.A."        
        for line in response.text.strip().split("\n")[1:]:  #[1:] evita que se impriman los encabezados
            try:
                l = line.strip().split("\t")                   
                query_name = l[2]                                
                partner_name = l[3]                
                combined_score = l[5]
                output = (query_name + "-" + partner_name +"(" + combined_score +")")
                interaction.append(output)                
                if partner_name != None:
                    Interactors.append(partner_name)                               
            except:
                pass
                output = "Genes no encontrados en STRING (¿ncRNAs?)"
                interaction.append(output)  

        #print(Rs + " " + str(len(Interactors)) + " Interactores encontrados" + "\n")
        interaction = str(interaction)
        interaction = interaction.replace("\n", ";") 
        file4.write(Rs + "\t" + str(len(Interactors)) + "\t" + interaction + "\n")
STRING()    