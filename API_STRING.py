#!/usr/bin/env python3.8  
##################################################################
import sys
sys.path.append("/usr/local/lib/python2.7/dist-packages")
from Bio import Entrez
import re 
import requests ## python -m pip install requests
from gprofiler import GProfiler
##################################################################

rs_ID = input("Inserte rs: ")
### a rs10953105 tiene 2 genes uno codificante y el otro no codificante
Entrez.email = "adolfo.rh@postqyf.uchile.cl" 
handle = Entrez.efetch(db="snp", id=rs_ID, retmode="text")
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
if len(GENE_ID) == 2:
    print("\t".join([Variant_Type, GENE_ID[0], NAME[0], GENE_ID[1], NAME[1]]))
elif len(GENE_ID) == 1:
    print("\t".join([Variant_Type, GENE_ID[0], NAME[0]]))
elif len(GENE_ID) > 2:
    print("Revisar script, aumentar capacidad de busquenda a " + len(GENE_ID))
elif len(GENE_ID) == 0:
    print(rs_ID + " no presenta un gen asociado en dbSNP buscar en otras bases de datos")
handle.close() 

##############################  STRING  #########################################

if len(NAME) != 0:
   
    string_api_url = "https://string-db.org/api"
    output_format = "tsv"
    method = "interaction_partners" 

    request_url = "/".join([string_api_url, output_format, method]) ## Construct URL

    my_genes = list()
    params = {}
    my_genes = NAME[0]
    params = {
    "identifiers" : my_genes, # your protein
    "species" : 9606, # species NCBI identifier 
    "required_score" : 700
    }
    if len(NAME) == 2:
        my_genes = [NAME[0], NAME[1]]
        params = {
        "identifiers" : "%0d".join(my_genes), # your protein
        "species" : 9606, # species NCBI identifier 
        "required_score" : 700
        }
                   
    print(my_genes)
    print(len(NAME))
 
    response = requests.post(request_url, data=params) ## Call STRING
    
    Interactors= list()

    for line in response.text.strip().split("\n"):
        try:
            l = line.strip().split("\t")
            query_ensp = l[0]
            query_name = l[2]
            partner_ensp = l[1]
            partner_name = l[3]
            combined_score = l[5]
            print("\t".join([query_ensp, query_name, partner_name, combined_score]))
            if partner_name != None:
                Interactors.append(partner_name)
        except:
            pass
            print("Not found in String")
    print (len(Interactors[1:]))

    with open(rs_ID+'_interactors.txt', 'w') as interactors:
        for item in Interactors[1:]:
            interactors.write("%s\n" % item)
    interactors.close()

    query_range = list()

    pathways_interactors = list()

    for i in range(1, len(Interactors)):
        query_range.append("Query"+str(i))
    
    d = dict( zip( Interactors[1:], Interactors[1:]))
    pathways_interactors_out= open(rs_ID+".tsv","w+")
    gp = GProfiler(return_dataframe=True)
    pathways_interactors= gp.profile(organism = 'hsapiens',
            query = Interactors[1:], combined = False, no_iea = True, ordered = True, sources = ["GO:MF","REAC"])     
    pathways_interactors.to_csv(pathways_interactors_out, sep = "\t")
    pathways_interactors_out.close()

    #pathways_interactors_out2= open(rs_ID+"_full.tsv","w+")
    #gp = GProfiler(return_dataframe=True)
    #pathways_interactors2= gp.profile(organism = 'hsapiens',
    #        query = Interactors[1:])
    #pathways_interactors2.to_csv(pathways_interactors_out2, sep = "\t")
    #pathways_interactors_out2.close()



    