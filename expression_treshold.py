#!/usr/bin/env python3
import pandas as pd  
input_file = input("Full file name\n")
Count_file = pd.read_csv(input_file, sep="\t", header=0, index_col= 0).drop(columns="gene_id")
drop_rows = []
for INDEX in Count_file.index:
    Boolean = Count_file.loc[[INDEX]].ge(0.1).sum(axis=1) >= len(Count_file.columns)*0.2
    if Boolean[0] == False:
        drop_rows.append(INDEX)


print(str(len(drop_rows)) + " Transcritos eliminados")
Count_file = Count_file.drop(drop_rows)
print(str(len(Count_file.index)) + " Transcritos se mantienen")
Count_file.reset_index(inplace=True)
Count_file = Count_file.rename(columns = {'index':'transcript_id'})
Count_file['transcript_id'].to_csv("expressed_transcripts.txt", sep="\t",header=True, index=False)
drop_rows = pd.DataFrame(drop_rows,columns=['transcript_id'])
drop_rows.to_csv("removed_transcripts.txt", sep="\t",header=True, index=False)