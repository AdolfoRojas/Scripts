#!/usr/bin/env python3
import os
import sys
import argparse
import telebot
from os import walk
TOKEN = os.environ.get("TOKEN")
tb = telebot.TeleBot(TOKEN)
chatid = os.environ.get("chatID")
chatid_2 = os.environ.get("Vinicius_chatID")
#os.system("scp (archivo) leticia@:(ruta)/(archivo")

parser = argparse.ArgumentParser()
#parser.add_argument("echo", help="echo the string you use here")
parser.add_argument("--R1", help="Read 1 input file")
parser.add_argument("--R2", help="Read 2 input file")
parser.add_argument("--O", help="Output directory")
parser.add_argument("--Ref", help="Reference Genome")
parser.add_argument("--T", help="Threads number", default="6", type=str)
parser.add_argument("--Trim_poly_G", help = "Trim poly-G when using NextSeq/NovaSeq", action="store_true")
parser.add_argument("--Call", help = "Run GATK tools in pipeline", action="store_true")
args = parser.parse_args()
print(args.T)
print(args.R1)
print(args.R2)
print(args.O)
sample_name = args.R1.split("/")[-1].split("_1.fastq")[0]
tb.send_message(chatid, "Pipeline BWA-GATK iniciado en muestra: " + sample_name)
#tb.send_message(chatid_2, "Pipeline BWA-GATK iniciado en muestra: " + sample_name)
GATK = "java -jar /media/storage2/software/gatk-4.1.8.1/gatk-package-4.1.8.1-local.jar"
Genoma = args.O + "genome.fasta"  
Picard = "java -jar /media/storage2/software/Picard/picard.jar"
print("fastqc -t "+ args.T + " " + args.R1 +" "+ args.R2 + " -o " + args.O + "fastqc/ -q")

if sample_name == "output_CS15":
    os.system("mkdir -p " + args.O + "fastQC1")
    os.system("mkdir -p " + args.O + "fastQC2")
    if os.path.isfile(args.O + "fastQC2/trimmed_R1_" + sample_name + "_fastqc.html") == False and os.path.isfile(args.O + "fastQC2/trimmed_R2_" + sample_name + "_fastqc.html") == False:
        os.system("fastqc -t 20 " + args.R1 +" "+ args.R2 + " -o " + args.O + "fastQC1/ -q")
        if os.path.isfile(args.O + "fastQC1/" + sample_name + "_2_fastqc.html") == True and os.path.isfile(args.O + "fastQC1/" + sample_name + "_1_fastqc.html") == True:
            os.system("mkdir -p " + args.O + "Trimmed")
            trim_samp1 = "trimmed_R1_"+sample_name+".fastq.gz"
            trim_samp2 = "trimmed_R2_"+sample_name+".fastq.gz"
            os.system("fastp --in1 "+args.R1+" --in2 "+args.R2+" --adapter_sequence=GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGCTCATTATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAA --unpaired1 " + args.O + "Trimmed/un_1_"+sample_name+".fastq.gz —unpaired2 " + args.O + "Trimmed/un_2_"+sample_name+".fastq.gz --out1 " + args.O + "Trimmed/" + trim_samp1 + " --out2 " + args.O + "Trimmed/" + trim_samp2 + " --thread 16 --length_required 40 -h " + args.O + "Trimmed/report_trim_"+sample_name+".html")
            args.R1 = args.O + "Trimmed/" + trim_samp1
            args.R2 = args.O + "Trimmed/" + trim_samp2
            print(args.R1)
            print(args.R1)
            os.system("fastqc -t 20 " + args.R1 + " " + args.R2 + " -o " + args.O + "fastQC2/ -q")

### Trimming de colas de poliG generadas por **NovaSeq y NextSeq**
#    print("fastp --trim_poly_g --in1 " + args.R1 + " --in2 " + args.R2 + " --out1 Trimmed_" + args.R1 + " --out2 Trimmed_" + args.R2)

#trimmomatic PE                                               #fastp
#-threads                                                       # --thread
#-phred33                                                     # por defecto
#R1.fastq R2.fastq                                            # --in1 --in2
#R1_trimmomatic.fastq R1_trimmomatic_unpaired.fastq            # --out1    --unpaired1
#R2_trimmomatic.fastq R2_trimmomatic_unpaired.fastq            # --out2    --unpaired2
#ILLUMINACLIP:contaminants.fasta:2:30:10                       # --detect_adapter_for_pe
#LEADING:5                                                       # --cut_front --cut_front_window_size 1 --cut_front_mean_quality 5
#TRAILING:5                                                    # --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 5
#SLIDINGWINDOW:4:28                                           # --cut_right --cut_right_window_size 4 --cut_right_mean_quality 28
#MINLEN:$minLen                                                # -l or --length_required 40
#$ trimmomaticPE -threads $cores -phred33 R1.fastq R2.fastq R1_trimmomatic.fastq R1_trimmomatic _unpaired.fastq  R2_trimmomatic.fastq R2_trimmomatic _unpaired.fastq ILLUMINACLIP:contaminants.fasta:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:$quality MINLEN:$minLen
if os.path.isfile(Genoma) == False:
    os.system("cp "+ args.Ref + " " + Genoma)
    print("bwa index -a bwtsw " + args.Ref +" -p " + Genoma)
    os.system("bwa index -a bwtsw " + args.Ref +" -p " + Genoma)    
    #print("picard-tools CreateSequenceDictionary R= " + args.Ref + " O= " + Genoma + ".dict")
    #os.system("picard-tools CreateSequenceDictionary R= " + args.Ref + " O= " + Genoma + ".dict")
    os.system("samtools faidx " + Genoma)
    os.system(GATK + " CreateSequenceDictionary -R " + Genoma) ### Generar diccionario del genoma

if os.path.isfile(args.O + sample_name + "_RG.bam") == False:
    ## **Alineamiento a genoma de referencia**
    print("bwa mem -v 1 -M -t " + args.T + " " + Genoma + " " + args.R1 + " " + args.R2 + " > " + args.O + "aligned_" + sample_name + ".sam")
    os.system("bwa mem -v 1 -M -t " + args.T + " " + Genoma + " " + args.R1 + " " + args.R2 + " > " + args.O + "aligned_" + sample_name + ".sam")
    ### Convertir archivo .sam a .bam
    print("samtools view -S -b " + args.O + "aligned_" + sample_name + ".sam > " + args.O + "aligned_" + sample_name + ".bam")
    os.system("samtools view -S -b " + args.O + "aligned_" + sample_name + ".sam > " + args.O + "aligned_" + sample_name + ".bam")    
    ### Ordenar lecturas por posición (Sort by coordinate)
    print(Picard + " SortSam VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate INPUT= " + args.O + "aligned_" + sample_name + ".bam OUTPUT= " + args.O + "aligned_" + sample_name + "_sorted.bam")
    os.system(Picard + " SortSam VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate INPUT= " + args.O + "aligned_" + sample_name + ".bam OUTPUT= " + args.O + "aligned_" + sample_name + "_sorted.bam")
    ### Marcar duplicados (Mark duplicates)
    print(Picard + " MarkDuplicates VALIDATION_STRINGENCY=SILENT I= " + args.O + "aligned_" + sample_name + "_sorted.bam O= " + args.O + "aligned_dups_removed_" + sample_name + ".bam REMOVE_DUPLICATES=true M=metrics")
    os.system(Picard + " MarkDuplicates VALIDATION_STRINGENCY=SILENT I= " + args.O + "aligned_" + sample_name + "_sorted.bam O= " + args.O + "aligned_dups_removed_" + sample_name + ".bam REMOVE_DUPLICATES=true M=metrics")
    ### Añadir campo ReadGroup (AddOrReplaceReadGroups)
    print(Picard + " AddOrReplaceReadGroups VALIDATION_STRINGENCY=SILENT I= " + args.O + "aligned_dups_removed_" + sample_name + ".bam O= " + args.O + sample_name + "_RG.bam SO=coordinate RGLB=lib_1 RGPL=illumina RGPU=barcode_1 RGSM=sample_1 CREATE_INDEX=true")
    os.system(Picard + " AddOrReplaceReadGroups VALIDATION_STRINGENCY=SILENT I= " + args.O + "aligned_dups_removed_" + sample_name + ".bam O= " + args.O + sample_name + "_RG.bam SO=coordinate RGLB=lib_1 RGPL=illumina RGPU=barcode_1 RGSM=" + sample_name + " CREATE_INDEX=true")


## **Recalibración de bases**
### Descargar base de datos dbSNP desde NCBI (versión GRCh38)
#print("wget https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.38.gz")
### Calcular recalibración de bases
#print("java -jar /media/storage2/software/gatk-4.1.8.1/gatk-package-4.1.8.1-local.jar BaseRecalibrator -R " + args.O + "genome -I " + args.O + sample_name + "_RG.bam -O " + args.O + sample_name + ".ALN.grp")  # --known-sites GCF_000001405.38
#os.system("java -jar /media/storage2/software/gatk-4.1.8.1/gatk-package-4.1.8.1-local.jar BaseRecalibrator -R " + args.O + "genome -I " + args.O + sample_name + "_RG.bam -O " + args.O + sample_name + ".ALN.grp")
### Aplicar recalibración de bases a los datos 
#print("java -jar /media/storage2/software/gatk-4.1.8.1/gatk-package-4.1.8.1-local.jar ApplyBQSR -R " + args.Ref + " -I " + args.O + sample_name + "_RG.bam --bqsr-recal-file " + args.O + sample_name + ".ALN.grp -O " + args.O + sample_name + ".ALN.BSQR.bam")
#os.system("java -jar /media/storage2/software/gatk-4.1.8.1/gatk-package-4.1.8.1-local.jar ApplyBQSR -R " + args.Ref + " -I " + args.O + sample_name + "_RG.bam --bqsr-recal-file " + args.O + sample_name + ".ALN.grp -O " + args.O + sample_name + ".ALN.BSQR.bam")

if args.Call == True:
    if os.path.isfile(args.O + sample_name + "_filtered_INDELs.vcf") == False:
        ## **Llamado de variantes**
        print(GATK + " HaplotypeCaller -R " + Genoma + " -I " + args.O + sample_name + "_RG.bam -O " + args.O + sample_name + ".ALIGNER.HC.raw.vcf -stand-call-conf 50 --native-pair-hmm-threads " + args.T) # -dbsnp GCF_000001405.38
        os.system(GATK + " HaplotypeCaller -R " + Genoma + " -I " + args.O + sample_name + "_RG.bam -O " + args.O + sample_name + ".ALIGNER.HC.raw.vcf -stand-call-conf 50 --native-pair-hmm-threads " + args.T) 
        ### Seleccionar SNPs 
        print(GATK + " SelectVariants -R  " + Genoma + " -V " + args.O + sample_name + ".ALIGNER.HC.raw.vcf -select-type SNP -O " + args.O + sample_name + ".ALIGNER.HC.raw.SNPs.vcf")
        os.system(GATK + " SelectVariants -R  " + Genoma + " -V " + args.O + sample_name + ".ALIGNER.HC.raw.vcf -select-type SNP -O " + args.O + sample_name + ".ALIGNER.HC.raw.SNPs.vcf")
        ### Filtrar SNPs usando parámetros estándar de GATK
        print(GATK + " VariantFiltration -R " + Genoma + " -V " + args.O + sample_name + ".ALIGNER.HC.raw.SNPs.vcf --filter-expression \"QD < 2.0 || FS > 60.0 || MQ <40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" --filter-name \"1_snp_filter\" -O " + args.O + sample_name + "_filtered_SNPs.vcf")
        os.system(GATK + " VariantFiltration -R " + Genoma + " -V " + args.O + sample_name + ".ALIGNER.HC.raw.SNPs.vcf --filter-expression \"QD < 2.0 || FS > 60.0 || MQ <40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" --filter-name \"1_snp_filter\" -O " + args.O + sample_name + "_filtered_SNPs.vcf")
        ### Seleccionar InDels
        print(GATK + " SelectVariants -R " + Genoma + " -V " + args.O + sample_name + ".ALIGNER.HC.raw.vcf -select-type INDEL -O " + args.O + sample_name + ".ALIGNER.HC.raw.INDELs.vcf")
        os.system(GATK + " SelectVariants -R " + Genoma + " -V " + args.O + sample_name + ".ALIGNER.HC.raw.vcf -select-type INDEL -O " + args.O + sample_name + ".ALIGNER.HC.raw.INDELs.vcf")
        ### Filtrar InDels usando parámetros estándar de GATK
        print(GATK + " VariantFiltration -R " + Genoma + " -V " + args.O + sample_name + ".ALIGNER.HC.raw.INDELs.vcf --filter-expression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\" --filter-name \"1_indel_filter\" -O " + args.O + sample_name + "_filtered_INDELs.vcf")
        os.system(GATK + " VariantFiltration -R " + Genoma + " -V " + args.O + sample_name + ".ALIGNER.HC.raw.INDELs.vcf --filter-expression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\" --filter-name \"1_indel_filter\" -O " + args.O + sample_name + "_filtered_INDELs.vcf")

os.system("rm " + args.O + "aligned*")
filenames = next(walk(args.R2.split(sample_name)[0]), (None, None, []))[2]
filenames2 = next(walk(args.O), (None, None, []))[2]
totales = str(sum('2.fastq' in s for s in filenames))
terminados = str(sum('_RG.bai' in s for s in filenames2))

tb.send_message(chatid, sample_name + " Finalizado " + terminados + " de " + totales)
#tb.send_message(chatid_2, sample_name + " Finalizado " + terminados + " de " + totales)

#  ./Variant_calling_bwa-mem_GATK.py 
# --Ref /media/storage/vinicius/vitis/chardonnay/VvChar04_v1.fasta
# --R1 /media/storage/vinicius/vitis/clones/discoVini/Chardonnay/output_CH1_1.fastq.gz 
# --R2 /media/storage/vinicius/vitis/clones/discoVini/Chardonnay/output_CH1_2.fastq.gz 
# --O chardonnay/ 
# --T 20