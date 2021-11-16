#!/usr/bin/env python3
import os
import sys
import argparse
import telebot
from os import walk
import pandas as pd
TOKEN = os.environ.get("TOKEN")
tb = telebot.TeleBot(TOKEN)
chatid = os.environ.get("chatID")
chatid_2 = os.environ.get("Vinicius_chatID")

parser = argparse.ArgumentParser()
parser.add_argument("--R1", help="Read 1 input file")
parser.add_argument("--R2", help="Read 2 input file")
parser.add_argument("--O", help="Output directory")
parser.add_argument("--Ref", help="Reference Genome")
parser.add_argument("--T", help="Threads number", default="6", type=str)
parser.add_argument("--QC1", help = "Do first fastQC", action="store_true")
parser.add_argument("--Trim", help = "Trim fastq", action="store_true")
parser.add_argument("--Trim_adapters", help="Auto trimming adapters", action="store_true")
parser.add_argument("--Trim_poly_G", help = "Trim poly-G when using NextSeq/NovaSeq", action="store_true")
parser.add_argument("--Trim_quality", help="Trim bad quality Reads", action="store_true")
parser.add_argument("--Trim_fasta_sequences", help="Adds fasta sequences to trim", default="", type=str)
parser.add_argument("--Map", help = "Map to genome", action="store_true")
parser.add_argument("--Homo_sapiens", help="Run BQSR/VQSR, need knowsites input file (Homo sapiens only)", action="store_true")
parser.add_argument("--Call", help = "Run GATK tools in pipeline", action="store_true")
parser.add_argument("--Ploidy", help="Ploidy to call", default="2", type=str)
parser.add_argument("--Cohort_mode", help = "Run Joint Call and downstream tools", action="store_true")
args = parser.parse_args()
print(args.T)
print(args.R1)
print(args.R2)
print(args.O)

GATK = "java -jar /media/storage2/software/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar"
Picard = "java -jar /media/storage2/software/Picard/picard.jar"
Genoma = args.O + "genome.fasta"

if args.Cohort_mode != True:     
    sample_name = args.R1.split("/")[-1].split("_1.fastq")[0]
    print("fastqc -t "+ args.T + " " + args.R1 +" "+ args.R2 + " -o " + args.O + "fastQC1/ -q")

    if args.Trim_fasta_sequences != "":
        args.Trim_fasta_sequences = " --adapter_fasta " + args.Trim_fasta_sequences
    polyG = ""
    if args.Trim_poly_G == True:
        polyG = " --trim_poly_g"
    adapters =""
    if args.Trim_adapters == True:
        adapters = " --detect_adapter_for_pe"
    quality_trim = ""
    if args.Trim_quality == True:
        quality_trim = " --cut_right --cut_right_window_size 2 --cut_right_mean_quality 30 --trim_tail1=76 --trim_tail2=76 --trim_front1=16 --trim_front2=16" #--cut_front --cut_front_window_size 1 --cut_front_mean_quality 28 --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 30
        #

    if args.QC1 == True:
        os.system("mkdir -p " + args.O + "fastQC1")    
        if os.path.isfile(args.O + "fastQC2/trimmed_R1_" + sample_name + "_fastqc.html") == False and os.path.isfile(args.O + "fastQC2/trimmed_R2_" + sample_name + "_fastqc.html") == False:
            os.system("fastqc -t " + args.T + " " + args.R1 +" "+ args.R2 + " -o " + args.O + "fastQC1/ -q")
    if args.QC1 != True:
        if args.Trim == True:
            os.system("mkdir -p " + args.O + "Trimmed")
            os.system("mkdir -p " + args.O + "fastQC2")
            trim_samp1 = "trimmed_R1_"+sample_name+".fastq.gz"
            trim_samp2 = "trimmed_R2_"+sample_name+".fastq.gz"  
            print("fastp --in1 "+args.R1+" --in2 "+args.R2+" --unpaired1 " + args.O + "Trimmed/un_1_"+sample_name+".fastq.gz —unpaired2 " + args.O + "Trimmed/un_2_"+sample_name+".fastq.gz --out1 " + args.O + "Trimmed/" + trim_samp1 + " --out2 " + args.O + "Trimmed/" + trim_samp2 + " --thread 16 --length_required 40 -h " + args.O + "Trimmed/report_trim_"+sample_name+".html" + adapters + polyG + quality_trim + args.Trim_fasta_sequences)       
            os.system("fastp --in1 "+args.R1+" --in2 "+args.R2+" --unpaired1 " + args.O + "Trimmed/un_1_"+sample_name+".fastq.gz —unpaired2 " + args.O + "Trimmed/un_2_"+sample_name+".fastq.gz --out1 " + args.O + "Trimmed/" + trim_samp1 + " --out2 " + args.O + "Trimmed/" + trim_samp2 + " --thread 16 --length_required 40 -h " + args.O + "Trimmed/report_trim_"+sample_name+".html" + adapters + polyG + quality_trim + args.Trim_fasta_sequences)        
            args.R1 = args.O + "Trimmed/" + trim_samp1
            args.R2 = args.O + "Trimmed/" + trim_samp2
            print(args.R1)
            print(args.R1)
            os.system("fastqc -t " + args.T + " " + args.R1 + " " + args.R2 + " -o " + args.O + "fastQC2/ -q")        
        if os.path.isfile(args.O + "Trimmed/" + "trimmed_R1_"+sample_name+".fastq.gz") == True and args.Trim != True:
            trim_samp1 = "trimmed_R1_"+sample_name+".fastq.gz"
            trim_samp2 = "trimmed_R2_"+sample_name+".fastq.gz"
            args.R1 = args.O + "Trimmed/" + trim_samp1
            args.R2 = args.O + "Trimmed/" + trim_samp2
            print(args.R1)
            print(args.R1)        
        if os.path.isfile(args.O + "genome.dict") == False:
            os.system("cp "+ args.Ref + " " + Genoma)
            print("bwa index -a bwtsw " + args.Ref +" -p " + Genoma)
            os.system("bwa index -a bwtsw " + args.Ref +" -p " + Genoma)    
            os.system("samtools faidx " + Genoma)
            os.system(GATK + " CreateSequenceDictionary -R " + Genoma) ### Generar diccionario del genoma
        if args.Map == True:
            #tb.send_message(chatid, "Pipeline BWA-GATK iniciado en muestra: " + sample_name)
            #tb.send_message(chatid_2, "Pipeline BWA-GATK iniciado en muestra: " + sample_name)
            if os.path.isfile(args.O + sample_name + "_RG.bam") == False:
                ## **Alineamiento a genoma de referencia**
                print("bwa mem -v 1 -M -t " + args.T + " " + Genoma + " " + args.R1 + " " + args.R2 + " > " + args.O + "aligned_" + sample_name + ".sam")
                os.system("bwa mem -v 1 -M -t " + args.T + " " + Genoma + " " + args.R1 + " " + args.R2 + " > " + args.O + "aligned_" + sample_name + ".sam")
                ### Convertir archivo .sam a .bam
                print("samtools view -S -b " + args.O + "aligned_" + sample_name + ".sam > " + args.O + "aligned_" + sample_name + ".bam")
                os.system("cat  " + args.O + "aligned_" + sample_name + ".sam > " + args.O + "aligned_" + sample_name + ".bam")    
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
            if args.Homo_sapiens == True:
                if os.path.isfile(args.O + sample_name + ".ALN.BSQR.bam") == False: 
                    ### Descargar base de datos dbSNP desde NCBI (versión GRCh38)
                    #print("wget https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.38.gz")
                    known_sites = " --known-sites /media/storage/datasets/GATK_known_sites/Homo_sapiens/"
                    ### Calcular recalibración de bases
                    print(GATK + " BaseRecalibrator -R " + Genoma + " -I " + args.O + sample_name + "_RG.bam -O " + args.O + sample_name + ".ALN.grp" + known_sites +"00-All.vcf.gz"+known_sites+"Homo_sapiens_assembly38.dbsnp138.vcf.gz"+known_sites+"Homo_sapiens_assembly38.known_indels.vcf.gz"+known_sites+"Mills_and_1000G_gold_standard.indels.hg38.vcf.gz")
                    os.system(GATK + " BaseRecalibrator -R " + Genoma + " -I " + args.O + sample_name + "_RG.bam --verbosity WARNING -O " + args.O + sample_name + ".ALN.grp" + known_sites+"human_9606_b151_GRCh38p7_V2.vcf.gz")
                    ### Aplicar recalibración de bases a los datos 
                    os.system(GATK + " ApplyBQSR -R " + Genoma + " -I " + args.O + sample_name + "_RG.bam --bqsr-recal-file " + args.O + sample_name + ".ALN.grp --verbosity WARNING -O " + args.O + sample_name + ".ALN.BSQR.bam")                
            if args.Call == True:
                if os.path.isfile(args.O + sample_name + ".HC.raw.g.vcf.gz") == False:
                    ## **Llamado de variantes**
                    if os.path.isfile(args.O + sample_name + ".ALN.BSQR.bam") == True:
                        os.system(GATK + " HaplotypeCaller -R " + Genoma + " -I " + args.O + sample_name + ".ALN.BSQR.bam -O " + args.O + sample_name + ".HC.raw.g.vcf.gz -ERC GVCF --verbosity WARNING --native-pair-hmm-threads " + args.T + " --sample-ploidy " + args.Ploidy) # -dbsnp GCF_000001405.38
                    else:
                        os.system(GATK + " HaplotypeCaller -R " + Genoma + " -I " + args.O + sample_name + "_RG.bam -O " + args.O + sample_name + ".HC.raw.g.vcf.gz -ERC GVCF --verbosity WARNING --native-pair-hmm-threads " + args.T + " --sample-ploidy " + args.Ploidy) 
            os.system("rm " + args.O + "aligned*")
            #filenames = next(walk(args.O), (None, None, []))[2]
            filenames2 = next(walk(args.O), (None, None, []))[2]
            #totales = str(sum('2.fastq' in s for s in filenames))
            terminados = str(sum('HC.raw.g.vcf.gz.tbi' in s for s in filenames2))
            #tb.send_message(chatid, sample_name + " Listo, finalizados: " + terminados)
            #tb.send_message(chatid_2, sample_name + " Listo, finalizados: " + terminados)
else:    
    Output_dir_files = pd.DataFrame(next(walk(args.O), (None, None, []))[2], columns={"file"})
    gVCF_files = Output_dir_files.loc[Output_dir_files.file.str.contains("HC.raw.g.vcf.gz")].loc[~Output_dir_files.file.str.contains(".tbi")]
    os.system(Picard + " ScatterIntervalsByNs REFERENCE=" + Genoma+" OUTPUT=" + Genoma + ".interval_list")
    GenomicsDBImport_command = GATK +" GenomicsDBImport --genomicsdb-workspace-path "+ args.O + "GenomicsDB/ --intervals " + Genoma + ".interval_list"
    for i in gVCF_files.file: GenomicsDBImport_command+= " -V " + args.O + i
    os.system(GenomicsDBImport_command)
    raw_joint_vcf = args.O + os.path.basename(os.getcwd()) +"_n" + str(len(gVCF_files.file)) +".GenotypeGVCFs.raw.vcf.gz"
    os.system(GATK + " GenotypeGVCFs -R " + Genoma + " -V gendb://"+ args.O + "GenomicsDB/ --sample-ploidy "+ args.Ploidy +" -O " + raw_joint_vcf)
    if args.Homo_sapiens == True:
        known_sites2 = "/media/storage/datasets/GATK_known_sites/Homo_sapiens/"
        SNPs_resources = " --resource:hapmap,known=false,training=true,truth=true,prior=15.0 " + known_sites2 + "bundle_GATK/hapmap_3.3.hg38.vcf.gz --resource:omni,known=false,training=true,truth=false,prior=12.0 " + known_sites2 + "bundle_GATK/1000G_omni2.5.hg38.vcf.gz --resource:1000G,known=false,training=true,truth=false,prior=10.0 " + known_sites2 + "bundle_GATK/1000G_phase1.snps.high_confidence.hg38.vcf.gz --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 " + known_sites2 + "human_9606_b151_GRCh38p7_V2.vcf.gz"
        INDELs_resources = " --resource:mills,known=false,training=true,truth=true,prior=12.0 " + known_sites2 + "bundle_GATK/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 " + known_sites2 + "human_9606_b151_GRCh38p7_V2.vcf.gz"
        os.system(GATK + " VariantRecalibrator -R " + Genoma + " -V " + raw_joint_vcf + SNPs_resources  + " -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP --max-gaussians 4 -O "+args.O+"output.SNPs.recal --tranches-file "+args.O+"output.SNPs.tranches --rscript-file "+args.O+"output.SNPs.plots.R")
        os.system(GATK + " ApplyVQSR -R " + Genoma + " -V " + raw_joint_vcf + " --truth-sensitivity-filter-level 99.0 --tranches-file "+args.O+"output.SNPs.tranches --recal-file "+args.O+"output.SNPs.recal -mode SNP -O "+args.O+"recalibrated_snps_raw_indels.vcf.gz")
        os.system(GATK + " VariantRecalibrator -R " + Genoma + " -V "+args.O+"recalibrated_snps_raw_indels.vcf.gz" + INDELs_resources  + " -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode INDEL --max-gaussians 4 -O "+args.O+"output.INDELs.recal --tranches-file "+args.O+"output.INDELs.tranches --rscript-file "+args.O+"output.INDELs.plots.R")
        os.system(GATK + " ApplyVQSR -R " + Genoma + " -V "+args.O+"recalibrated_snps_raw_indels.vcf.gz --truth-sensitivity-filter-level 99.0 --tranches-file "+args.O+"output.INDELs.tranches --recal-file "+args.O+"output.INDELs.recal -mode INDEL -O "+args.O + os.path.basename(os.getcwd()) +"_n" + str(len(gVCF_files.file)) +".GenotypeGVCFs.VQSR_filtered.vcf.gz")
        filtered_SNPs = args.O + os.path.basename(os.getcwd()) +"_n" + str(len(gVCF_files.file)) +".GenotypeGVCFs.VQSR_filtered.SNPs.vcf.gz"
        filtered_INDELs = args.O + os.path.basename(os.getcwd()) +"_n" + str(len(gVCF_files.file)) +".GenotypeGVCFs.VQSR_filtered.INDELs.vcf.gz"
        os.system(GATK + " SelectVariants -R  " + Genoma + " -V " + args.O + os.path.basename(os.getcwd()) +"_n" + str(len(gVCF_files.file)) +".GenotypeGVCFs.VQSR_filtered.vcf.gz -select-type SNP -O " + filtered_SNPs)
        os.system(GATK + " SelectVariants -R " + Genoma + " -V " + args.O + os.path.basename(os.getcwd()) +"_n" + str(len(gVCF_files.file)) +".GenotypeGVCFs.VQSR_filtered.vcf.gz -select-type INDEL -O " + filtered_INDELs)
        #$GATK CollectVariantCallingMetrics -I bizamaExoma_n33.GenotypeGVCFs.VQSR_filtered.vcf.gz --DBSNP /media/storage/datasets/GATK_known_sites/Homo_sapiens/human_9606_b151_GRCh38p7_V2.vcf.gz -SD genome.dict -O Variant_calling_metrics
    else:
        raw_joint_SNPs_vcf = args.O + os.path.basename(os.getcwd()) +"_n" + str(len(gVCF_files.file)) +".GenotypeGVCFs.raw.SNPs.vcf.gz"
        raw_joint_INDELs_vcf = args.O + os.path.basename(os.getcwd()) +"_n" + str(len(gVCF_files.file)) +".GenotypeGVCFs.raw.INDELs.vcf.gz"
        filtered_SNPs = args.O + os.path.basename(os.getcwd()) +"_n" + str(len(gVCF_files.file)) +".GenotypeGVCFs.hard_filtered.SNPs.vcf.gz"
        filtered_INDELs = args.O + os.path.basename(os.getcwd()) +"_n" + str(len(gVCF_files.file)) +".GenotypeGVCFs.hard_filtered.INDELs.vcf.gz"
        filtered_vcf = args.O + os.path.basename(os.getcwd()) +"_n" + str(len(gVCF_files.file)) +".GenotypeGVCFs.hard_filtered.all.vcf.gz"
        SNPs_filters = ' -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"'
        INDELs_filters = ' -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "FS > 200.0" --filter-name "FS200" -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"'
        ### Seleccionar SNPs         
        os.system(GATK + " SelectVariants -R  " + Genoma + " -V " + raw_joint_vcf + " -select-type SNP -O " + raw_joint_SNPs_vcf)
        ### Filtrar SNPs usando parámetros estándar de GATK        
        os.system(GATK + " VariantFiltration -R " + Genoma + " -V " + raw_joint_SNPs_vcf + SNPs_filters + " -O " + filtered_SNPs)
        ### Seleccionar InDels        
        os.system(GATK + " SelectVariants -R " + Genoma + " -V " + raw_joint_vcf + " -select-type INDEL -O " + raw_joint_INDELs_vcf)
        ### Filtrar InDels usando parámetros estándar de GATK        
        os.system(GATK + " VariantFiltration -R " + Genoma + " -V " + raw_joint_INDELs_vcf + INDELs_filters + " -O " + filtered_INDELs)
        #os.system(Picard + " MergeVcfs I= "+filtered_INDELs+" I= "+filtered_SNPs+" O= "+filtered_vcf)
    #tb.send_message(chatid, sample_name + " Listo, finalizados: " + terminados)
    #tb.send_message(chatid_2, sample_name + " Listo, finalizados: " + terminados)
    #java -jar /media/storage2/software/snpEff/snpEff.jar -q -noLog -nodownload hg38 bizamaExoma_n33.GenotypeGVCFs.VQSR_filtered.vcf.gz | vcf-sort | gzip > bizamaExoma_n33.GenotypeGVCFs.VQSR_filtered.snpEff.vcf.gz
    #java -jar /media/storage2/software/snpEff/SnpSift.jar annotate -clinvar Annotated_snpEff_filtered_snps2.vcf | gzip > Annotated_snpEff-Clinvar_filtered_snps.vcf.gz
    #java -jar /media/storage2/software/snpEff/snpEff.jar -q -noLog -nodownload hg38 1574_T_FFPE_.FilterMutectCalls.filtered.vcf.gz | vcf-sort | gzip > 1574_T_FFPE_.FilterMutectCalls.filtered.snpEff.vcf.gz &