#!/usr/bin/env python3
import os
import sys
import argparse
import telebot
from os import walk
import pandas as pd
parser = argparse.ArgumentParser()
parser.add_argument("--mode", help="Action mode", default=1, type=int)
args = parser.parse_args()

Output_dir_files = pd.DataFrame(next(walk("../"), (None, None, []))[2], columns={"file"})
gVCF_files = Output_dir_files.loc[Output_dir_files.file.str.contains(".BSQR.bam")].copy()

if args.mode == 1:
    for i in gVCF_files.file.sort_values()[0:9]:
        print(i)
        print("$GATK Mutect2 --verbosity WARNING -R ../genome.fasta -I ../"+i+" --native-pair-hmm-threads 40 -germline-resource GATK_resources/af-only-gnomad.hg38.vcf.gz -pon GATK_resources/1000g_pon.hg38.vcf.gz --f1r2-tar-gz "+i.replace(".ALN.BSQR.bam", "")+"_f1r2.tar.gz -O "+i.replace(".ALN.BSQR.bam", "")+".unfiltered.vcf.gz")
        os.system("$GATK Mutect2 --verbosity WARNING -R ../genome.fasta -I ../"+i+" --native-pair-hmm-threads 40 -germline-resource GATK_resources/af-only-gnomad.hg38.vcf.gz -pon GATK_resources/1000g_pon.hg38.vcf.gz --f1r2-tar-gz "+i.replace(".ALN.BSQR.bam", "")+"_f1r2.tar.gz -O "+i.replace(".ALN.BSQR.bam", "")+".unfiltered.vcf.gz")

if args.mode == 2:
    for i in gVCF_files.file.sort_values()[9:20]:
        print(i)
        print("$GATK Mutect2 --verbosity WARNING -R ../genome.fasta -I ../"+i+" --native-pair-hmm-threads 40 -germline-resource GATK_resources/af-only-gnomad.hg38.vcf.gz -pon GATK_resources/1000g_pon.hg38.vcf.gz --f1r2-tar-gz "+i.replace(".ALN.BSQR.bam", "")+"_f1r2.tar.gz -O "+i.replace(".ALN.BSQR.bam", "")+".unfiltered.vcf.gz")
        os.system("$GATK Mutect2 --verbosity WARNING -R ../genome.fasta -I ../"+i+" --native-pair-hmm-threads 40 -germline-resource GATK_resources/af-only-gnomad.hg38.vcf.gz -pon GATK_resources/1000g_pon.hg38.vcf.gz --f1r2-tar-gz "+i.replace(".ALN.BSQR.bam", "")+"_f1r2.tar.gz -O "+i.replace(".ALN.BSQR.bam", "")+".unfiltered.vcf.gz")

if args.mode == 3:
    for i in gVCF_files.file.sort_values()[20:33]:
        print(i)
        print("$GATK Mutect2 --verbosity WARNING -R ../genome.fasta -I ../"+i+" --native-pair-hmm-threads 40 -germline-resource GATK_resources/af-only-gnomad.hg38.vcf.gz -pon GATK_resources/1000g_pon.hg38.vcf.gz --f1r2-tar-gz "+i.replace(".ALN.BSQR.bam", "")+"_f1r2.tar.gz -O "+i.replace(".ALN.BSQR.bam", "")+".unfiltered.vcf.gz")
        os.system("$GATK Mutect2 --verbosity WARNING -R ../genome.fasta -I ../"+i+" --native-pair-hmm-threads 40 -germline-resource GATK_resources/af-only-gnomad.hg38.vcf.gz -pon GATK_resources/1000g_pon.hg38.vcf.gz --f1r2-tar-gz "+i.replace(".ALN.BSQR.bam", "")+"_f1r2.tar.gz -O "+i.replace(".ALN.BSQR.bam", "")+".unfiltered.vcf.gz")

#os.system("$GATK GenomicsDBImport -R ../genome.fasta -L ../intervals.interval_list --genomicsdb-workspace-path pon_db -V normal1.vcf.gz -V normal2.vcf.gz -V normal3.vcf.gz")