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
args = parser.parse_args()

print(args.R1)
print(args.R2)

sample_name = args.R1.split("/")[-1].split("_R1_001.fastq")[0]
tb.send_message(chatid, "Trimming iniciado en muestra: " + sample_name)
#tb.send_message(chatid_2, "Trimming iniciado en muestra: " + sample_name)

if os.path.isfile(args.O + "fastQC2/trimmed_R1_" + sample_name + "_fastqc.html") == False and os.path.isfile(args.O + "fastQC2/trimmed_R2_" + sample_name + "_fastqc.html") == False:
    os.system("fastqc -t 20 " + args.R1 +" "+ args.R2 + " -o " + args.O + "fastQC1/ -q")
    if os.path.isfile(args.O + "fastQC1/" + sample_name + "_R1_001_fastqc.html") == True and os.path.isfile(args.O + "fastQC1/" + sample_name + "_R2_001_fastqc.html") == True:
        os.system("fastp --in1 "+args.R1+" --in2 "+args.R2+" --trim_poly_g --adapter_sequence=GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTGTGAAATCTCGTAT --unpaired1 un_1_"+sample_name+".fastq.gz —unpaired2 un_2_"+sample_name+".fastq.gz --out1 trimmed_R1_"+sample_name+".fastq.gz --out2 trimmed_R2_"+sample_name+".fastq.gz --thread 16 --length_required 40 -h report_trim_"+sample_name+".html")
        #os.system("fastp --in1 "+args.R1+" --in2 "+args.R2+" --trim_poly_g --trim_poly_x --adapter_sequence=GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTGTGAAATCTCGTAT --unpaired1 un_1_"+sample_name+".fastq.gz —unpaired2 un_2_"+sample_name+".fastq.gz --out1 trimmed_R1_"+sample_name+".fastq.gz --out2 trimmed_R2_"+sample_name+".fastq.gz --thread 16 --length_required 40 -h report_trim_"+sample_name+".html")
        trim_samp1 = "trimmed_R1_"+sample_name+".fastq.gz"
        trim_samp2 = "trimmed_R2_"+sample_name+".fastq.gz"
        os.system("fastqc -t 20 " + trim_samp1 +" "+ trim_samp2 + " -o " + args.O + "fastQC2/ -q")


#os.system("fastp --in1 "+args.R1+" --in2 "+args.R2+" --trim_poly_g                # inputs y argumento para habilitar remocion de polyG             
#--adapter_sequence=GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTGTGAAATCTCGTAT     ### adaptador encontrado a remover
#--unpaired1 un_1 "+sample_name+".fastq.gz —unpaired2 un_2 "+sample_name+".fastq.gz # Reads huerfanos
#--out1 trimmed_R1_"+sample_name+".fastq.gz --out2 trimmed_R2_"+sample_name+".fastq.gz  # fastq outputs
#--thread 16 --length_required 40 -h report_trim_"+sample_name+".html")  # nucleos, reportes y tamaño minimo de read