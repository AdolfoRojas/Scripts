# **Llamado de variantes desde SRA a GATK 4**

### Descargar secuencias de SRA en servidor
~~~
fastq-dump nombrerunsra --split-files  ###divide en read 1 y 2 para el análisis de control de calidad 
~~~
## **Control de Calidad**
### Control de calidad con Fast QC 
~~~
fastqc archivo.fastq -o directorio_salida
~~~
### Trimming de colas de poliG generadas por **NovaSeq y NextSeq**
~~~
fastp --trim_poly_g --in1 *_1.fastq --in2 *_2.fastq --out1 t*_1.fastq --out2 t*_2.fastq
~~~
### Control de calidad post trimming 
~~~
fastqc t*.fastq -o directorio_salida
~~~
## **Genoma de referencia** 
### GRCh38 desde NCBI
~~~
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz

gzip -d GRCh38_latest_genomic.fna.gz 
~~~
### Generar index del genoma
~~~
bwa index GRCh38_latest_genomic.fa
~~~
### Generar diccionario del genoma
~~~
picard-tools CreateSequenceDictionary R=GRCh38_latest_genomic.fna O=genome.dict
~~~
### Concatenar las corridas de cada individuo por R1 y R2
~~~
cat t*_1.fastq t*_1.fastq > allruns_muestra_1.fastq
cat t*_2.fastq t*_2.fastq > allruns_muestra_2.fastq
~~~
## **Alineamiento a genoma de referencia**
~~~
bwa mem -v 1 -M -t 30 GRCh38_latest_genomic.fa directorio/allruns_muestra_1.fastq directorio/allruns_muestra_1.fastq > directorio/aligned_muestra.sam
~~~
### Convertir archivo .sam a .bam
~~~
samtools view -S -b aligned_muestra.sam > aligned_muestra.bam
~~~
### Ordenar lecturas por posición (Sort by coordinate)
~~~
java -Xmx2g -Djava.io.tmpdir='pwd'/tmp -jar/media/storage2/software/Picard/picard.jar SortSam VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate INPUT=aligned_muestra.bam OUTPUT=aligned_muestra_HG01488.bam TMP_DIR='pwd'/tmp
~~~
### Marcar duplicados (Mark duplicates)
~~~
java -jar /media/storage2/software/Picard/picard.jar MarkDuplicates VALIDATION_STRINGENCY=SILENT I=aligned_sorted_muestra.bam O=aligned_dups_removed_muestra.bam REMOVE_DUPLICATES=true M=metrics
~~~
### Añadir campo ReadGroup (AddOrReplaceReadGroups)
~~~
java -jar /media/storage2/software/Picard/picard.jar AddOrReplaceReadGroups VALIDATION_STRINGENCY=SILENT I=aligned_dups_removed2_muestra.bam O=muestra_RG.bam SO=coordinate RGLB=lib_1 RGPL=illumina RGPU=barcode_1 RGSM=sample_1 CREATE_INDEX=true
~~~
## **Recalibración de bases**

### Descargar base de datos dbSNP desde NCBI (versión GRCh38)
~~~
wget https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.38.gz
~~~
### Calcular recalibración de bases
~~~
java -jar /media/storage2/software/gatk-4.1.8.1/gatk-package-4.1.8.1-local.jar BaseRecalibrator -R GRCh38_latest_genomic.fa -I /directorio/muestra_RG.bam --known-sites GCF_000001405.38 -O muestra.ALN.grp 
~~~
### Aplicar recalibración de bases a los datos 
~~~
java -jar /media/storage2/software/gatk-4.1.8.1/gatk-package-4.1.8.1-local.jar ApplyBQSR -R GRCh38_latest_genomic.fa -I /directorio/muestra_RG.bam --bqsr-recal-file muestra.ALN.grp -O muestra.ALN.BSQR.bam
~~~
## **Llamado de variantes**
~~~
java -jar /media/storage2/software/gatk-4.1.8.1/gatk-package-4.1.8.1-local.jar HaplotypeCaller -R GRCh38_latest_genomic.fa -I muestra.ALN.BQSR.bam -dbsnp GCF_000001405.38 -O muestra.ALIGNER.HC.raw.vcf -stand-call-conf 50
~~~
### Filtrado de variantes usando archivo .bed para regiones codificantes y lncRNA de interés (archivo modificado con notación NC).

~~~
vcftools --vcf muestra.ALIGNER.HC.raw.vcf --bed coding_lncRNA.bed --out muestra_lnc_cod.vcf --recode --keep-INFO-all
~~~

### Seleccionar SNPs 
~~~
java -jar /media/storage2/software/gatk-4.1.8.1/gatk-package-4.1.8.1-local.jar SelectVariants -R GRCh38_latest_genomic.fa -V muestra_lnc_cod.vcf.recode.vcf -select-type SNP -O muestra_lnc_cod_snp.vcf
~~~
### Filtrar SNPs usando parámetros estándar de GATK
~~~
java -jar /media/storage2/software/gatk-4.1.8.1/gatk-package-4.1.8.1-local.jar VariantFiltration -R GRCh38_latest_genomic.fa -V muestra_lnc_cod_snp.vcf --filter-expression "QD < 2.0 || FS > 60.0 || MQ <40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "1_snp_filter" -O muestra_filtered_snps.vcf
~~~
### Seleccionar InDels
~~~
java -jar /media/storage2/software/gatk-4.1.8.1/gatk-package-4.1.8.1-local.jar SelectVariants -R GRCh38_latest_genomic.fa -V muestra_lnc_cod.vcf.recode.vcf -select-type INDEL -O muestra_lnc_cod_indels.vcf
~~~
### Filtrar InDels usando parámetros estándar de GATK
~~~
java -jar /media/storage2/software/gatk-4.1.8.1/gatk-package-4.1.8.1-local.jar VariantFiltration -R GRCh38_latest_genomic.fa -V muestra_lnc_cod_indels.vcf --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filter-name "1_indel_filter" -O muestra_filtered_indels.vcf
~~~

## **Anotación de variantes** 

### Descargar base de datos de snpEff
~~~
java -jar /media/storage2/software/snpEff/snpEff.jar download -v hg38 
~~~
### Anotación de variantes SNPs (archivo .vcf y .bed con notación de chr)
~~~
va -jar /media/storage2/software/snpEff/snpEff.jar hg38 -interval coding_and_lncRNA_original.bed muestra_filtered_snp.vcf | vcf-sort > muestra_annotated_snpEff_filtered_snps.vcf
~~~
### Anotación de variantes InDels (archivo .vcf y .bed con notación de chr)
~~~
java -jar /media/storage2/software/snpEff/snpEff.jar hg38 -interval coding_lncRNA_original.bed hg38 muestra_filtered_indels.vcf | vcf-sort > muestra_annotated_snpEff_filtered_indels.vcf
~~~
### Anotación de variantes por Clinvar

### Descargar base de datos clinvar para GRCh38
~~~
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.md5
~~~
### Anotación de variantes SNPs por clinvar
~~~
java -jar /media/storage2/software/snpEff/SnpSift.jar annotate clinvar.vcf.gz muestra_annotated_snpEff_filtered_snps.vcf > muestra_annotated_clinvar_filtered_snps.vcf
~~~
### Anotación de variantes InDels por clinvar
~~~
java -jar /media/storage2/software/snpEff/SnpSift.jar annotate clinvar.vcf.gz muestra_annotated_snpEff_filtered_indels.vcf > muestra_annotated_clinvar_filtered_indels.vcf
~~~