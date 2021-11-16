# **Imputacion de Variantes**

### Filtro de individuos y variantes por missing data se remueven individuos con mas de 1% y variantes con mas de 10%
~~~
plink --bfile dbgap_cimba_b_c1 --geno 0.01 --mind 0.01 --make-bed --out dbgap_cimba_b_c1_non_missing # --geno 0.01 --mind 0.005
~~~
### Separacion en archivos plink por cromosoma
~~~
for chr in {1..22} 
do 
    plink --bfile dbgap_cimba_b_c1_non_missing --chr $chr --make-bed --out dbgap_cimba_chr_${chr} 
    plink --bfile dbgap_cimba_chr_${chr} --list-duplicate-vars ids-only suppress-first 
    plink --bfile dbgap_cimba_chr_${chr} -exclude plink.dupvar --make-bed --out dbgap_cimba_chr_${chr}.DuplicatesRemoved 
done
~~~
### Division de recursos por servidor

#### taruca 1-8 10T 4J
#### sol 15-22 8T 2J
#### jagua 9-14 8T 3J

### Comando para separar archivos en carpetas por servidor
~~~
for chr in {15..22} # Ejemplo Sol
do 
    mv dbgap_cimba_chr_${chr}.DuplicatesRemoved* Sol/
    mv genetic_map_chr${chr}_combined_b37.txt Sol/
done
~~~
### Generacion de comando para SHAPEIT, para la ejecucion en paralelo
~~~
#-M genetic_map_chr${chr}_combined_b37.txt 
for chr in {1..22} # Ejemplo Taruca
do
  echo "-B dbgap_cimba_chr_${chr}.DuplicatesRemoved -M chr${chr}.b37.gmap.gz -O cimba_chr${chr}.phased --states 300 -T 8 --output-log Shapeit_cimba_chr${chr}.phased" >> myCommands.txt
done

cat myCommands.txt | xargs -P3 -n12 shapeit & # P4 indica 4 trabajos a la vez n12 la cantidad de argumentos del comando
rm myCommands.txt
~~~
## **Correccion cromosoma 22** un individuo (IDB001413) presento una missingness superior a 10% para este cromosoma por lo que debe ser removido de este y los otros cromosomas. 
~~~
plink1.9 --bfile dbgap_cimba_chr_22.DuplicatesRemoved --geno --mind 0.1 --make-bed --out dbgap_cimba_chr_22.DuplicatesRemoved_nonmissing

cat *.irem | cut -d$'\t' -f 2 > individuos_excluidos.txt ## IDB001413       IDB001413

shapeit -B dbgap_cimba_chr_22.DuplicatesRemoved_nonmissing -M genetic_map_chr22_combined_b37.txt -O cimba_chr22.phased --states 200 -T 8 --output-log Shapeit_cimba_chr22.phased &sc
~~~
### Convertir SHAPEIT Files a VCF
~~~
for chr in {9..14} # Ejemplo Jagua
do
    shapeit -convert --input-haps  cimba_chr${chr}.phased --exclude-ind individuos_excluidos.txt --output-vcf  cimba_chr${chr}.phased.vcf
done
~~~
## **Minimac 4**
### Descomprimir VCFs de referencia
~~~
for chr in {1..22}
do
    gunzip ${chr}.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz
done
~~~
## Realizar imputacion
~~~
for chr in {1..22}
do
  echo "--refHaps ../m3vcf_files/${chr}.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf --haps cimba_chr${chr}.phased.vcf --prefix cimba_chr${chr}.imputed --cpus 7" >> myCommands.txt
done

cat myCommands.txt | xargs -P4 -n8 minimac4 & 
rm myCommands.txt
~~~
## Transformar a PLINK data
~~~
for chr in {1..22}
do
  awk -F "\t" '{ if(($7 >= 0.3) && ($5 >= 0.01)) { print } }' cimba_chr${chr}.imputed.info | cut -d$'\t' -f 1 | tail -n +2 > good_variants_chr${chr}
  vcftools --gzvcf cimba_chr${chr}.imputed.dose.vcf.gz --snps good_variants_chr${chr} --recode --recode-INFO-all --stdout | vcf-sort > filtered_cimba_chr${chr}.imputed.dose.recode.vcf & 
done

for chr in {1..22}
do
  DosageConvertor --vcfDose filtered_cimba_chr${chr}.imputed.dose.recode.vcf --info cimba_chr${chr}.imputed.info --prefix plink_files/CIMBA_imputed_chr${chr} --type plink --format 1 &
done

for chr in {1..22}
do
  awk '!x[$0]++' CIMBA_imputed_chr${chr}.plink.map > CIMBA_imputed_chr${chr}_v2.plink.map
  gunzip CIMBA_imputed_chr${chr}.plink.dosage.gz
  ./plink2 --import-dosage CIMBA_imputed_chr${chr}.plink.dosage format=1 --map CIMBA_imputed_chr${chr}_v2.plink.map --fam CIMBA_imputed_chr${chr}.plink.fam --make-bed --out CIMBA_chr${chr}.Imputed &
done


for chr in {2..22} 
do
	echo "CIMBA_chr${chr}.Imputed.bed CIMBA_chr${chr}.Imputed.bim CIMBA_chr${chr}.Imputed.fam" >> unir_cromosomas.txt
done

plink --noweb --make-bed --bfile CIMBA_chr1.Imputed --merge-list unir_cromosomas.txt --out CIMBA_Imputed
~~~
## **Añadir info caso-control**
~~~
grep -v "#" phs001321.v1.pht006198.v1.p1.c1.OncoArray_CIMBA_Sample_Attributes.GRU.txt| grep -v "ANALYTE_TYPE"| cut -d$'\t' -f 2,6| tail -n +2 > phenotypes_CIMBA.txt

cat phenotypes_CIMBA.txt | awk '{print $1,$0}' | grep "Y" > cases_pheno_CIMBA.txt
sed -i 's/Y//g' cases_pheno_CIMBA.txt 

plink --bfile CIMBA_Imputed --make-pheno cases_pheno_CIMBA.txt  '*' --make-bed --out CIMBA_pheno_added

grep -E "1$" phs001321.v1.pht006197.v1.p1.c1.OncoArray_CIMBA_Subject_Phenotypes.GRU.txt| grep -v "#" | awk '{print $2,$2}' > non_CEU.txt
~~~

## **Añadir info de sexo de individuos**
~~~
cat CIMBA_pheno_added.fam | awk '{print $1,$2, "2"}' > sex_CIMBA.txt

plink --bfile CIMBA_pheno_added --update-sex sex_CIMBA.txt --remove non_CEU.txt --make-bed --out CIMBA_EUR_pheno_sex_added
~~~

## **Control de calidad Target Data** 
~~~
plink --bfile CIMBA_EUR_pheno_sex_added --maf 0.01 --hwe 1e-6 --geno 0.1 --write-snplist --out EUR_CIMBA.QC
plink --bfile CIMBA_EUR_pheno_sex_added --extract EUR_CIMBA.QC.snplist --mind 0.1 --make-just-fam --out EUR_CIMBA.QC

#plink --bfile CIMBA_EUR_pheno_sex_added --keep EUR_CIMBA.QC.fam --extract EUR_CIMBA.QC.snplist --indep-pairwise 200 50 0.25 --out EUR_CIMBA.QC

plink --bfile CIMBA_EUR_pheno_sex_added --keep EUR_CIMBA.QC.fam --extract EUR_CIMBA.QC.snplist --het --out EUR_CIMBA.QC

R
dat <- read.table("EUR_CIMBA.QC.het", header=T) # Read in the EUR.het file, specify it has header
m <- mean(dat$F) # Calculate the mean  
s <- sd(dat$F) # Calculate the SD
valid <- subset(dat, F <= m+3*s & F >= m-3*s) # Get any samples with F coefficient within 3 SD of the population mean
write.table(valid[,c(1,2)], "EUR_CIMBA.valid.sample", quote=F, row.names=F) # print FID and IID for valid samples
q() # exit R
~~~
#### **Merge variants**
~~~
R
# Read in bim file
bim <- read.table("CIMBA_EUR_pheno_sex_added.bim")
colnames(bim) <- c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2")
# Read in QCed SNPs
qc <- read.table("info_snp_tab", header = T, stringsAsFactors = F)
qc <- qc[c(1,2,3,4,5)]
info <- merge(bim, qc, by.x = c("CHR", "BP", "B.A1", "B.A2"), by.y = c("chr", "pos", "a0", "a1"))
match_variants <-info$SNP
write.table(match_variants, "EUR_CIMBA_GWAS.shared_variants", quote = F, row.names = F, col.names = F)
q() # exit R
~~~
##### **--EUR_CGEMS_GWAS.shared_variants**

~~~
plink --bfile CIMBA_EUR_pheno_sex_added --extract EUR_CIMBA_GWAS.shared_variants --keep EUR_CIMBA.valid.sample --rel-cutoff 0.125 --out EUR_CIMBA.QC
plink --bfile CIMBA_EUR_pheno_sex_added --make-bed --keep EUR_CIMBA.QC.rel.id --out EUR_CIMBA.QC --extract EUR_CIMBA_GWAS.shared_variants
~~~

Minimac3 --refHaps HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz --processReference --prefix HRC.r1-1.GRCh37.wgs.mac5.sites