
~~~
plink1.9 --keep-fam Denmark_samples --mind 0.01 --geno 0.01 --bfile dbgap_bcac_oncoarray_c8 --make-bed --out Drive_denmark

for chr in {1..22} 
do 
    plink1.9 --bfile Drive_denmark --chr $chr --make-bed --out Drive_denmark_chr_${chr} 
    plink1.9 --bfile Drive_denmark_chr_${chr} --list-duplicate-vars ids-only suppress-first 
    plink1.9 --bfile Drive_denmark_chr_${chr} -exclude plink.dupvar --make-bed --out Drive_denmark_chr_${chr}.DuplicatesRemoved 
done
~~~

# Phasing main process
~~~
for chr in {1..22} # Ejemplo Taruca (hecho en Taruca)
do
  echo "-B plink_files/Drive_denmark_chr_${chr}.DuplicatesRemoved -M ../1000GP_Phase3/genetic_map_chr${chr}_combined_b37.txt -O phased_data_last/Drive_denmark_chr_${chr}.phased --states 300 -T 10 --burn 10 --main 50 --prune 10 --output-log phased_data_last/Shapeit_Drive_denmark_chr_${chr}.phased" >> myCommands.txt
done

cat myCommands.txt | xargs -P6 -n18 shapeit & # P4 indica 4 trabajos a la vez n12 la cantidad de argumentos del comando
rm myCommands.txt
~~~
# Haplotype Harmonizer
~~~
for i in {1..22} 
do     
  echo "-jar GenotypeHarmonizer-1.4.23/GenotypeHarmonizer.jar --inputType SHAPEIT2 --input Denmark_data/phased_data_last/Drive_denmark_chr_$i.phased --forceChr $i --keep --output Denmark_data/strand_change_phased/Denmark_strand_change.chr$i --refType VCF --ref ALL.chr$i.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz --update-reference-allele --update-id" >> GenotypeHarmonizer_commands.txt
done
cat GenotypeHarmonizer_commands.txt | xargs -P10 -n17 java & # P4 indica 4 trabajos a la vez n12 la cantidad de argumentos del comando
rm GenotypeHarmonizer_commands.txt
~~~
# Imputation main process
~~~
python3
import pandas as pd
import os 
lengths_chr = pd.read_csv("../Chromosomes_length.csv")
os.system("rm IMPUTE2_commands.txt")
for i in range(1,23):
  print(i)
  Chromosome_length = int(lengths_chr.loc[lengths_chr.Chromosome == str(i),"Total length (bp)"])
  print(Chromosome_length)
  chunk = 1
  end_site = 0
  while end_site < Chromosome_length:
    end_site+=5000000
    if end_site > Chromosome_length:
      end_site = Chromosome_length
    os.system("echo -use_prephased_g -known_haps_g strand_change_phased/Denmark_strand_change.chr"+str(i)+".haps -filt_rules_l 'EUR<0.001' -h ../1000GP_Phase3/1000GP_Phase3_chr"+str(i)+".hap.gz -int "+str((chunk-1)*5000000+1)+" "+str(end_site)+" -l ../1000GP_Phase3/1000GP_Phase3_chr"+str(i)+".legend.gz -m ../1000GP_Phase3/genetic_map_chr"+str(i)+"_combined_b37.txt -Ne 20000 -pgs_miss -buffer 500 -k_hap 800 -o Imputed_data/Drive_denmark_chr_"+str(i)+".chunk0"+str(chunk)+".impute2 >> IMPUTE2_commands.txt") ## Harmonizer applied
    chunk+=1

quit()

cat IMPUTE2_commands.txt |xargs -P60 -n23 ./impute2
rm IMPUTE2_commands.txt
~~~
## Impute CHR concat
~~~
python3
from os import walk
import os
import pandas as pd
Output_dir_files = pd.DataFrame(next(walk("./"), (None, None, []))[2], columns={"file"})
impute_files = Output_dir_files.loc[Output_dir_files.file.str.contains(".impute2$")].copy()

for cromosoma in range(1,23):
  print(cromosoma)
  chr_chunks =  impute_files.file.loc[Output_dir_files.file.str.contains("chr_" + str(cromosoma))].tolist()
  print(chr_chunks)
  comando = "cat"
  for i in range(1,60):
    chunk_file = "Drive_denmark_chr_"+ str(cromosoma) +".chunk0"+ str(i) +".impute2"
    if chunk_file in chr_chunks:
      comando+= " " + chunk_file
  comando += " > joined_chr/Drive_denmark_chr_"+ str(cromosoma) +".complete.impute2"
  print(comando)
  os.system(comando)


for i in range(1,23):
    os.system("sed -i 's/---/"+str(i)+"/g' joined_chr/Drive_denmark_chr_"+str(i)+".complete.impute2")

quit()
~~~
## Impute info concat
~~~
python3
from os import walk
import os
import pandas as pd
Output_dir_files = pd.DataFrame(next(walk("./"), (None, None, []))[2], columns={"file"})
impute_files = Output_dir_files.loc[Output_dir_files.file.str.contains(".impute2_info$")].copy()

for cromosoma in range(1,23):
  print(cromosoma)
  chr_chunks =  impute_files.file.loc[Output_dir_files.file.str.contains("chr_" + str(cromosoma))].tolist()
  print(chr_chunks)
  comando = "cat"
  for i in range(1,60):
    chunk_file = "Drive_denmark_chr_"+ str(cromosoma) +".chunk0"+ str(i) +".impute2_info"
    if chunk_file in chr_chunks:
      comando+= " " + chunk_file
  comando += " > joined_chr/Drive_denmark_chr_"+ str(cromosoma) +".complete.impute2_info"
  print(comando)
  os.system(comando)

for i in range(1,23):
    os.system("sed -i 's/---/"+str(i)+"/g' joined_chr/Drive_denmark_chr_"+str(i)+".complete.impute2_info")

quit()
~~~    

cat *.impute2_info |awk '!x[$0]++'|awk -F " " '{ if(($7 >= 0.8)) { print $2} }'| tail -n +2 > pass_INFO_filter.txt
cd joined_chr
~~~
python3
import os
os.system("mkdir -p plink_data")
for i in range(1,23):
  os.system("mv Drive_denmark_chr_"+str(i)+".complete.impute2 Drive_denmark_chr_"+str(i)+".complete.impute2.gen")
  os.system("./gtool -G --g Drive_denmark_chr_"+str(i)+".complete.impute2.gen --s Drive_Denmark.phased.sample --ped plink_data/Denmark_chr_"+str(i)+".ped --map plink_data/Denmark_chr_"+str(i)+".map") #--inclusion pass_INFO_filter.txt
  os.system("plink --file plink_data/Denmark_chr_"+str(i)+" --make-bed --out plink_data/DRIVE_Denmark_chr_"+str(i)+"_plink &")

quit()
~~~
~~~
cd plink_data
for chr in {2..22} 
do
	echo "DRIVE_Denmark_chr_${chr}_plink.bed DRIVE_Denmark_chr_${chr}_plink.bim DRIVE_Denmark_chr_${chr}_plink.fam" >> unir_cromosomas.txt
done

plink --noweb --make-bed --bfile DRIVE_Denmark_chr_1_plink --merge-list unir_cromosomas.txt --out Drive_Denmark_Imputed
~~~
awk '{print $2,$1}' cases_Drive_Denmark > cases_Drive_Denmark2
~~~

plink --bfile Drive_Denmark_Imputed --make-pheno cases_Drive_Denmark2 '*' --make-bed --out Drive_Denmark_pheno_added

plink --bfile Drive_Denmark_pheno_added --geno 0.02 --write-snplist --out Demmark_DRIVE.QC #--hwe 1e-6

R
# Read in bim file
bim <- read.table("Drive_Denmark_pheno_added.bim")
colnames(bim) <- c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2")  #608479 sin filtro
good_snp <- read.delim("Demmark_DRIVE.QC.snplist", header=FALSE)
bim <- bim[bim$SNP %in% good_snp$V1,]
# Read in QCed SNPs
qc <- read.table("info_snp_tab", header = T, stringsAsFactors = F)
qc <- qc[c(1,2,3,4,5)]
bim$B.A2 <- sapply(strsplit(as.character(bim$SNP),':'), "[", 3)
bim$B.A1 <- sapply(strsplit(as.character(bim$SNP),':'), "[", 4)
info <- merge(bim, qc, by.x = c("CHR", "BP", "B.A1", "B.A2"), by.y = c("chr", "pos", "a0", "a1"))
match_variants <-info$SNP
write.table(match_variants, "Drive_Denmark.shared_variants", quote = F, row.names = F, col.names = F)
q() # exit R
n
wc -l Drive_Denmark.shared_variants

plink --bfile Drive_Denmark_pheno_added --extract Drive_Denmark.shared_variants --mind 0.01 --make-just-fam --out Demmark_DRIVE.QC #--extract Demmark_DRIVE.QC.snplist

plink --bfile Drive_Denmark_pheno_added --keep Demmark_DRIVE.QC.fam --extract Drive_Denmark.shared_variants --het --out Demmark_DRIVE.QC

R
dat <- read.table("Demmark_DRIVE.QC.het", header=T) # Read in the EUR.het file, specify it has header
m <- mean(dat$F) # Calculate the mean  
s <- sd(dat$F) # Calculate the SD
valid <- subset(dat, F <= m+3*s & F >= m-3*s) # Get any samples with F coefficient within 3 SD of the population mean
write.table(valid[,c(1,2)], "Demmark_DRIVE.valid.sample", quote=F, row.names=F) # print FID and IID for valid samples
q() # exit R
n


plink --bfile Drive_Denmark_pheno_added --make-bed --keep Demmark_DRIVE.valid.sample --out EUR_Drive_Denmark.QC --extract Drive_Denmark.shared_variants

R
bim <- read.table("EUR_Drive_Denmark.QC.bim")
colnames(bim) <- c("CHR", "SNP", "CM", "BP", "ALT_old", "REF_old")  #608479 sin filtro
# Read in QCed SNPs
bim$ALT_New <- sapply(strsplit(as.character(bim$SNP),':'), "[", 4)
bim$REF_New <- sapply(strsplit(as.character(bim$SNP),':'), "[", 3)
bim$CHR <- NULL
bim$CM <- NULL
bim$BP <- NULL
write.table(bim, "reemplazo_alelos.txt", quote = F, row.names = F, col.names = F)
q()
n

plink --bfile EUR_Drive_Denmark.QC --update-alleles  reemplazo_alelos.txt --make-bed --out EUR_Drive_Denmark.QC2

R
load("../../../../2do_Objetivo/solo_pareados/2_entorno_validacion_cruzada.RData")
library(bigsnpr)
#Sys.setenv(LANG="en_US.UTF-8")
system("rm *.QC.rds")
system("rm *.QC.bk")
snp_readBed("EUR_Drive_Denmark.QC2.bed") # Read from bed/bim/fam, it generates .bk and .rds files.
obj.bigSNP2 <- snp_attach("EUR_Drive_Denmark.QC2.rds") # Attach the "bigSNP" object in R session
# Get aliases for useful slots
G2   <- obj.bigSNP2$genotypes
CHR2 <- obj.bigSNP2$map$chromosome
POS2 <- obj.bigSNP2$map$physical.pos
y2   <- obj.bigSNP2$fam$affection - 1
NCORES <- 20 # nb_cores()

pred_inf <- big_prodVec(G2, beta_inf, ind.row = ind.test, ind.col = info_snp$`_NUM_ID_`, ncores = NCORES)
AUCBoot(pred_inf, y[ind.test])