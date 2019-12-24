#!/usr/bin/perl

qx/screen -d -R GWAS/;

sub Directorios {

system ("mkdir -p ../Temp");
system ("mkdir -p ../IDs");
system ("mkdir -p ../Resultados/GWAS");
system ("rm ../Temp/*");
}

Directorios();

print "Nombre del archivo en que trabajara\n\n";
$Name = <>;
chomp ($Name);
print "Numero de columna que trabajara\n\n";
$Col = <>;
chomp ($Col);

sub GO_Validados {
open (ENTRADA1, ">../Temp/$Name.in");
qx/chmod a+r $Name/;
$Entrada= qx/grep $Name | cut -d$\'\t' -f $Col | grep -Eo "rs[0-9]+" | sort | uniq/;
print ENTRADA1 "$Entrada";
close(ENTRADA1);

open (VALIDADOS, "../Temp/$Name.in");
open (SALIDA1, ">../Resultados/GWAS/$Name.tsv");
print SALIDA1 "CancerRelated?\tGeneID\tGeneName\tPathways->\n";
while ($RS_ID=<VALIDADOS>) {
    chomp ($RS_ID);
    $Data1 = qx/esearch -db snp -query $RS_ID | efetch -format docsum | xtract -pattern DocumentSummary -element GENE_ID NAME | uniq/;
    chomp ($Data1);
    @Gene_ID= split(/\s+/,$Data1);
    $KEGG_data_GO = qx/curl http:\/\/rest.kegg.jp\/get\/hsa:$Gene_ID[0] | grep -Eo "hsa[0-9]+[^\]]+"/;
    $KEGG_data_GO=~ s/\n/\t/ig;
    $KEGG_data_GO=~ s/hsa[0-9]+\t//ig;
    if ($KEGG_data_GO=~ m/hsa052[0-9][0-9]/ig){
        $Cancer= "Yes";
    }else{
        if ($KEGG_data_GO=~ m/hsa\w+/ig){
            $Cancer= "No";
        } else {
            $Cancer= "NA";
        }
    }
    print SALIDA1 "$Cancer\t$Gene_ID[0]\t$Gene_ID[1]\t$KEGG_data_GO\n";
    $Cancer="";
    }
close (VALIDADOS);