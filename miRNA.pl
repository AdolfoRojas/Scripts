#!/usr/bin/perl

print "Nombre del miRNA que trabajara\n\n";
$Name = <>;
chomp ($Name);

sub GO_Validados {
open (ENTRADA1, ">Validados.in");
qx/chmod a+r *csv/;
$Entrada= qx/grep Tarbase *.csv | cut -d , -f 6 |sort | uniq | grep -Eo "[A-Z]+[0-9]+"/;
print ENTRADA1 "$Entrada";
close(ENTRADA1);

open (VALIDADOS, "Validados.in");
open (SALIDA1, ">$Name\_validated_targets.tsv");
print SALIDA1 "CancerRelated?\tGeneID\tGeneName\tPathways->\n";
while ($ID_Ensembl=<VALIDADOS>) {
    chomp ($ID_Ensembl);
    $Data1 = qx/esearch -db gene -query $ID_Ensembl | efetch -format docsum | xtract -pattern DocumentSummary -element Name Id/;
    chomp ($Data1);
    @Gene_ID= split(/\s+/,$Data1);
    $KEGG_data_GO = qx/curl http:\/\/rest.kegg.jp\/get\/hsa:$Gene_ID[1] | grep -Eo "hsa[0-9]+[^\]]+"/;
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
    print SALIDA1 "$Cancer\t$Gene_ID[1]\t$Gene_ID[0]\t$KEGG_data_GO\n";
    $Cancer="";
    }
close (VALIDADOS);
close (SALIDA1);
}
GO_Validados();

sub GO_Predichos {
open (ENTRADA2, ">Predichos.in");
qx/chmod a+r *xls/;
$Entrada= qx/cut -d$\'\t' -f 2,4,5 *xls | grep -v Gene/;
print ENTRADA2 "$Entrada";
close(ENTRADA2);

open (PREDICHOS, "Predichos.in");
open (SALIDA2, ">$Name\_predicted_targets.tsv");
print SALIDA2 "CancerRelated?\tGeneID\tGeneName\tPredScore\tPathways->\n";
while ($IDs=<PREDICHOS>) {
    chomp ($IDs);
    @Gene_ID2= split(/\s+/,$IDs);
    $KEGG_data_GO2 = qx/curl http:\/\/rest.kegg.jp\/get\/hsa:$Gene_ID2[1] | grep -Eo "hsa[0-9]+[^\]]+"/;
    $KEGG_data_GO2=~ s/\n/\t/ig;
    $KEGG_data_GO2=~ s/hsa[0-9]+\t//ig;
    if ($KEGG_data_GO2=~ m/hsa052[0-9][0-9]/ig){
        $Cancer= "Yes";
    }else{
        if ($KEGG_data_GO2=~ m/hsa\w+/ig){
            $Cancer= "No";
        } else {
            $Cancer= "NA";
        }
    }
    print SALIDA2 "$Cancer\t$Gene_ID2[1]\t$Gene_ID2[2]\t$Gene_ID2[0]\t$KEGG_data_GO2\n";
    $Cancer="";
    }
close (PREDICHOS);
close (SALIDA2);
}

GO_Predichos();