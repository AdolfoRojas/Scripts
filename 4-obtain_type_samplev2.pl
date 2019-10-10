#!/usr/bin/perl

$input_doc=$ARGV[0];
chomp ($input_doc);
open(DOC, "$input_doc");
open (FILTERED, ">Filtered_$input_doc");
open(SALIDA, ">New_$input_doc");
while($DATA=<DOC>){
    chomp ($DATA);
if ($DATA=~ m/(Metadata-->)/i){
    $Type_sample="Clasificacion";
}else{if ($DATA=~ m/(mcf7)|(mcf-7)/i){
    $Type_sample="Breast cancer cell line";
}else{if ($DATA=~ m/(Xenograft)|(PRJNA310000)|(PDX)/i){
    $Type_sample="Xenograft tumor";
}else{if ($DATA=~ m/(MDA\S*MB\S*231)|(MB\S*231)/i){
    $Type_sample="Breast cancer cell line";
} else{if ($DATA=~ m/(SKBR3)/i){
    $Type_sample="Breast cancer cell line";
}else{if ($DATA=~ m/(PRJNA378161)|(PRJNA454715)|(PRJNA272156)|(PRJNA142565)|(primary tumour)|(fresh)|(PRJNA540861)|(PRJNA531611)|(PRJNA229096)|(PRJNA170122)|(PRJNA314351)|(PRJNA292118)|(PRJEB13586)|(PRJNA454288)|(PRJNA305054)|(PRJNA327871)|(PRJNA193673)|(PRJEB15096)|(Solid tumor)|(FFPE)|(PRJNA172761)|(PRJNA484546)|(PRJNA497449)|(breast cancer tissue)|(primary tumor)|(breast tumor tissue)|(Primary Breast)|(patient)|(tumor sample)|(biopsy)|(PRJNA317919)|(primary breast tumor)/i){
    $Type_Tumoral="Tumoral Sample";##############
}else{if ($DATA=~ m/(MCF10CA)/i){
    $Type_sample="Breast cancer cell line";
}else{if ($DATA=~ m/(CTC)|(circulating tumor cell from blood)|(PRJNA225522)|(PRJNA304018)|(circulating tumor cell)|(PRJNA343124)/){
    $Type_sample="Circulating Tumor cells";
}else{if ($DATA=~ m/(mixture)/i){
    $Type_sample="Breast cancer cell line";
}else{if ($DATA=~ m/(MCF10a)/i){
    $Type_sample="Breast cancer cell line";
}else{if ($DATA=~ m/(miRNA)/i){
    $Type_sample="miRNA";
}else{if ($DATA=~ m/(T-*47D)/i){
    $Type_sample="Breast cancer cell line";
}else{if ($DATA=~ m/(ZR-*75-*\w*)/i){
    $Type_sample="Breast cancer cell line";
}else{if ($DATA=~ m/(breast cancer cell line)|(PRJNA170832)|(PRJNA471825)|(PRJNA450921)|(PRJNA417666)|(PRJNA371468)|(PRJNA276634)|(bt474)|(SRP110816)|(PRJNA388881)|(PRJNA225955)|(PRJNA247397)|(PRJNA247397)|(PRJNA379957)|(PRJNA354957)|(PRJNA188684)|(SRP133085)|(PRJNA427773)|(PRJNA515297)|(PRJNA515910)|(PRJNA384523)|(PRJNA231850)|(PRJNA263155)|(PRJNA142887)|(PRJNA543438)|(PRJNA358776)|(PRJNA338851)|(PRJEB25825)|(PRJNA423078)|(PRJNA432071)|(PRJNA447949)|(PRJNA445823)|(PRJNA488269)|(PRJNA414329)|(Immortalized)|(PRJNA318813)|(PRJEB30617)|(PRJEB29886)|(PRJEB29951)|(1833)|(line)|(BT474 cells)|(PRJNA401827)|(CAL51)|(PRJNA428150)|(PRJNA471883)|(PRJNA335350)|(PRJNA322427)|(PRJNA319220)|(PRJNA318513)/i){
    $Type_sample="Breast cancer cell line";
}else{if ($DATA=~ m/(Healthy)|(primary normal)/i){
    $Type_sample="Healthy control";
}else{if ($DATA=~ m/(MDA-MB-*\w*)/i){
    $Type_sample="Breast cancer cell line";
}else{if ($DATA=~ m/(SUM*\d*)/i){
    $Type_sample="Breast cancer cell line";
}else{if ($DATA=~ m/(breast cancer serum)/i){
    $Type_sample="Breast cancer serum";
}else{if ($DATA=~ m/(TNBC)|(PRJNA396019)|(PRJNA485429)|(Triple negative breast cancer)|(PRJEB33799)|(PRJNA172756)|(PRJNA545898)|(PRJNA396880)/i){
    $Type_Tumoral="Tumoral Sample";###################
}else{if ($DATA=~ m/(Primary lung cancer tumor)|(PRJNA421096)|(reduction mammoplasty)|(PRJNA181310)|(PRJNA318622)|(PRJNA147521)/i){
    $Type_sample="Other Tissues related to breast cancer";
}else{if ($DATA=~ m/(Regulatory T cells)|(350777)|(PRJNA483877)|(PRJNA380940)|(PRJNA345006)|(PRJNA302603)/i){
    $Type_sample="Study of Immunity cells";
}else{if ($DATA=~ m/(PRJNA515531)|(PRJNA488462)|(PRJNA497539)|(PRJNA527860)|(PRJEB26804)|(PRJNA306880)|(PRJEB27891)/i){
    $Type_sample="Treatments";
}else{
    $Type_sample="Sin Clasificar";
}}}}}}}}}}}}}}}}}}}}}}
print SALIDA "$Type_sample\t$DATA\n";
$Type_sample="";
if ($Type_Tumoral=""){}else{
print FILTERED "$Type_Tumoral\t$DATA\n";
}
$Type_Tumoral="";
}
close (DOC);
close(SALIDA);
close (FILTERED);

