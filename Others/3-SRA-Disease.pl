#!/usr/bin/perl
system ("mkdir -p Temp");
system ("mkdir -p IDs");
system ("mkdir -p Resultados");
system ("rm Temp/*");

print "Ingrese enfermedad a buscar\n";
$disease=<>;
print "Pasos N°1 y N°2, Busqueda y obtencion de IDs SRA\n\n";
	chomp ($disease);
open (disease, ">Temp\/$disease") or die "Error en escritura de IDs (diseasees)";
## Obtencion de IDs a partir del input
$datos = qx/esearch -db sra -query '$disease AND Human [ORGN] AND public' | efetch -format docsum | xtract -pattern DocumentSummary -ACC \@acc -block DocumentSummary -element "&ACC"/;
## Generar formato de lista para el Output
$datos=~ s/\t/\n/g;
## Filtrar IDs que no sean del tipo SRX o ERX
$datos=~ s/[SED]R[PZRAS]\d+\n//gi;
print "Procesando $disease\t ";
## Guardar IDs a un archivo correspondiente al Input
print disease "$datos\n";
print " Finalizado\n\n";
## Restaurar valor nulo de la Variaable en caso de no obtener resultados con el Input
$datos="";
close (disease);

## Eliminar IDs repetidas y generar archivo unico
system ("cat Temp\/\*|sort|uniq>Temp\/all.txt");

print "Recopilando IDs antiguas\n\n";

system ("cat Temp\/all.txt IDs\/ID_disease | sort | uniq -u > IDs\/IDs_Nuevos");
$IDs_Nuevos = qx/wc -l IDs\/IDs_Nuevos/;
chomp ($IDs_Nuevos);
print "$IDs_Nuevos\n\n";

open (SALIDA, ">Temp\/Salida.txt");
## Obtener infromacion de las IDs recopiladas
$muestra= qx/epost -db sra -input IDs\/IDs_Nuevos -format acc | esummary -format runinfo -mode xml/;
print SALIDA "$muestra\n";
$muestra="";
close (SALIDA);

print "Paso N°3 Filtrado preeliminar\n\n";

open (INPUT, "Temp\/Salida.txt") or die "in";
open (OUT, ">Temp\/Filtrados.txt") or die "Error Filtrado";
## Cambio de delimitador segun el formato de archivo 
$/='<Row>';
while ($datos=<INPUT>){
	chomp ($datos);
	if ($datos =~ m/\<ScientificName\>Homo sapiens\<\/ScientificName\>/i){; # Filtrar las que no sean de Homo sapiens
		if ($datos =~ m/biosample/i) {
				if ($datos =~ m/(metagenome)|(metagenomic)|(microbiota)/i) { # filtrar Otros tipos de estudio no pertinentes
				} else {if ($datos=~ m/\<Consent\>public\<\/Consent\>/i) { # Filtrar para que los datos sean publicos
				if ($datos=~ m/\<LibraryStrategy\>RNA-Seq\<\/LibraryStrategy\>/i) {
				print OUT "$datos\n";
				
}}}}}}
close (INPUT);
close (OUT);

print "Paso N°4 Extraccion de datos y metadatos; generacion de tabla\n\n";

open (TABLA, "Temp\/Filtrados.txt") or die "Error Tabla";
open (SALIDA, ">>Muestras.tsv");
open (UTILIZADOS, ">>IDs\/ID_disease");
$/="</Row>";
## Imprimir nombres de las Columnas
print SALIDA "1_Run	ReleaseDate	LoadDate	AssemblyName	Spots	Bases	Spots with mates	avgLength	Size(MB)	Dowload	Experiment	Library	Strategy	Selection	Source	Layout	InsertSize	InsertDev	Platform	Model	SRAstudy	Bioproject	ProjectID	Sample	BioSample	SampleType	TaxID	ScientificName	SampleName	Sex	Tumor	Center	Submission	Acceso	RunHash	ReadHash	Project Description	Sample Metadata-->\n";
	while ($objetos= <TABLA>){;
		chomp ($objetos);
	if ($objetos =~ /<Experiment>\s*["-]?(\w+\s*\w*)"?/i) {
				$Experiment = $1; } else {
					$Experiment = "";}
# #Cada de los siguientes if genera la extraccion de informaccion a partir de el archivo de cada ID (linea 26 Aprox)
	if ($objetos =~ /<BioSample>\s*["-]?(\w+\s*\w*)"?/i) { 
			 $BioSample = $1; } else {
				$BioSample = "";}
## Obtencion de la informacion de la muestra relacionada al ID evaluado
	$metadatos= qx/esearch -db biosample -query $BioSample | efetch -format docsum | xtract -pattern DocumentSummary -block Attribute -element Attribute/;
	$metadatos=~ s/\n//g;
	if ($metadatos =~ m/(metagenome)|(metagenomic)|(microbiota)/i) { # Segundo filtrado en caso de falla del anterior
	} else { # En caso de no presentar palabras del filtro se continua extrayendo mas informacion
	if ($objetos =~ /<Run>\d*\s*["-]?(\w+\s*\w*)"?/i) {
	      $Run = $1; } else {
 			   $Run = "";}
	if ($objetos =~ /<ReleaseDate>\s*["-]?(\w+-\w+-\w+\s*\d+:\d+:\d+)"?/i) {
        $ReleaseDate = $1; } else {
					$ReleaseDate = "";}
    if ($objetos =~ /<LoadDate>\s*["-]?(\w+-\w+-\w+\s*\d+:\d+:\d+)"?/i) {
			  $LoadDate = $1; } else {
					$LoadDate = "";}
	if ($objetos =~ /<AssemblyName>\s*["-]?(\w+\.*\w*)"?/i) {
				$AssemblyName = $1; } else {
					$AssemblyName = "";}
	if ($objetos =~ /<spots>\s*["-]?(\d+)"?/i) {
	      $spots = $1; } else {
					$spots = "";}
	if ($objetos =~ /<bases>\s*["-]?(\d+)"?/i) {
	      $bases = $1; } else {
					$bases = "";}
	if ($objetos =~ /<spots_with_mates>\s*["-]?(\d+)"?/i) {
	      $spots_with_mates = $1; } else {
					$spots_with_mates = "";}
	if ($objetos =~ /<avgLength>\s*["-]?(\d+)"?/i) {
	      $avgLength = $1; } else {
					$avgLength = "";}
	if ($objetos =~ /<size_MB>\s*["-]?(\d+)"?/i) {
	      $size_MB = $1; } else {
					$size_MB = "";}
	if ($objetos =~ /<download_path>\s*["-]?(\w+\s*\w*\S*\d+)"?/i) {
			  $download_path = $1; } else {
					$download_path = "";}
	if ($objetos =~ /<LibraryName>\s*["-]?(\w+\.*\w*)"?/i) {
				$LibraryName = $1; } else {
					$LibraryName = "";}
	if ($objetos =~ /<LibraryStrategy>\s*["-]?(\w+\s*\w*)"?/i) {
				$LibraryStrategy= $1; } else {
					$LibraryStrategy = "";}
	if ($objetos =~ /<LibrarySelection>\s*["-]?(\w+\s*\w*)"?/i) {
				$LibrarySelection = $1; } else {
					$LibrarySelection = "";}
	if ($objetos =~ /<LibrarySource>\s*["-]?(\w+\s*\w*)"?/i) {
				$LibrarySource = $1; } else {
					$LibrarySource = "";}
	if ($objetos =~ /<LibraryLayout>\s*["-]?(\w+\s*\w*)"?/i) {
				$LibraryLayout = $1; } else {
					$LibraryLayout = "";}
	if ($objetos =~ /<InsertSize>\s*["-]?(\w+\s*\w*)"?/i) {
				$InsertSize = $1; } else {
					$InsertSize = "";}
	if ($objetos =~ /<InsertDev>\s*["-]?(\w+\s*\w*)"?/i) {
				$InsertDev = $1; } else {
					$InsertDev = "";}
	if ($objetos =~ /<Platform>\s*["-]?(\w+\s*\w*)"?/i) {
				$Platform = $1; } else {
					$Platform = "";}
	if ($objetos =~ /<Model>\s*["-]?(\w+\s*\w*)"?/i) {
				$Model = $1; } else {
					$Model = "";}
	if ($objetos =~ /<SRAStudy>\s*["-]?(\w+\s*\w*)"?/i) {
				$SRAStudy = $1; } else {
					$SRAStudy = "";}
	if ($objetos =~ /<BioProject>\s*["-]?(\w+\s*\w*)"?/i) {
				$BioProject = $1; } else {
					$BioProject = "";}
	## Obtencion de la informacion del proyecto
	$Project_description= qx/esearch -db bioproject -query $BioProject | efetch -format docsum | xtract -pattern DocumentSummary -block Project_Description -element Project_Description/;
	$Project_description=~ s/\n//g;
	if ($objetos =~ /<ProjectID>\s*["-]?(\w+\s*\w*)"?/i) {
				$ProjectID = $1; } else {
					$ProjectID = "";}
	if ($objetos =~ /<Sample>\s*["-]?(\w+\s*\w*)"?/i) {
				$Sample = $1; } else {
					$Sample = "";}
	if ($objetos =~ /<SampleType>\s*["-]?(\w+\s*\w*)"?/i) {
				$SampleType = $1; } else {
					$SampleType = "";}
	if ($objetos =~ /<TaxID>\s*["-]?(\w+\s*\w*)"?/i) {
				$TaxID = $1; } else {
					$TaxID = "";}
	if ($objetos =~ /<ScientificName>\s*["-]?(\w+\s*\w*)"?/i) {
				$ScientificName = $1; } else {
					$ScientificName = "";}
	if ($objetos =~ /<SampleName>\s*["-]?(\w+\s*\w*)"?/i) {
		  $SampleName = $1; } else {
						$SampleName = "";}
	if ($objetos =~ /<Sex>\s*["-]?(\w+\s*\w*)"?/i) {
	      $Sex = $1; } else {
					$Sex = "";}
	if ($objetos =~ /<Tumor>\s*["-]?(\w+\s*\w*)"?/i) {
		  	$Tumor = $1; } else {
					$Tumor = "";}
	if ($objetos =~ /<CenterName>\s*["-]?(\w+\s*\w*)"?/i) {
	      $CenterName = $1; } else {
					$CenterName = "";}
	if ($objetos =~ /<Submission>\s*["-]?(\w+\s*\w*)"?/i) {
				$Submission = $1; } else {
					$Submission = "";}
	if ($objetos =~ /<Consent>\s*["-]?(\w+\s*\w*)"?/i) {
				$Consent = $1; } else {
					$Consent = "";}
	if ($objetos =~ /<RunHash>\s*["-]?(\w+\s*\w*)"?/i) {
				$RunHash = $1; } else {
					$RunHash = "";}
	if ($objetos =~ /<ReadHash>\s*["-]?(\w+\s*\w*)"?/i) {
				$ReadHash = $1; } else {
					$ReadHash = "";}
print SALIDA "$Run	$ReleaseDate	$LoadDate	$AssemblyName	$spots	$bases	$spots_with_mates	$avgLength	$size_MB	$download_path	$Experiment	$LibraryName	$LibraryStrategy	$LibrarySelection	$LibrarySource	$LibraryLayout	$InsertSize	$InsertDev	$Platform	$Model	$SRAStudy	$BioProject	$ProjectID	$Sample	$BioSample	$SampleType	$TaxID	$ScientificName	$SampleName	$Sex	$Tumor	$CenterName	$Submission	$Consent	$RunHash	$ReadHash	$Project_description	$metadatos\n";
$metadatos="";
$Project_description="";
	}
	print UTILIZADOS "$Experiment\n";
	}
close (TABLA);
close (SALIDA);
close (UTILIZADOS);
system ("sort Muestras.tsv | uniq > Resultados/Tabla_Muestras.tsv");
system ("wc -l Resultados/Tabla_Muestras.tsv");
system ("less -S Resultados/Tabla_Muestras.tsv");