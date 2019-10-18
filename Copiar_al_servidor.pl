#!/usr/bin/perl
print "Ingrese Archivo a copiar\n";
$file1=<>;
chomp ($file1);
print "Ingrese ruta y nombre de destino\n";
$file2=<>;
chomp ($file2); 
system ("scp -P 1313 $file1 adolforojas\@200.89.65.156:/home/adolforojas/$file2");
