#!/usr/bin/perl

$cmd = join ' ', @ARGV;

#java to use
$JAVA = "/home/uec-00/shared/production/software/java/default/bin/java -Xmx7G";
$libpath = "/home/uec-00/shared/production/software/uecgatk/2013-11-20/lib";
$uecgatk = "/home/uec-00/shared/production/software/uecgatk/2013-11-20/uecgatk.jar";

#dep libs, last one is our actual inhouse code.
$JARLIBS = "$libpath/GenomeAnalysisTK-2.7-4-g6f46d11.jar:$libpath/batik-anim.jar:$libpath/batik-css.jar:$libpath/batik-dom.jar:$libpath/batik-ext.jar:$libpath/batik-parser.jar:$libpath/batik-svg-dom.jar:$libpath/batik-util.jar:$libpath/batik-xml.jar:$libpath/charts4j-1.2.jar:$libpath/commons-math-2.2.jar:$libpath/heatMap.jar:$libpath/sam-1.61.jar:$libpath/tuple.jar:$libpath/UscKeck.jar:$libpath/xml-apis-ext.jar:$libpath/xml-apis.jar:$libpath/genomeLibs.jar:$uecgatk";


$execString = "$JAVA -classpath $JARLIBS org.broadinstitute.sting.gatk.CommandLineGATK $cmd ";
print "$execString\n";
system($execString);