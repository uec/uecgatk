#!/usr/bin/perl

$cmd = join ' ', @ARGV;

#java to use
$JAVA = "/home/uec-00/shared/production/software/java/default/bin/java -Xmx15G";

#dep libs, last one is our actual inhouse code.
$JARLIBS = "/home/uec-00/shared/production/software/uecgatk/default/lib/GenomeAnalysisTK.jar:/home/uec-00/shared/production/software/uecgatk/default/lib/AnalyzeCovariates.jar:/home/uec-00/shared/production/software/uecgatk/default/lib/StingUtils.jar:/home/uec-00/shared/production/software/uecgatk/default/lib/batik-anim.jar:/home/uec-00/shared/production/software/uecgatk/default/lib/batik-css.jar:/home/uec-00/shared/production/software/uecgatk/default/lib/batik-dom.jar:/home/uec-00/shared/production/software/uecgatk/default/lib/batik-ext.jar:/home/uec-00/shared/production/software/uecgatk/default/lib/batik-parser.jar:/home/uec-00/shared/production/software/uecgatk/default/lib/batik-svg-dom.jar:/home/uec-00/shared/production/software/uecgatk/default/lib/batik-util.jar:/home/uec-00/shared/production/software/uecgatk/default/lib/batik-xml.jar:/home/uec-00/shared/production/software/uecgatk/default/lib/charts4j-1.2.jar:/home/uec-00/shared/production/software/uecgatk/default/lib/commons-math-1.1.jar:/home/uec-00/shared/production/software/uecgatk/default/lib/heatMap.jar:/home/uec-00/shared/production/software/uecgatk/default/lib/sam-1.26.jar:/home/uec-00/shared/production/software/uecgatk/default/lib/tuple.jar:/home/uec-00/shared/production/software/uecgatk/default/lib/UscKeck.jar:/home/uec-00/shared/production/software/uecgatk/default/lib/xml-apis-ext.jar:/home/uec-00/shared/production/software/uecgatk/default/lib/xml-apis.jar:/home/uec-00/shared/production/software/uecgatk/default/lib/genomeLibs.jar:/home/uec-00/shared/production/software/uecgatk/default/uecgatk.jar";



$execString = "$JAVA -classpath $JARLIBS org.broadinstitute.sting.gatk.CommandLineGATK $cmd ";
print "$execString\n";
system($execString);
