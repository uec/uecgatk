#!/usr/bin/env bash
#
# Run following command ./README.sh


# Make sure data files are in current directory
export GENOMEFA=./hg18_unmasked.plusContam.fa; # Must also have fa.fai index file
export BAM=./IMR90GNOMe-815C5ABXX-81B71ABXX-NODUPS.calmd.bam # Must also have bam.bami file

# Required Java libraries should all be here
export UECROOT=../..
export UECLIBS=$UECROOT/lib
export CLASSPATH=$UECROOT/uecgatk.jar:$UECLIBS/GenomeAnalysisTK.jar:$UECLIBS/tuple.jar:$UECLIBS/GenomeAnalysisTK.jar:$UECLIBS/genomeLibs.jar:$UECLIBS/StingUtils.jar:$UECLIBS/UscKeck.jar:$UECLIBS/sam-1.26.jar:$UECLIBS/charts4j-1.2.jar:$UECLIBS/biojava-live.jar:$UECLIBS/commons-math-1.1.jar;


# Documentation (very rudimentary) for walkers can will appear by running the walkers without arguments:
# java org.broadinstitute.sting.gatk.CommandLineGATK -T  GnomeSeqFeatureAlignment
# java org.broadinstitute.sting.gatk.CommandLineGATK -T  GnomeSeqAutocorrByRead

# Make a bed file that covers 2200 bp on either side of selected elements (input elements are ctcf.chr11.hg18.gff)
java -Xmx13000m  GtfFlankedWindows 2200 ./ctcf.chr11.hg18.gff ; grep -v 'track' ctcf.chr11.hg18.gff.2200bp.gtf  > ctcf.chr11.hg18.gff.2200bp.gtf.temp ; mv ctcf.chr11.hg18.gff.2200bp.gtf.temp ctcf.chr11.hg18.gff.2200bp.gtf ; ./gffToBed.pl ctcf.chr11.hg18.gff.2200bp.gtf ; 

# Output percent methylation across all reads (one line for each element)
java -Xmx13000m org.broadinstitute.sting.gatk.CommandLineGATK -T  GnomeSeqFeatureAlignment -R $GENOMEFA -I $BAM --min_mapping_quality_score 30  -cph  --minCT 3 --intervals ctcf.chr11.hg18.gff.2200bp.bed --elementGff ctcf.chr11.hg18.gff --windSize 2000 -et NO_ET -nt 1 --downscaleCols 1000 --outPrefix ALIGN_2000_1000_mapq30_minct3_mincontext0.90_ctcf --minContextFracReadsMatching 0.899  -pats WCG -pats GCH ;

# Output percent CG/GC pairs within read (SAMEREAD) and across all reads (DIFFERENTREAD)
# 00 = WCG unmethylated, GCH unmethylated
# 01 = WCG unmethylated, GCH methylated
# 10 = WCG methylated,   GCH unmethylated
# 11 = WCG methylated,   GCH methylated
java -Xmx13000m org.broadinstitute.sting.gatk.CommandLineGATK -T  GnomeSeqAutocorrByRead -R $GENOMEFA -I $BAM --min_mapping_quality_score 30  -cph   --intervals ctcf.chr11.hg18.gff.2200bp.bed --elementGff ctcf.chr11.hg18.gff --windSize 20 --featWindSize 2000 -et NO_ET -nt 1 --downscaleCols 1000 --outPrefix ALIGN_2000_1000_mapq30_minct3_mincontext0.90_ctcfPairs  ; 
# It outputs some extra ones we don't need
rm -f *ctcfPairs-WCG*


# Remove temp files
rm -f *2200bp*




# --- Unimportant notes ---
#$UECLIBS/biojava-live_1.6/apps-live.jar:$UECLIBS/biojava-live_1.6/bytecode.jar:$UECLIBS/biojava-live_1.6/commons-cli.jar:$UECLIBS/biojava-live_1.6/commons-collections-2.1.jar:$UECLIBS/biojava-live_1.6/commons-dbcp-1.1.jar:$UECLIBS/biojava-live_1.6/commons-pool-1.1.jar:$UECLIBS/biojava-live_1.6/demos-live.jar:$UECLIBS/biojava-live_1.6/jgrapht-jdk1.5.jar:$UECLIBS/biojava-live_1.6/junit-4.4.jar;
#/home/uec-00/shared/production/software/picard/default/picard-tools-1.40.zip ; 

