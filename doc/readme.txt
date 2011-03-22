to run your own walkers with gatk:
make sure the gatk.jar is in your classpath and it will detect your walker at runtime. you must always call gatk's main().
you can call your walker by name as follows (I have a walker named "ZackTest", src is ZackTestWalker.java):

java -classpath GenomeAnalysisTK.jar org.broadinstitute.sting.gatk.CommandLineGATK -T ZackTest -R exampleFASTA.fasta -I exampleBAM.bam -o output.txt


for help, running:

java -classpath GenomeAnalysisTK.jar org.broadinstitute.sting.gatk.CommandLineGATK

with no options should show you all available walkers, including the ones in uecgatk 