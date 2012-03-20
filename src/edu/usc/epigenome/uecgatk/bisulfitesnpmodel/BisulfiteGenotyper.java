package edu.usc.epigenome.uecgatk.bisulfitesnpmodel;


import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMSequenceDictionary;

import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.codecs.vcf.SortingVCFWriter;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.codecs.vcf.VCFFilterHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.DownsampleType;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.datasources.reads.SAMDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.BAQMode;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.Downsample;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyper;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyper.UGStatistics;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import edu.usc.epigenome.uecgatk.YapingWalker.verboseWriter;
import edu.usc.epigenome.uecgatk.bisulfitesnpmodel.BisulfiteSNPGenotypeLikelihoodsCalculationModel.MethylSNPModel;
import edu.usc.epigenome.uecgatk.bisulfitesnpmodel.BaseUtilsMore.*;
import edu.usc.epigenome.uecgatk.YapingWriter.FormatWriterBase;
import edu.usc.epigenome.uecgatk.YapingWriter.NOMeSeqReads;
import edu.usc.epigenome.uecgatk.YapingWriter.SortingCpgReadsWriter;
import edu.usc.epigenome.uecgatk.YapingWriter.SortingFormatWriterBase;
import edu.usc.epigenome.uecgatk.YapingWriter.cpgReads;
import edu.usc.epigenome.uecgatk.YapingWriter.cpgReadsWriterImp;
import edu.usc.epigenome.uecgatk.YapingWriter.SortingNOMeSeqReadsWriter;
import edu.usc.epigenome.uecgatk.YapingWriter.NOMeSeqReadsWriterImp;

/*
 * Bis-SNP/BisSNP: It is a genotyping and methylation calling in bisulfite treated 
 * massively parallel sequencing (Bisulfite-seq and NOMe-seq) on Illumina platform
 * Copyright (C) <2011>  <Yaping Liu: lyping1986@gmail.com>

 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
 * A Bisulfite genotyper. Works for single-sample data right now. 
 */

@BAQMode(QualityMode = BAQ.QualityMode.OVERWRITE_QUALS, ApplicationTime = BAQ.ApplicationTime.ON_INPUT)
@Reference(window=@Window(start=-500,stop=500))
@By(DataSource.REFERENCE)
@Downsample(by=DownsampleType.NONE)
public class BisulfiteGenotyper extends LocusWalker<BisulfiteVariantCallContext, BisulfiteGenotyper.ContextCondition> implements TreeReducible<BisulfiteGenotyper.ContextCondition> {


    @ArgumentCollection private static BisulfiteArgumentCollection BAC = new BisulfiteArgumentCollection();
 
    private static boolean autoEstimateC = false;
    private static boolean secondIteration = false;
    
    private static int MAXIMUM_CACHE_FOR_OUTPUT_VCF = 3000000;
    
    private static long COUNT_CACHE_FOR_OUTPUT_VCF = 0;
    
    private static long COUNT_CACHE_FOR_OUTPUT_READS = 0;
    
    private static int MAXIMUM_CACHE_FOR_OUTPUT_READS = 300000;
    
    protected TcgaVCFWriter writer = null;
    
    protected SortingTcgaVCFWriter multiThreadWriter = null;
    
    protected TcgaVCFWriter additionalWriterForDefaultTcgaMode = null;
    
    protected SortingTcgaVCFWriter multiAdditionalWriterForDefaultTcgaMode = null;
    
    protected FormatWriterBase readsWriter = null; //only works with single core right now
    
    protected SortingFormatWriterBase multiThreadCpgReadsWriter = null;
    
    protected TcgaVCFWriter verboseWriter = null;
    
    protected SAMFileWriter samWriter = null; //only works with single core right now
    
    private int SAMPLE_READS_MEAN_COVERAGE = 30;

    private BisulfiteGenotyperEngine BG_engine = null;
    
    private String argCommandline = "";
    
    //to record cytosine pattern methylation status estimated in the first iteration
    CytosineTypeStatus summary = null;

    /**
     * Inner class for collecting output statistics
     */
    public static class ContextCondition {
        /** The total number of passes examined -- i.e., the number of map calls */
        long nBasesVisited = 0;

        /** The number of bases that were potentially callable -- i.e., those not at excessive coverage or masked with N */
        long nBasesCallable = 0;

        /** The number of bases called confidently (according to user threshold), either ref or other */
        long nBasesCalledConfidently = 0;

        /** The number of bases for which calls were emitted */
        long nCallsMade = 0;

        /** The total number of extended events encountered */
        long nExtendedEvents = 0;
        
        /** The number of Cytosine bases called confidently (according to user threshold), either ref or other */
        long nCytosineBasesCalledConfidently = 0;
        
        /** The number of Cpg bases called confidently (according to user threshold), either ref or other */
        long nCpgBasesCalledConfidently = 0;
        
        /** The number of Chg bases called confidently (according to user threshold), either ref or other */
       // long nChgBasesCalledConfidently = 0;
        
        /** The number of Chh bases called confidently (according to user threshold), either ref or other */
       // long nChhBasesCalledConfidently = 0;
        
        long nCphBasesCalledConfidently = 0;
        
        /** The number of Gch bases called confidently (according to user threshold), either ref or other */
        long nGchBasesCalledConfidently = 0;
        
        /** The number of Gch bases called confidently (according to user threshold), either ref or other */
        long nCchBasesCalledConfidently = 0;
        
        /** The number of Gch bases called confidently (according to user threshold), either ref or other */
        long nWchBasesCalledConfidently = 0;
        
        /** The number of Gcg bases called confidently (according to user threshold), either ref or other */
        long nGcgBasesCalledConfidently = 0;
        
        /** The number of Hcg bases called confidently (according to user threshold), either ref or other */
        long nCcgBasesCalledConfidently = 0;
        
        /** The number of Wcg bases called confidently (according to user threshold), either ref or other */
        long nWcgBasesCalledConfidently = 0;
        
        
        /** The sum of methylation value of Cytosine bases called confidently (according to user threshold), either ref or other */
        double sumMethyCytosineBasesCalledConfidently = 0;
        
        /** The sum of methylation value of Cpg bases called confidently (according to user threshold), either ref or other */
        double sumMethyCpgBasesCalledConfidently = 0;
        
        /** The sum of methylation value of Chg bases called confidently (according to user threshold), either ref or other */
       // double sumMethyChgBasesCalledConfidently = 0;
        
        /** The sum of methylation value of Chh bases called confidently (according to user threshold), either ref or other */
       // double sumMethyChhBasesCalledConfidently = 0;
        
        double sumMethyCphBasesCalledConfidently = 0;
        
        /** The sum of methylation value of Gch bases called confidently (according to user threshold), either ref or other */
        double sumMethyGchBasesCalledConfidently = 0;
        
        /** The sum of methylation value of Gch bases called confidently (according to user threshold), either ref or other */
        double sumMethyCchBasesCalledConfidently = 0;
        
        /** The sum of methylation value of Gch bases called confidently (according to user threshold), either ref or other */
        double sumMethyWchBasesCalledConfidently = 0;
        
        /** The sum of methylation value of Gcg bases called confidently (according to user threshold), either ref or other */
        double sumMethyGcgBasesCalledConfidently = 0;
        
        /** The sum of methylation value of Hcg bases called confidently (according to user threshold), either ref or other */
        double sumMethyCcgBasesCalledConfidently = 0;
        
        /** The sum of methylation value of Wcg bases called confidently (according to user threshold), either ref or other */
        double sumMethyWcgBasesCalledConfidently = 0;

        double percentCallableOfAll()    { return (100.0 * nBasesCallable) / (nBasesVisited-nExtendedEvents); }
        double percentCalledOfAll()      { return (100.0 * nBasesCalledConfidently) / (nBasesVisited-nExtendedEvents); }
        double percentCalledOfCallable() { return (100.0 * nBasesCalledConfidently) / (nBasesCallable); }
        double percentMethyLevelOfC() { return (sumMethyCytosineBasesCalledConfidently) / (double)(nCytosineBasesCalledConfidently); }
        double percentMethyLevelOfCpg() { return (sumMethyCpgBasesCalledConfidently) / (double)(nCpgBasesCalledConfidently); }
       // double percentMethyLevelOfChh() { return (sumMethyChhBasesCalledConfidently) / (double)(nChhBasesCalledConfidently); }
       // double percentMethyLevelOfChg() { return (sumMethyChgBasesCalledConfidently) / (double)(nChgBasesCalledConfidently); }
        double percentMethyLevelOfCph() { return (sumMethyCphBasesCalledConfidently) / (double)(nCphBasesCalledConfidently); }
        
        double percentMethyLevelOfGch() { return (sumMethyGchBasesCalledConfidently) / (double)(nGchBasesCalledConfidently); }
        double percentMethyLevelOfCch() { return (sumMethyCchBasesCalledConfidently) / (double)(nCchBasesCalledConfidently); }
        double percentMethyLevelOfWch() { return (sumMethyWchBasesCalledConfidently) / (double)(nWchBasesCalledConfidently); }
        double percentMethyLevelOfGcg() { return (sumMethyGcgBasesCalledConfidently) / (double)(nGcgBasesCalledConfidently); }
        double percentMethyLevelOfCcg() { return (sumMethyCcgBasesCalledConfidently) / (double)(nCcgBasesCalledConfidently); }
        double percentMethyLevelOfWcg() { return (sumMethyWcgBasesCalledConfidently) / (double)(nWcgBasesCalledConfidently); }
        
        //record other kind of cytosine statistics status user provided
        HashMap<String, Double[]> otherCytosine = new HashMap<String, Double[]>();//value[0]: number of cytosine; value[1]: sumMethyLevel;
        
        void makeOtherCytosineMap(){
        	
        	if(!BAC.forceOtherCytosine.isEmpty()){
				
        		String[] tmpArray = BAC.forceOtherCytosine.split(";");
				for(String tmp : tmpArray){
					String[] tmp2 = tmp.split(":");
					String key = tmp2[0].toUpperCase();
					Double[] value = new Double[2];
					value[0] = 0.0;//value[0]: number of cytosine;
					value[1] = 0.0;//value[1]: sum of MethyLevel;
					otherCytosine.put(key,value);
					
				}
			}
        	else if(!BAC.autoEstimateOtherCytosine.isEmpty()){
        		String[] tmpArray = BAC.autoEstimateOtherCytosine.split(";");
        		for(String tmp : tmpArray){
					String[] tmp2 = tmp.split(":");
					String key = tmp2[0].toUpperCase();
					Double[] value = new Double[2];
					value[0] = 0.0;
					value[1] = 0.0;
					otherCytosine.put(key,value);
					
				}
        	}
        	else{
        		
        	}
        }
        
    }

    /**
     * Initialize the samples, output, and genotype calculation model
     *
     **/
    public void initialize() {

        Set<String> samples = new TreeSet<String>();
        //sometimes, BAM file also provided sample name, and it is different from user provided in the argument, then there will be an error~
        
        if ( BAC.ASSUME_SINGLE_SAMPLE != null ){
        	samples = SampleUtils.getSAMFileSamples(getToolkit().getSAMFileHeader());
        	if(!samples.isEmpty()){
        		System.out.println("sample name provided was masked by bam file header");
        	}
        	else{
        		samples.add(BAC.ASSUME_SINGLE_SAMPLE);
        	}
        }    
        else{
        	samples = SampleUtils.getSAMFileSamples(getToolkit().getSAMFileHeader());
        	System.out.println("samples provided: " + samples.toString());
        	if(samples.isEmpty()){
        		System.err.println("No sample name provided, program will automately provide the bam file header: " + samples.toString());
        		
        	}
        }
            
        
        //initiate BisulfiteGenotyperEngine
        BG_engine = new BisulfiteGenotyperEngine(getToolkit(), BAC, logger, samples);
        SAMSequenceDictionary refDict = getToolkit().getMasterSequenceDictionary();
        // initialize the header
        if(autoEstimateC){
        	if(secondIteration){
        		File outputVcfFile = new File(BAC.vfn1);
    			writer = new TcgaVCFWriter(outputVcfFile,refDict, false);
    			writer.setRefSource(getToolkit().getArguments().referenceFile.toString());
    			
    			writer.writeHeader(new VCFHeader(getHeaderInfo(), samples));
    			
    			if(getToolkit().getArguments().numberOfThreads > 1){
    				multiThreadWriter = new SortingTcgaVCFWriter(writer,MAXIMUM_CACHE_FOR_OUTPUT_VCF);
    				multiThreadWriter.enableDiscreteLoci(BAC.lnc);
    				if(BAC.ovd){
    					File outputVerboseFile = new File(BAC.fnovd);
        				verboseWriter = new TcgaVCFWriter(outputVerboseFile,refDict, false);
        				verboseWriter.writeHeader(new VCFHeader(getHeaderInfo(), samples));
    				}
    				
    			}
    			
        		if(BAC.OutputMode == BisulfiteGenotyperEngine.OUTPUT_MODE.DEFAULT_FOR_TCGA){
        			File outputAdditionalVcfFile = new File(BAC.vfn2);
        			additionalWriterForDefaultTcgaMode = new TcgaVCFWriter(outputAdditionalVcfFile,refDict, false);
        			
        			additionalWriterForDefaultTcgaMode.writeHeader(new VCFHeader(getHeaderInfo(), samples));
        			if(getToolkit().getArguments().numberOfThreads > 1){
        				multiAdditionalWriterForDefaultTcgaMode = new SortingTcgaVCFWriter(additionalWriterForDefaultTcgaMode,MAXIMUM_CACHE_FOR_OUTPUT_VCF);
        				multiAdditionalWriterForDefaultTcgaMode.enableDiscreteLoci(BAC.lnc);
        			}
        			
        		}
        			
        		if(BAC.orad){
            		File outputBamFile = new File(BAC.fnorad);
            		SAMFileWriterFactory samFileWriterFactory = new SAMFileWriterFactory();
            		samFileWriterFactory.setCreateIndex(true);
            		samWriter = samFileWriterFactory.makeBAMWriter(getToolkit().getSAMFileHeader(), false, outputBamFile);
            	}
        		if(BAC.fnocrd != null){
        			File outputReadsDetailFile = new File(BAC.fnocrd);
        			if(BAC.sequencingMode == MethylSNPModel.GM){
        				readsWriter = new NOMeSeqReadsWriterImp(outputReadsDetailFile);
        			}
        			else{
        				readsWriter = new cpgReadsWriterImp(outputReadsDetailFile);
        			}
        			
        			
        			if(getToolkit().getArguments().numberOfThreads > 1){
        				if(BAC.sequencingMode == MethylSNPModel.GM){
        					multiThreadCpgReadsWriter = new SortingNOMeSeqReadsWriter(readsWriter,MAXIMUM_CACHE_FOR_OUTPUT_READS);
            			}
            			else{
            				multiThreadCpgReadsWriter = new SortingCpgReadsWriter(readsWriter,MAXIMUM_CACHE_FOR_OUTPUT_READS);
            			}
        				multiThreadCpgReadsWriter.writeHeader(true);
        				//multiThreadCpgReadsWriter.enableDiscreteLoci(BAC.lnc);
        			}
        			else{
        				readsWriter.addHeader(true);
        			}
        		}
        		
        	}
        	else{
        		
        	}
        }
        else{
        	File outputVcfFile = new File(BAC.vfn1);
			writer = new TcgaVCFWriter(outputVcfFile,refDict, false);
			writer.setRefSource(getToolkit().getArguments().referenceFile.toString());
			writer.writeHeader(new VCFHeader(getHeaderInfo(), samples));
			
			if(getToolkit().getArguments().numberOfThreads > 1){
				multiThreadWriter = new SortingTcgaVCFWriter(writer,MAXIMUM_CACHE_FOR_OUTPUT_VCF);
				multiThreadWriter.enableDiscreteLoci(BAC.lnc);
			}
			
        	if(BAC.OutputMode == BisulfiteGenotyperEngine.OUTPUT_MODE.DEFAULT_FOR_TCGA){
    			File outputAdditionalVcfFile = new File(BAC.vfn2);
    			additionalWriterForDefaultTcgaMode = new TcgaVCFWriter(outputAdditionalVcfFile,refDict, false);
    			
    			additionalWriterForDefaultTcgaMode.writeHeader(new VCFHeader(getHeaderInfo(), samples));
    			
    			if(getToolkit().getArguments().numberOfThreads > 1){
    				multiAdditionalWriterForDefaultTcgaMode = new SortingTcgaVCFWriter(additionalWriterForDefaultTcgaMode,MAXIMUM_CACHE_FOR_OUTPUT_VCF);
    				multiAdditionalWriterForDefaultTcgaMode.enableDiscreteLoci(BAC.lnc);
    			}
    		}
        	
        	if(BAC.orad){
        		File outputBamFile = new File(BAC.fnorad);
        		SAMFileWriterFactory samFileWriterFactory = new SAMFileWriterFactory();
        		samFileWriterFactory.setCreateIndex(true);
        		samWriter = samFileWriterFactory.makeBAMWriter(getToolkit().getSAMFileHeader(), false, outputBamFile);
        	}
        	if(BAC.fnocrd != null){
    			File outputReadsDetailFile = new File(BAC.fnocrd);
    			if(BAC.sequencingMode == MethylSNPModel.GM){
    				readsWriter = new NOMeSeqReadsWriterImp(outputReadsDetailFile);
    			}
    			else{
    				readsWriter = new cpgReadsWriterImp(outputReadsDetailFile);
    			}
    			
    			if(getToolkit().getArguments().numberOfThreads > 1){
    				if(BAC.sequencingMode == MethylSNPModel.GM){
    					multiThreadCpgReadsWriter = new SortingNOMeSeqReadsWriter(readsWriter,MAXIMUM_CACHE_FOR_OUTPUT_READS);
        			}
        			else{
        				multiThreadCpgReadsWriter = new SortingCpgReadsWriter(readsWriter,MAXIMUM_CACHE_FOR_OUTPUT_READS);
        			}
    				multiThreadCpgReadsWriter.writeHeader(true);
    				//multiThreadCpgReadsWriter.enableDiscreteLoci(BAC.lnc);
    			}
    			else{
    				readsWriter.addHeader(true);
    			}
    		}
        }
        //in the first iteration, initiate CytosineTypeStatus
        if(!secondIteration)
        	summary = new CytosineTypeStatus(BAC);
        
    }

    /**
     * get VCF header for the output VCF file
     *
     **/
    private Set<VCFHeaderLine> getHeaderInfo() {
        Set<VCFHeaderLine> headerInfo = new HashSet<VCFHeaderLine>();

        if ( !BAC.NO_SLOD )
            headerInfo.add(new VCFInfoHeaderLine(VCFConstants.STRAND_BIAS_KEY, 1, VCFHeaderLineType.Float, "Strand Bias"));
        
        headerInfo.add(new VCFInfoHeaderLine(BisulfiteVCFConstants.CYTOSINE_TYPE, -1, VCFHeaderLineType.String, "Cytosine Type"));
        headerInfo.add(new VCFInfoHeaderLine(BisulfiteVCFConstants.GENOTYPE_TYPE, 1, VCFHeaderLineType.String, "Genotype Type"));
        headerInfo.add(new VCFInfoHeaderLine(BisulfiteVCFConstants.C_STRAND_KEY, 1, VCFHeaderLineType.String, "Cytosine in negative strand"));
        headerInfo.add(new VCFInfoHeaderLine(BisulfiteVCFConstants.NUMBER_OF_C_KEY, 1, VCFHeaderLineType.Integer, "number of C in this Cytosine position"));
        headerInfo.add(new VCFInfoHeaderLine(BisulfiteVCFConstants.NUMBER_OF_T_KEY, 1, VCFHeaderLineType.Integer, "number of T in this Cytosine position"));
        headerInfo.add(new VCFInfoHeaderLine(BisulfiteVCFConstants.CYTOSINE_METHY_VALUE, 1, VCFHeaderLineType.Float, "Methylation value in this Cytosine position"));
        
        //check in dbSNP or not
        if ( BAC.dbsnp.isBound() )
            headerInfo.add(new VCFInfoHeaderLine(VCFConstants.DBSNP_KEY, 0, VCFHeaderLineType.Flag, "dbSNP Membership"));
        
        

        // FORMAT and INFO fields
       // headerInfo.addAll(VCFUtils.getSupportedHeaderStrings(VCFConstants.GENOTYPE_POSTERIORS_KEY));
   

        // FILTER fields
        if ( BAC.STANDARD_CONFIDENCE_FOR_EMITTING < BAC.STANDARD_CONFIDENCE_FOR_CALLING )
            headerInfo.add(new VCFFilterHeaderLine(UnifiedGenotyperEngine.LOW_QUAL_FILTER_NAME, "Low quality"));
     
     // Program commandLine fields
        headerInfo.add(new VCFHeaderLine(BisulfiteVCFConstants.PROGRAM_ARGS, argCommandline));
        
        return headerInfo;
    }


    /**
     * Compute at a given locus.
     *
     * @param tracker the meta data tracker
     * @param refContext the reference base
     * @param rawContext contextual information around the locus
     * @return the VariantCallContext object
     */
    public BisulfiteVariantCallContext map(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext) {
    	
        CytosineTypeStatus cts = null;
        
        if(secondIteration){
            cts = summary.clone();
        }
        else{
        	cts = new CytosineTypeStatus(BAC);

        	if(BAC.orad)
        		downsamplingBamFile(rawContext);

        }

        
        BG_engine.setAutoParameters(autoEstimateC, secondIteration);
        BG_engine.setCytosineTypeStatus(cts);

        
        //calculation LikelihoodsAndGenotypes for this loci
    	return BG_engine.calculateLikelihoodsAndGenotypes(tracker, refContext, rawContext);
    }

    /**
     * Initiate statistics object.
     */
    public ContextCondition reduceInit() { 
    	ContextCondition initiated = new ContextCondition();
    	initiated.makeOtherCytosineMap();
    	return initiated; 
    	
    }

    public ContextCondition treeReduce(ContextCondition lhs, ContextCondition rhs) {
        lhs.nBasesCallable += rhs.nBasesCallable;
        lhs.nBasesCalledConfidently += rhs.nBasesCalledConfidently;
        lhs.nBasesVisited += rhs.nBasesVisited;
        lhs.nCallsMade += rhs.nCallsMade;
        lhs.nCytosineBasesCalledConfidently += rhs.nCytosineBasesCalledConfidently;
        if(BAC.sequencingMode == MethylSNPModel.BM){
        	lhs.nCpgBasesCalledConfidently += rhs.nCpgBasesCalledConfidently;
            // lhs.nChhBasesCalledConfidently += rhs.nChhBasesCalledConfidently;
            // lhs.nChgBasesCalledConfidently += rhs.nChgBasesCalledConfidently;
             lhs.nCphBasesCalledConfidently += rhs.nCphBasesCalledConfidently;
             lhs.sumMethyCpgBasesCalledConfidently += rhs.sumMethyCpgBasesCalledConfidently;
             // lhs.sumMethyChgBasesCalledConfidently += rhs.sumMethyChgBasesCalledConfidently;
             // lhs.sumMethyChhBasesCalledConfidently += rhs.sumMethyChhBasesCalledConfidently;
              lhs.sumMethyCphBasesCalledConfidently += rhs.sumMethyCphBasesCalledConfidently;
        }
        else if(BAC.sequencingMode == MethylSNPModel.GM){
        	lhs.nGchBasesCalledConfidently += rhs.nGchBasesCalledConfidently;
        	lhs.nCchBasesCalledConfidently += rhs.nCchBasesCalledConfidently;
        	lhs.nWchBasesCalledConfidently += rhs.nWchBasesCalledConfidently;
            lhs.nGcgBasesCalledConfidently += rhs.nGcgBasesCalledConfidently;
            lhs.nCcgBasesCalledConfidently += rhs.nCcgBasesCalledConfidently;
            lhs.nWcgBasesCalledConfidently += rhs.nWcgBasesCalledConfidently;
            lhs.sumMethyGchBasesCalledConfidently += rhs.sumMethyGchBasesCalledConfidently;
            lhs.sumMethyCchBasesCalledConfidently += rhs.sumMethyCchBasesCalledConfidently;
            lhs.sumMethyWchBasesCalledConfidently += rhs.sumMethyWchBasesCalledConfidently;
            lhs.sumMethyGcgBasesCalledConfidently += rhs.sumMethyGcgBasesCalledConfidently;
            lhs.sumMethyCcgBasesCalledConfidently += rhs.sumMethyCcgBasesCalledConfidently;
            lhs.sumMethyWcgBasesCalledConfidently += rhs.sumMethyWcgBasesCalledConfidently;
        }
        
        
        
        lhs.sumMethyCytosineBasesCalledConfidently += rhs.sumMethyCytosineBasesCalledConfidently;
        
        
        if(!rhs.otherCytosine.isEmpty()){
        	for(String key : rhs.otherCytosine.keySet()){
            	Double[] rhsValue = rhs.otherCytosine.get(key);
            	Double[] lhsValue = new Double[2];
            	if(lhs.otherCytosine.containsKey(key)){
            		lhsValue = lhs.otherCytosine.get(key);
            		lhsValue[0] += rhsValue[0];
            		lhsValue[1] += rhsValue[1];
            		lhs.otherCytosine.put(key,lhsValue);
            	}
            	
            }
        }
        
        
        return lhs;
    }

    
    /**
     * calculate statistics number in each reduce steps.
     */
    public ContextCondition reduce(BisulfiteVariantCallContext value, ContextCondition sum) {
    	// the base vistited
        sum.nBasesVisited++;

        // not call the locus when there is no coverage
        if ( value == null )
            return sum;

        System.err.println(value.vc.getGenotypes().size() + "\t" + value.vc.getGenotypes().getSampleNames());
        sum.nBasesCallable++;

        // the base was confidently callable
        sum.nBasesCalledConfidently += value.confidentlyCalled ? 1 : 0;
        
        CytosineTypeStatus ctsFirst = value.cts.get(value.cts.keySet().toArray()[0]);
     // the cytosine base was confidently callable
        sum.nCytosineBasesCalledConfidently += ctsFirst.isC ? 1 : 0;
        
        if(BAC.sequencingMode == MethylSNPModel.BM){
        	sum.nCpgBasesCalledConfidently += ctsFirst.isCpg ? 1 : 0;
            //sum.nChgBasesCalledConfidently += value.cts.isChg ? 1 : 0;
           // sum.nChhBasesCalledConfidently += value.cts.isChh ? 1 : 0;
            sum.nCphBasesCalledConfidently += ctsFirst.isCph ? 1 : 0;
            sum.sumMethyCpgBasesCalledConfidently += ctsFirst.isCpg ? ctsFirst.cpgMethyLevel : 0;
            //sum.sumMethyChgBasesCalledConfidently += ctsFirst.isChg ? ctsFirst.chgMethyLevel : 0;
           // sum.sumMethyChhBasesCalledConfidently += ctsFirst.isChh ? ctsFirst.chhMethyLevel : 0;
            sum.sumMethyCphBasesCalledConfidently += ctsFirst.isCph ? ctsFirst.cphMethyLevel : 0;
        }
        else if(BAC.sequencingMode == MethylSNPModel.GM){
        	  sum.nGchBasesCalledConfidently += ctsFirst.isGch ? 1 : 0;
        	  sum.nCchBasesCalledConfidently += ctsFirst.isCch ? 1 : 0;
        	  sum.nWchBasesCalledConfidently += ctsFirst.isWch ? 1 : 0;
              sum.nGcgBasesCalledConfidently += ctsFirst.isGcg ? 1 : 0;
              sum.nCcgBasesCalledConfidently += ctsFirst.isCcg ? 1 : 0;
              sum.nWcgBasesCalledConfidently += ctsFirst.isWcg ? 1 : 0;
              sum.sumMethyGchBasesCalledConfidently += ctsFirst.isGch ? ctsFirst.gchMethyLevel : 0;
              sum.sumMethyCchBasesCalledConfidently += ctsFirst.isCch ? ctsFirst.cchMethyLevel : 0;
              sum.sumMethyWchBasesCalledConfidently += ctsFirst.isWch ? ctsFirst.wchMethyLevel : 0;
              sum.sumMethyGcgBasesCalledConfidently += ctsFirst.isGcg ? ctsFirst.gcgMethyLevel : 0;
              sum.sumMethyCcgBasesCalledConfidently += ctsFirst.isCcg ? ctsFirst.ccgMethyLevel : 0;
              sum.sumMethyWcgBasesCalledConfidently += ctsFirst.isWcg ? ctsFirst.wcgMethyLevel : 0;
        }
        
      

        sum.sumMethyCytosineBasesCalledConfidently += ctsFirst.isC ? ctsFirst.cytosineMethyLevel : 0;
        
      
       // if(ctsFirst.isCpg){
       // 	System.err.println(value.vc.getChr() + "\t" + value.vc.getStart() + "\t" + sum.sumMethyCpgBasesCalledConfidently + "\t" + sum.nCpgBasesCalledConfidently + "\t" + ctsFirst.cpgMethyLevel);
      //  }
        
        //other cytosine provided.
        if(!BAC.forceOtherCytosine.isEmpty() || !BAC.autoEstimateOtherCytosine.isEmpty()){
        	for(String key : sum.otherCytosine.keySet()){
        		if(ctsFirst.cytosineListMap.containsKey(key)){
        			Double[] cytosineStatus = ctsFirst.cytosineListMap.get(key);
        			Double[] values = sum.otherCytosine.get(key);
        			values[0] += Double.compare(cytosineStatus[3], 1.0) == 0 ? 1: 0; //cytosineStatus[3]==1 means it is this type of cytosine pattern. otherwise, not..
        			values[1] += Double.compare(cytosineStatus[3], 1.0) == 0 ? cytosineStatus[2]: 0;
        			sum.otherCytosine.put(key,values); 			
        		}
        		
			}
        	
		}
    	
       
        // can't make a confident variant call here
        if ( value.vc == null)
            return sum;

        //System.out.println(value.vc.getStart());
        
        try {
            // actually making a call
            sum.nCallsMade++;
            if(BAC.sequencingMode == MethylSNPModel.GM){
    			if(BAC.fnocrd != null && ((autoEstimateC && secondIteration) || (!autoEstimateC && !secondIteration))){
    				boolean positiveStrand = true;
    				String sampleContext = "";
    				for(String key : ctsFirst.cytosineListMap.keySet()){
    					if(key.equalsIgnoreCase("C-1")){
    						continue;
    					}
    					Double[] cytosineStatus = ctsFirst.cytosineListMap.get(key);
    					if(Double.compare(cytosineStatus[3], 1.0) == 0){
    						if(cytosineStatus[0] > cytosineStatus[1]){
    							positiveStrand = true;
    						}
    						else{
    							positiveStrand = false;
    						}
    						String[] tmpKey = key.split("-");
    						if(sampleContext.equals("")){
    							sampleContext = tmpKey[0];
    						}
    						else{
    							sampleContext = sampleContext + "," + tmpKey[0];
    						}
    							
    						
    					}
    				}
    	        	if(!sampleContext.equals("")){
    	        		readsDetailReportForNOMeSeq(value.rawContext, positiveStrand, getToolkit().getArguments().numberOfThreads > 1, sampleContext, value.refBase);	
    	        	}
    	        	
    	        }
    		}
            else{
            	if(ctsFirst.isCpg){
        			if(BAC.fnocrd != null && ((autoEstimateC && secondIteration) || (!autoEstimateC && !secondIteration))){
        	        	if(ctsFirst.cytosineListMap.get("CG-1")[0] > ctsFirst.cytosineListMap.get("CG-1")[1]){
        	        		readsDetailReport(value.rawContext, true, getToolkit().getArguments().numberOfThreads > 1, value.refBase);
        	        	}
        	        	else{
        	        		readsDetailReport(value.rawContext, false, getToolkit().getArguments().numberOfThreads > 1, value.refBase);
        	        	}
        				
        	        	
        	        }
        		}
            }
            
            
            if(autoEstimateC){
            	if(secondIteration){ //in auto estimate mode, don't output in the first iteration	
            		if(BAC.OutputMode == BisulfiteGenotyperEngine.OUTPUT_MODE.EMIT_ALL_CYTOSINES){ // only output homozygous cytosine
            			if(ctsFirst.isC){
            				if(getToolkit().getArguments().numberOfThreads > 1){
                				multiThreadWriter.add(value.vc, value.refBase.getBase());
                				COUNT_CACHE_FOR_OUTPUT_VCF++;
                				if(COUNT_CACHE_FOR_OUTPUT_VCF % MAXIMUM_CACHE_FOR_OUTPUT_VCF ==0 ){
                					//writer.writeFlush();
                					//System.err.println("flush");
                			//		synchronized(multiThreadWriter){
                						multiThreadWriter.writerFlush();
                			//			COUNT_CACHE_FOR_OUTPUT_VCF=0;
                			//		}
                					System.gc();
                					
                				}
                			}
                			else{
                				writer.add(value.vc);
                			}
            			}
            				
            		}
            		else if(BAC.OutputMode == BisulfiteGenotyperEngine.OUTPUT_MODE.EMIT_ALL_CPG){ // only output homozygous cpg
            			if(ctsFirst.isCpg){
            				if(getToolkit().getArguments().numberOfThreads > 1){
                				multiThreadWriter.add(value.vc, value.refBase.getBase());
                				COUNT_CACHE_FOR_OUTPUT_VCF++;
                				if(COUNT_CACHE_FOR_OUTPUT_VCF % MAXIMUM_CACHE_FOR_OUTPUT_VCF ==0 ){
                					//writer.writeFlush();
                					//System.err.println("flush");
                			//		synchronized(multiThreadWriter){
                						multiThreadWriter.writerFlush();
                			//			COUNT_CACHE_FOR_OUTPUT_VCF=0;
                			//		}
                					System.gc();
                					
                				}
                			}
                			else{
                				writer.add(value.vc);
                			}
            			}
            				
            		}
            		else if(BAC.OutputMode == BisulfiteGenotyperEngine.OUTPUT_MODE.EMIT_VARIANTS_ONLY){ // only output variants
            			if(value.isVariant()){
            				if(getToolkit().getArguments().numberOfThreads > 1){
                				multiThreadWriter.add(value.vc, value.refBase.getBase());
                				COUNT_CACHE_FOR_OUTPUT_VCF++;
                				if(COUNT_CACHE_FOR_OUTPUT_VCF % MAXIMUM_CACHE_FOR_OUTPUT_VCF ==0 ){
                					//writer.writeFlush();
                					//System.err.println("flush");
                			//		synchronized(multiThreadWriter){
                						multiThreadWriter.writerFlush();
                			//			COUNT_CACHE_FOR_OUTPUT_VCF=0;
                			//		}
                					System.gc();
                					
                				}
                			}
                			else{
                				writer.add(value.vc);
                			}
            			}
            		}
            		else if(BAC.OutputMode == BisulfiteGenotyperEngine.OUTPUT_MODE.EMIT_HET_SNPS_ONLY){ // only output variants
            			if(value.isHetSnp()){
            				if(getToolkit().getArguments().numberOfThreads > 1){
                				multiThreadWriter.add(value.vc, value.refBase.getBase());
                				COUNT_CACHE_FOR_OUTPUT_VCF++;
                				if(COUNT_CACHE_FOR_OUTPUT_VCF % MAXIMUM_CACHE_FOR_OUTPUT_VCF ==0 ){
                					//writer.writeFlush();
                					//System.err.println("flush");
                			//		synchronized(multiThreadWriter){
                						multiThreadWriter.writerFlush();
                			//			COUNT_CACHE_FOR_OUTPUT_VCF=0;
                			//		}
                					System.gc();
                					
                				}
                			}
                			else{
                				writer.add(value.vc);
                			}
            			}
            		}
            		else if(BAC.OutputMode == BisulfiteGenotyperEngine.OUTPUT_MODE.DEFAULT_FOR_TCGA){ // only output variants
            			if(ctsFirst.isCpg){
            				if(getToolkit().getArguments().numberOfThreads > 1){
                				multiThreadWriter.add(value.vc, value.refBase.getBase());
                				COUNT_CACHE_FOR_OUTPUT_VCF++;
                				if(COUNT_CACHE_FOR_OUTPUT_VCF % MAXIMUM_CACHE_FOR_OUTPUT_VCF ==0 ){
                					//writer.writeFlush();
                					//System.err.println("flush");
                			//		synchronized(multiThreadWriter){
                						multiThreadWriter.writerFlush();
                			//			COUNT_CACHE_FOR_OUTPUT_VCF=0;
                			//		}
                					System.gc();
                					
                				}
                			}
                			else{
                				writer.add(value.vc);
                			}
            			}
            			if(value.isVariant()){
            				if(getToolkit().getArguments().numberOfThreads > 1){
                				multiAdditionalWriterForDefaultTcgaMode.add(value.vc, value.refBase.getBase());
                				COUNT_CACHE_FOR_OUTPUT_VCF++;
                				if(COUNT_CACHE_FOR_OUTPUT_VCF % MAXIMUM_CACHE_FOR_OUTPUT_VCF ==0 ){
                					//writer.writeFlush();
                					//System.err.println("flush");
                			//		synchronized(multiThreadWriter){
                					multiAdditionalWriterForDefaultTcgaMode.writerFlush();
                			//			COUNT_CACHE_FOR_OUTPUT_VCF=0;
                			//		}
                					System.gc();
                					
                				}
                			}
            				else{
            					additionalWriterForDefaultTcgaMode.add(value.vc);
            				}
            			}
            			
            		}
            		else{
            			if(getToolkit().getArguments().numberOfThreads > 1){
            				multiThreadWriter.add(value.vc, value.refBase.getBase());
            				
            				COUNT_CACHE_FOR_OUTPUT_VCF++;
            				if(COUNT_CACHE_FOR_OUTPUT_VCF % MAXIMUM_CACHE_FOR_OUTPUT_VCF ==0 ){
            					//writer.writeFlush();
            					//System.err.println("flush");
            			//		synchronized(multiThreadWriter){
            						multiThreadWriter.writerFlush();
            			//			COUNT_CACHE_FOR_OUTPUT_VCF=0;
            			//		}
            					System.gc();
            					
            				}
            				if(BAC.ovd){
            					verboseWriter.add(value.vc);
            				}
            				
            			}
            			else{
            				writer.add(value.vc);
            				
            			}
            			
            		}
            		
            	}
            }
            else{
            	if(BAC.OutputMode == BisulfiteGenotyperEngine.OUTPUT_MODE.EMIT_ALL_CYTOSINES){ // only output homozygous cytosine
        			if(ctsFirst.isC){
        				if(getToolkit().getArguments().numberOfThreads > 1){
            				multiThreadWriter.add(value.vc, value.refBase.getBase());
            				COUNT_CACHE_FOR_OUTPUT_VCF++;
            				if(COUNT_CACHE_FOR_OUTPUT_VCF % MAXIMUM_CACHE_FOR_OUTPUT_VCF ==0 ){
            					//writer.writeFlush();
            					//System.err.println("flush");
            			//		synchronized(multiThreadWriter){
            						multiThreadWriter.writerFlush();
            			//			COUNT_CACHE_FOR_OUTPUT_VCF=0;
            			//		}
            					System.gc();
            					
            				}
            			}
            			else{
            				writer.add(value.vc);
            			}
        				
        			}
        				
        		}
        		else if(BAC.OutputMode == BisulfiteGenotyperEngine.OUTPUT_MODE.EMIT_ALL_CPG){ // only output homozygous cpg
        			if(ctsFirst.isCpg){
        				if(getToolkit().getArguments().numberOfThreads > 1){
            				multiThreadWriter.add(value.vc, value.refBase.getBase());
            				COUNT_CACHE_FOR_OUTPUT_VCF++;
            				if(COUNT_CACHE_FOR_OUTPUT_VCF % MAXIMUM_CACHE_FOR_OUTPUT_VCF ==0 ){
            					//writer.writeFlush();
            					//System.err.println("flush");
            			//		synchronized(multiThreadWriter){
            						multiThreadWriter.writerFlush();
            			//			COUNT_CACHE_FOR_OUTPUT_VCF=0;
            			//		}
            					System.gc();
            					
            				}
            			}
            			else{
            				writer.add(value.vc);
            			}
        			}
        				
        		}
        		else if(BAC.OutputMode == BisulfiteGenotyperEngine.OUTPUT_MODE.EMIT_VARIANTS_ONLY){ // only output variants
        			if(value.isVariant()){
        				if(getToolkit().getArguments().numberOfThreads > 1){
            				multiThreadWriter.add(value.vc, value.refBase.getBase());
            				COUNT_CACHE_FOR_OUTPUT_VCF++;
            				if(COUNT_CACHE_FOR_OUTPUT_VCF % MAXIMUM_CACHE_FOR_OUTPUT_VCF ==0 ){
            					//writer.writeFlush();
            					//System.err.println("flush");
            			//		synchronized(multiThreadWriter){
            						multiThreadWriter.writerFlush();
            			//			COUNT_CACHE_FOR_OUTPUT_VCF=0;
            			//		}
            					System.gc();
            					
            				}
            			}
            			else{
            				writer.add(value.vc);
            			}
        			}
        		}
        		else if(BAC.OutputMode == BisulfiteGenotyperEngine.OUTPUT_MODE.EMIT_HET_SNPS_ONLY){ // only output variants
        			if(value.isHetSnp()){
        				if(getToolkit().getArguments().numberOfThreads > 1){
            				multiThreadWriter.add(value.vc, value.refBase.getBase());
            				COUNT_CACHE_FOR_OUTPUT_VCF++;
            				if(COUNT_CACHE_FOR_OUTPUT_VCF % MAXIMUM_CACHE_FOR_OUTPUT_VCF ==0 ){
            					//writer.writeFlush();
            					//System.err.println("flush");
            			//		synchronized(multiThreadWriter){
            						multiThreadWriter.writerFlush();
            			//			COUNT_CACHE_FOR_OUTPUT_VCF=0;
            			//		}
            					System.gc();
            					
            				}
            			}
            			else{
            				writer.add(value.vc);
            			}
        			}
        		}
        		else if(BAC.OutputMode == BisulfiteGenotyperEngine.OUTPUT_MODE.DEFAULT_FOR_TCGA){ // only output variants
        			if(ctsFirst.isCpg){
        				if(getToolkit().getArguments().numberOfThreads > 1){
            				multiThreadWriter.add(value.vc, value.refBase.getBase());
            				COUNT_CACHE_FOR_OUTPUT_VCF++;
            				if(COUNT_CACHE_FOR_OUTPUT_VCF % MAXIMUM_CACHE_FOR_OUTPUT_VCF ==0 ){
            					//writer.writeFlush();
            					//System.err.println("flush");
            			//		synchronized(multiThreadWriter){
            						multiThreadWriter.writerFlush();
            			//			COUNT_CACHE_FOR_OUTPUT_VCF=0;
            			//		}
            					System.gc();
            					
            				}
            			}
            			else{
            				writer.add(value.vc);
            			}
        			}
        			if(value.isVariant()){
        				
        				if(getToolkit().getArguments().numberOfThreads > 1){
            				multiAdditionalWriterForDefaultTcgaMode.add(value.vc, value.refBase.getBase());
            				COUNT_CACHE_FOR_OUTPUT_VCF++;
            				if(COUNT_CACHE_FOR_OUTPUT_VCF % MAXIMUM_CACHE_FOR_OUTPUT_VCF ==0 ){
            					//writer.writeFlush();
            					//System.err.println("flush");
            			//		synchronized(multiThreadWriter){
            					multiAdditionalWriterForDefaultTcgaMode.writerFlush();
            			//			COUNT_CACHE_FOR_OUTPUT_VCF=0;
            			//		}
            					System.gc();
            					
            				}
            			}
        				else{
        					additionalWriterForDefaultTcgaMode.add(value.vc);
        				}
        			}
        			
        		}
        		else{
        			if(getToolkit().getArguments().numberOfThreads > 1){
        				multiThreadWriter.add(value.vc, value.refBase.getBase());
        				COUNT_CACHE_FOR_OUTPUT_VCF++;
        				if(COUNT_CACHE_FOR_OUTPUT_VCF % MAXIMUM_CACHE_FOR_OUTPUT_VCF ==0 ){
        					//writer.writeFlush();
        					//System.err.println("flush");
        			//		synchronized(multiThreadWriter){
        						multiThreadWriter.writerFlush();
        			//			COUNT_CACHE_FOR_OUTPUT_VCF=0;
        			//		}
        					System.gc();
        					
        				}
        			}
        			else{
        				writer.add(value.vc);
        			}
        			
        		}
            }
            
        } catch (IllegalArgumentException e) {
            throw new IllegalArgumentException(e.getMessage() + "; this is often caused by using the --assume_single_sample_reads argument with the wrong sample name");
        }

        return sum;
    }

    public void onTraversalDone(ContextCondition sum) {
    	String fn = BAC.vfn1 + ".BisSNPMethySummarizeList.txt"; //write summarize cytosine methylation pattern information into a file. This is easy for BisSNP API to use.	
    	PrintWriter outWriter = null;
    	String outLine = null;
		try {
			outWriter = new PrintWriter(new File(fn));
			 
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    	
    	logger.info(String.format("Visited bases                                %d", sum.nBasesVisited));
        logger.info(String.format("Callable bases                               %d", sum.nBasesCallable));
        logger.info(String.format("Confidently called bases                     %d", sum.nBasesCalledConfidently));
        logger.info(String.format("%% callable bases of all loci                 %3.3f", sum.percentCallableOfAll()));
        logger.info(String.format("%% confidently called bases of all loci       %3.3f", sum.percentCalledOfAll()));
        logger.info(String.format("%% confidently called bases of callable loci  %3.3f", sum.percentCalledOfCallable()));
        logger.info(String.format("Actual calls made                            %d", sum.nCallsMade));
        logger.info(String.format("%% Methylation level of Cytosine loci       %3.3f", sum.percentMethyLevelOfC()));

        outLine = "C-1:" + sum.percentMethyLevelOfC() + "\t" + sum.nCytosineBasesCalledConfidently  + "\n";
        
        //logger.info(String.format("%% Methylation level of CHH loci       %3.3f", sum.percentMethyLevelOfChh()));
       // logger.info(String.format("%% Methylation level of CHG loci       %3.3f", sum.percentMethyLevelOfChg()));
        logger.info(String.format("%% number of Cytosine loci       %d", sum.nCytosineBasesCalledConfidently));
        if(BAC.sequencingMode == MethylSNPModel.BM){
        	logger.info(String.format("%% Methylation level of CpG loci       %3.3f", sum.percentMethyLevelOfCpg()));
        	logger.info(String.format("%% number of CpG loci       %d", sum.nCpgBasesCalledConfidently));
            logger.info(String.format("%% Methylation level of CpH loci       %3.3f", sum.percentMethyLevelOfCph()));
        	
            logger.info(String.format("%% number of CpH loci       %d", sum.nCphBasesCalledConfidently));
            outLine = outLine + "CG-1:" + sum.percentMethyLevelOfCpg() + "\t" + sum.nCpgBasesCalledConfidently  + "\n";
            outLine = outLine + "CH-1:" + sum.percentMethyLevelOfCph() + "\t" + sum.nCphBasesCalledConfidently  + "\n";
        }
        
     //   logger.info(String.format("%% sum of Cytosine methy       %3.3f", sum.sumMethyCytosineBasesCalledConfidently));
    //    logger.info(String.format("%% sum of CpG methy       %3.3f", sum.sumMethyCpgBasesCalledConfidently));
     //   logger.info(String.format("%% sum of CpH methy       %3.3f", sum.sumMethyCphBasesCalledConfidently));
       // logger.info(String.format("%% number of CHH loci       %d", sum.nChhBasesCalledConfidently));
       // logger.info(String.format("%% number of CHG loci       %d", sum.nChgBasesCalledConfidently));
        else if(BAC.sequencingMode == MethylSNPModel.GM){
        	logger.info(String.format("%% Methylation level of GCH loci       %3.3f", sum.percentMethyLevelOfGch()));
        	logger.info(String.format("%% number of GCH loci       %d", sum.nGchBasesCalledConfidently));
        	logger.info(String.format("%% Methylation level of CCH loci       %3.3f", sum.percentMethyLevelOfCch()));
        	logger.info(String.format("%% number of CCH loci       %d", sum.nCchBasesCalledConfidently));
        	logger.info(String.format("%% Methylation level of WCH loci       %3.3f", sum.percentMethyLevelOfWch()));
        	logger.info(String.format("%% number of WCH loci       %d", sum.nWchBasesCalledConfidently));
        	
        	logger.info(String.format("%% Methylation level of WCG loci       %3.3f", sum.percentMethyLevelOfWcg()));
            logger.info(String.format("%% number of WCG loci       %d", sum.nWcgBasesCalledConfidently));
            logger.info(String.format("%% Methylation level of CCG loci       %3.3f", sum.percentMethyLevelOfCcg()));
            logger.info(String.format("%% number of CCG loci       %d", sum.nCcgBasesCalledConfidently));
            logger.info(String.format("%% Methylation level of GCG loci       %3.3f", sum.percentMethyLevelOfGcg()));
            logger.info(String.format("%% number of GCG loci       %d", sum.nGcgBasesCalledConfidently));
            
            outLine = outLine + "GCH-2:" + sum.percentMethyLevelOfGch() + "\t" + sum.nGchBasesCalledConfidently  + "\n";
            outLine = outLine + "CCH-2:" + sum.percentMethyLevelOfCch() + "\t" + sum.nCchBasesCalledConfidently  + "\n";
            outLine = outLine + "WCH-2:" + sum.percentMethyLevelOfWch() + "\t" + sum.nWchBasesCalledConfidently  + "\n";
            
            outLine = outLine + "WCG-2:" + sum.percentMethyLevelOfWcg() + "\t" + sum.nWcgBasesCalledConfidently  + "\n";
            outLine = outLine + "CCG-2:" + sum.percentMethyLevelOfCcg() + "\t" + sum.nCcgBasesCalledConfidently  + "\n";
            outLine = outLine + "GCG-2:" + sum.percentMethyLevelOfGcg() + "\t" + sum.nGcgBasesCalledConfidently  + "\n";
        }
        if(!BAC.forceOtherCytosine.isEmpty() || !BAC.autoEstimateOtherCytosine.isEmpty()){
        	for(String key : sum.otherCytosine.keySet()){
        		Double[] values = sum.otherCytosine.get(key);
        		String[] tmp = key.split("-");
        		String name = tmp[0];
        //		logger.info(String.format("%% Methylation sum of %s loci       %3.3f", name, values[1]));
        		logger.info(String.format("%% Methylation level of %s loci       %3.3f", name, values[1]/values[0]));
        		logger.info(String.format("%% number of %s loci       %d", name, values[0].intValue()));
        		outLine = outLine + key + ":" + values[1]/values[0] + "\t" + values[0].intValue()  + "\n";
			}
		}
        
        summary.cytosineMethyLevel = Double.isNaN(sum.percentMethyLevelOfC()) ? 0 : sum.percentMethyLevelOfC();
        if(BAC.sequencingMode == MethylSNPModel.BM){
        	summary.cpgMethyLevel = Double.isNaN(sum.percentMethyLevelOfCpg()) ? 0 : sum.percentMethyLevelOfCpg();
            summary.cphMethyLevel = Double.isNaN(sum.percentMethyLevelOfCph()) ? 0 : sum.percentMethyLevelOfCph();
           // summary.chhMethyLevel = Double.isNaN(sum.percentMethyLevelOfChh()) ? 0 : sum.percentMethyLevelOfChh();
          //  summary.chgMethyLevel = Double.isNaN(sum.percentMethyLevelOfChg()) ? 0 : sum.percentMethyLevelOfChg();
        }
        else if(BAC.sequencingMode == MethylSNPModel.GM){
        	summary.gchMethyLevel = Double.isNaN(sum.percentMethyLevelOfGch()) ? 0 : sum.percentMethyLevelOfGch();
        	summary.cchMethyLevel = Double.isNaN(sum.percentMethyLevelOfCch()) ? 0 : sum.percentMethyLevelOfCch();
        	summary.wchMethyLevel = Double.isNaN(sum.percentMethyLevelOfWch()) ? 0 : sum.percentMethyLevelOfWch();
            summary.gcgMethyLevel = Double.isNaN(sum.percentMethyLevelOfGcg()) ? 0 : sum.percentMethyLevelOfGcg();
            summary.ccgMethyLevel = Double.isNaN(sum.percentMethyLevelOfCcg()) ? 0 : sum.percentMethyLevelOfCcg();
            summary.wcgMethyLevel = Double.isNaN(sum.percentMethyLevelOfWcg()) ? 0 : sum.percentMethyLevelOfWcg();
        }
        
        
        
        
        //copy cytosine statistics status for 2nd iteration to use
        for(String cytosineType : summary.cytosineListMap.keySet()){
        	String[] tmpKey = cytosineType.split("-");
			Double[] value = summary.cytosineListMap.get(cytosineType);
				if(tmpKey[0].equalsIgnoreCase("C")){
					value[2] = summary.cytosineMethyLevel;
				
				}
				
				else if(tmpKey[0].equalsIgnoreCase("CG")){
					value[2] = BAC.autoEstimateCpg ? summary.cpgMethyLevel : BAC.forceCpg;
					
				}
			//	else if(tmpKey[0].equalsIgnoreCase("CHH")){
			//		value[2] = BAC.autoEstimateChh ? summary.chhMethyLevel : BAC.forceChh;
			//	}
			//	else if(tmpKey[0].equalsIgnoreCase("CHG")){
			//		value[2] = BAC.autoEstimateChg ? summary.chgMethyLevel : BAC.forceChg;
			//	}
				else if(tmpKey[0].equalsIgnoreCase("CH")){
					value[2] = BAC.autoEstimateCph ? summary.cphMethyLevel : BAC.forceCph;
				}
				else if(tmpKey[0].equalsIgnoreCase("GCH")){
					value[2] = BAC.autoEstimateGch ? summary.gchMethyLevel : BAC.forceGch;
					
				}
				else if(tmpKey[0].equalsIgnoreCase("CCH")){
					value[2] = BAC.autoEstimateCch ? summary.cchMethyLevel : BAC.forceCch;
					
				}
				else if(tmpKey[0].equalsIgnoreCase("WCH")){
					value[2] = BAC.autoEstimateWch ? summary.wchMethyLevel : BAC.forceWch;
					
				}
				else if(tmpKey[0].equalsIgnoreCase("GCG")){
					value[2] = BAC.autoEstimateGcg ? summary.gcgMethyLevel : BAC.forceGcg;
					
				}
				else if(tmpKey[0].equalsIgnoreCase("CCG")){
					value[2] = BAC.autoEstimateCcg ? summary.ccgMethyLevel : BAC.forceCcg;
					
				}
				else if(tmpKey[0].equalsIgnoreCase("WCG")){
					value[2] = BAC.autoEstimateWcg ? summary.wcgMethyLevel : BAC.forceWcg;
					
				}
				else{
					if(!BAC.autoEstimateOtherCytosine.isEmpty()){
						String[] tmpArray = BAC.autoEstimateOtherCytosine.split(";");
						for(String tmp : tmpArray){
							String[] key = tmp.split(":");
							if(key[0].equalsIgnoreCase(cytosineType)){
								Double[] summaryStat = sum.otherCytosine.get(cytosineType);
								value[2] = summaryStat[1]/summaryStat[0];
							}
						}
					}
					if(!BAC.forceOtherCytosine.isEmpty()){
						String[] tmpArray = BAC.forceOtherCytosine.split(";");
						for(String tmp : tmpArray){
							String[] key = tmp.split(":");
							if(key[0].equalsIgnoreCase(cytosineType)){
								value[2] = Double.parseDouble(key[1]);
							}
						}
					}
				}
				if(value[2].isNaN())
					value[2] = 0.0;
				summary.cytosineListMap.put(cytosineType, value);
        }
        if(BAC.orad){
        	samWriter.close();
        }
        if(BAC.OutputMode == BisulfiteGenotyperEngine.OUTPUT_MODE.DEFAULT_FOR_TCGA){
        	//additionalWriterForDefaultTcgaMode.close();
        }
        if(getToolkit().getArguments().numberOfThreads > 1 && ((autoEstimateC && secondIteration) || (!autoEstimateC && !secondIteration))){
        	 multiThreadWriter.close();
        	 //verboseWriter.close();
        	 if(BAC.OutputMode == BisulfiteGenotyperEngine.OUTPUT_MODE.DEFAULT_FOR_TCGA){
     			multiAdditionalWriterForDefaultTcgaMode.close();
     		 }
 			
        }
        
        if(BAC.fnocrd != null){

			if((autoEstimateC && secondIteration) || (!autoEstimateC && !secondIteration)){
				if(getToolkit().getArguments().numberOfThreads > 1){
					multiThreadCpgReadsWriter.close();
				}
				else{
					readsWriter.close();
				}
				
			}
		}
       outWriter.println(outLine);
       outWriter.close();
    }
    
    //receive cytosine statistics status from main program
    public void setCytosineMethyStatus(CytosineTypeStatus sumCytosine) {
    	
    	summary = sumCytosine.clone();
    	if(BAC.sequencingMode == MethylSNPModel.BM){
    		BAC.forceCpg = BAC.autoEstimateCpg ? sumCytosine.cpgMethyLevel : BAC.forceCpg;
        	//BAC.forceChg = BAC.autoEstimateChg ? sumCytosine.chgMethyLevel : BAC.forceChg;
        	//BAC.forceChh = BAC.autoEstimateChh ? sumCytosine.chhMethyLevel : BAC.forceChh;
        	BAC.forceCph = BAC.autoEstimateCph ? sumCytosine.cphMethyLevel : BAC.forceCph;
		}
		else if(BAC.sequencingMode == MethylSNPModel.GM){
			BAC.forceGch = BAC.autoEstimateGch ? sumCytosine.gchMethyLevel : BAC.forceGch;
			BAC.forceCch = BAC.autoEstimateCch ? sumCytosine.cchMethyLevel : BAC.forceCch;
			BAC.forceWch = BAC.autoEstimateWch ? sumCytosine.wchMethyLevel : BAC.forceWch;
	    	BAC.forceGcg = BAC.autoEstimateGcg ? sumCytosine.gcgMethyLevel : BAC.forceGcg;
	    	BAC.forceCcg = BAC.autoEstimateCcg ? sumCytosine.ccgMethyLevel : BAC.forceCcg;
	    	BAC.forceWcg = BAC.autoEstimateWcg ? sumCytosine.wcgMethyLevel : BAC.forceWcg;
		}

    	for(String cytosineType : summary.cytosineListMap.keySet()){
			String[] tmpKey = cytosineType.split("-");
			Double[] value = summary.cytosineListMap.get(cytosineType);
				if(tmpKey[0].equalsIgnoreCase("C")){
					value[2] = sumCytosine.cytosineMethyLevel;
				
				}
				else if(tmpKey[0].equalsIgnoreCase("CG")){
					value[2] = BAC.autoEstimateCpg ? sumCytosine.cpgMethyLevel : BAC.forceCpg;
					
				}
			//	else if(tmpKey[0].equalsIgnoreCase("CHH")){
			//		value[2] = BAC.autoEstimateChh ? sumCytosine.chhMethyLevel : BAC.forceChh;
			//	}
			//	else if(tmpKey[0].equalsIgnoreCase("CHG")){
			//		value[2] = BAC.autoEstimateChg ? sumCytosine.chgMethyLevel : BAC.forceChg;
			//	}
				else if(tmpKey[0].equalsIgnoreCase("CH")){
					value[2] = BAC.autoEstimateCph ? summary.cphMethyLevel : BAC.forceCph;
				}
				else if(tmpKey[0].equalsIgnoreCase("GCH")){
					value[2] = BAC.autoEstimateGch ? sumCytosine.gchMethyLevel : BAC.forceGch;
					
				}
				else if(tmpKey[0].equalsIgnoreCase("CCH")){
					value[2] = BAC.autoEstimateCch ? sumCytosine.cchMethyLevel : BAC.forceCch;
					
				}
				else if(tmpKey[0].equalsIgnoreCase("WCH")){
					value[2] = BAC.autoEstimateWch ? sumCytosine.wchMethyLevel : BAC.forceWch;
					
				}
				else if(tmpKey[0].equalsIgnoreCase("GCG")){
					value[2] = BAC.autoEstimateGcg ? sumCytosine.gcgMethyLevel : BAC.forceGcg;
					
				}
				else if(tmpKey[0].equalsIgnoreCase("CCG")){
					value[2] = BAC.autoEstimateCcg ? sumCytosine.ccgMethyLevel : BAC.forceCcg;
					
				}
				else if(tmpKey[0].equalsIgnoreCase("WCG")){
					value[2] = BAC.autoEstimateWcg ? sumCytosine.wcgMethyLevel : BAC.forceWcg;
					
				}
				else{
					
					if(!BAC.forceOtherCytosine.isEmpty()){
						String[] tmpArray = BAC.forceOtherCytosine.split(";");
						for(String tmp : tmpArray){
							String[] key = tmp.split(":");
							if(key[0].equalsIgnoreCase(cytosineType)){
								value[2] = Double.parseDouble(key[1]);
							}
						}
					}
				}
				if(value[2].isNaN())
					value[2] = 0.0;
				
				summary.cytosineListMap.put(cytosineType, value);	

    	}
    	
    	
    }
    
    public void setAutoParameters(boolean autoEstimateC, boolean secondIteration, String argCommandline){
    	this.autoEstimateC = autoEstimateC;
    	this.secondIteration = secondIteration;
    	this.argCommandline = argCommandline;
    }
    
    public CytosineTypeStatus getCytosineMethyStatus() {
    	
		return summary.clone();
    }
    
    public TcgaVCFWriter getWriter(){
    	return writer;
    }
   
    public void setWriter(TcgaVCFWriter writer){
    	this.writer = writer;

    }
    /*
     * we use different downsampling strategy. e.g. when downsampling 10X, it randomly pick up s*10/r reads.(r is the mean coverage of the sample, we use 30 here for our sample, 
     * s is the total reads covered in this position) 
     * it does NOT treat reads like GATK, when reads number more than 10, cut off the reads number to 10, and keep the same when reads number lower than 10. 
     */
    public void downsamplingBamFile(AlignmentContext rawContext){
    	if(rawContext.hasReads()){
			String tag = "Xi";
			Integer coverageMarked = 0;
			//int covergaeLimit = getToolkit().getArguments().downsampleCoverage;
			int covergaeLimit = BAC.orcad;
			covergaeLimit = Math.max((covergaeLimit * rawContext.getBasePileup().depthOfCoverage())/SAMPLE_READS_MEAN_COVERAGE,1);
			//getToolkit().getArguments().downsampleCoverage = covergaeLimit;
			//covergaeLimit = (covergaeLimit * rawContext.getBasePileup().size())/SAMPLE_READS_MEAN_COVERAGE;
			//System.err.println("loc: " + rawContext.getLocation().getStart() + "\tcovergaeLimit: " + covergaeLimit + "\trawContext.getBasePileup().size(): " + rawContext.getBasePileup().size() + "\tdownsampleCoverage: " + getToolkit().getArguments().downsampleCoverage);
			
			//ReadBackedPileup downsampledPileup = rawContext.getBasePileup().getDownsampledPileup(covergaeLimit);
			ReadBackedPileup downsampledPileup = BisSNPUtils.getDownsampledPileup(rawContext.getBasePileup(), covergaeLimit);
			//if(rawContext.getBasePileup().size() < covergaeLimit){
			//	downsampledPileup = rawContext.getBasePileup();
			//}
			//else{
			//	downsampledPileup = rawContext.getBasePileup().getDownsampledPileup(covergaeLimit);
			//}
			// = rawContext.getBasePileup().getDownsampledPileup(covergaeLimit);
			for ( PileupElement p : rawContext.getBasePileup() ) {
				if(p.getRead().getIntegerAttribute(tag) != null){
					if(p.getRead().getIntegerAttribute(tag) == 2)
						//System.out.println("loc: " + rawContext.getLocation().getStart() + " tag: " + p.getRead().getIntegerAttribute(tag));
					if(p.getRead().getIntegerAttribute(tag) == 1)
						coverageMarked++;
				}
					
			}
			//System.out.println("loc: " + rawContext.getLocation().getStart() + " coverageMarked: " + coverageMarked);
			for ( PileupElement p : downsampledPileup ) {
				//System.out.println(p.toString());
				if(p.getRead().getIntegerAttribute(tag) != null){
					//System.out.println("loc: " + rawContext.getLocation().getStart() + " tag: " + p.getRead().getIntegerAttribute(tag));
				}
				if(coverageMarked >= covergaeLimit)
					break;
				if(p.getRead().getIntegerAttribute(tag) == null){
					samWriter.addAlignment(p.getRead());
    				p.getRead().setAttribute(tag, 1);
    				coverageMarked++;
				}
				
					
			}
			for ( PileupElement p : rawContext.getBasePileup() ) {
				if(p.getRead().getIntegerAttribute(tag) == null)
					p.getRead().setAttribute(tag, 2);
			}
			
		}

    }
    
    public void downsamplingBamFileLikeGATK(AlignmentContext rawContext){
    	if(rawContext.hasReads()){
			String tag = "Xi";
			Integer coverageMarked = 0;
			int covergaeLimit = getToolkit().getArguments().downsampleCoverage;
			ReadBackedPileup downsampledPileup = rawContext.getBasePileup().getDownsampledPileup(covergaeLimit);
			for ( PileupElement p : rawContext.getBasePileup() ) {
				if(p.getRead().getIntegerAttribute(tag) != null){
					if(p.getRead().getIntegerAttribute(tag) == 2)
						System.out.println("loc: " + rawContext.getLocation().getStart() + " tag: " + p.getRead().getIntegerAttribute(tag));
					if(p.getRead().getIntegerAttribute(tag) == 1)
						coverageMarked++;
				}
					
			}
			//System.out.println("loc: " + rawContext.getLocation().getStart() + " coverageMarked: " + coverageMarked);
			for ( PileupElement p : downsampledPileup ) {
				//System.out.println(p.toString());
				if(coverageMarked >= covergaeLimit)
					break;
				if(p.getRead().getIntegerAttribute(tag) == null){
					samWriter.addAlignment(p.getRead());
    				p.getRead().setAttribute(tag, 1);
    				coverageMarked++;
				}
				
					
			}
			for ( PileupElement p : rawContext.getBasePileup() ) {
				if(p.getRead().getIntegerAttribute(tag) == null)
					p.getRead().setAttribute(tag, 2);
			}
			
		}

    }
    
    public void readsDetailReport(AlignmentContext rawContext, boolean posStrand, boolean multiThread, ReferenceContext refContext){
    	if(rawContext.hasReads()){
			
			for ( PileupElement p : rawContext.getBasePileup() ) {
				GATKSAMRecordFilterStorage GATKrecordFilterStor = new GATKSAMRecordFilterStorage(p.getRead(), refContext, BAC);
				if(!GATKrecordFilterStor.isGoodBase(p.getOffset())){
					continue;
				}
				
				char strand;
				
				if(p.getRead().getReadNegativeStrandFlag()){
					strand = '-';
					if(!posStrand){
						char methyStatus;
						
						if(p.getBase()==BaseUtilsMore.G){
							methyStatus = 'm';
							cpgReads cr = new cpgReads(rawContext.getContig(), rawContext.getLocation().getStart(), methyStatus, p.getQual(), strand, p.getRead().getReadName());
							if(multiThread){
								multiThreadCpgReadsWriter.add(cr);
							}
							else{
								readsWriter.add(cr);
							}
							
						}
						else if(p.getBase()==BaseUtilsMore.A){
							methyStatus = 'u';
							cpgReads cr = new cpgReads(rawContext.getContig(), rawContext.getLocation().getStart(), methyStatus, p.getQual(), strand, p.getRead().getReadName());
							if(multiThread){
								multiThreadCpgReadsWriter.add(cr);
							}
							else{
								readsWriter.add(cr);
							}
						}

						
					}
				}
				else{
					strand = '+';
					if(posStrand){
						char methyStatus;
						if(p.getBase()==BaseUtilsMore.C){
							methyStatus = 'm';
							cpgReads cr = new cpgReads(rawContext.getContig(), rawContext.getLocation().getStart(), methyStatus, p.getQual(), strand, p.getRead().getReadName());
							if(multiThread){
								multiThreadCpgReadsWriter.add(cr);
							}
							else{
								readsWriter.add(cr);
							}
						}
						else if(p.getBase()==BaseUtilsMore.T){
							methyStatus = 'u';
							cpgReads cr = new cpgReads(rawContext.getContig(), rawContext.getLocation().getStart(), methyStatus, p.getQual(), strand, p.getRead().getReadName());
							if(multiThread){
								multiThreadCpgReadsWriter.add(cr);
							}
							else{
								readsWriter.add(cr);
							}
						}

					}
				}
				
	
			}

		}

    }
    
    public void readsDetailReportForNOMeSeq(AlignmentContext rawContext, boolean posStrand, boolean multiThread, String sampleContext, ReferenceContext refContext){
    	if(rawContext.hasReads()){
			String ref = "";
			if(posStrand){
				ref = ref + BaseUtilsMore.convertByteToString(refContext.getBase());
			}
			else{
				ref = ref + BaseUtilsMore.convertByteToString(BaseUtilsMore.iupacCodeComplement(refContext.getBase()));
			}
			GenomeLoc locPre = refContext.getGenomeLocParser().createGenomeLoc(refContext.getLocus().getContig(), refContext.getLocus().getStart()-1 );
			GenomeLoc locPost = refContext.getGenomeLocParser().createGenomeLoc(refContext.getLocus().getContig(), refContext.getLocus().getStart()+1 );
			
			if( refContext.getWindow().containsP(locPre) ){
				ReferenceContext tmpRef = new ReferenceContext(refContext.getGenomeLocParser(),locPre, refContext.getWindow(),refContext.getBases());
				if(posStrand){
					ref = BaseUtilsMore.convertByteToString(BaseUtilsMore.toIupacCodeNOMeSeqMode(tmpRef.getBase(),1)) + ref;
				}
				else{
					ref = ref + BaseUtilsMore.convertByteToString(BaseUtilsMore.toIupacCodeNOMeSeqMode(BaseUtilsMore.iupacCodeComplement(tmpRef.getBase()),3));
				}
				
			}
			if( refContext.getWindow().containsP(locPost) ){
				ReferenceContext tmpRef = new ReferenceContext(refContext.getGenomeLocParser(),locPost, refContext.getWindow(),refContext.getBases());
				if(posStrand){
					ref = ref + BaseUtilsMore.convertByteToString(BaseUtilsMore.toIupacCodeNOMeSeqMode(tmpRef.getBase(),3));
				}
				else{
					ref = BaseUtilsMore.convertByteToString(BaseUtilsMore.toIupacCodeNOMeSeqMode(BaseUtilsMore.iupacCodeComplement(tmpRef.getBase()),1)) + ref;
				}
			}
			
			
			
			for ( PileupElement p : rawContext.getBasePileup() ) {
				GATKSAMRecordFilterStorage GATKrecordFilterStor = new GATKSAMRecordFilterStorage(p.getRead(), refContext, BAC);
				if(!GATKrecordFilterStor.isGoodBase(p.getOffset())){
					continue;
				}
				
				char strand;
				
				COUNT_CACHE_FOR_OUTPUT_READS++;
				if(COUNT_CACHE_FOR_OUTPUT_READS % MAXIMUM_CACHE_FOR_OUTPUT_READS ==0 ){
					if(multiThread){
						//System.err.println("flush");
						multiThreadCpgReadsWriter.writerFlush();
					}
					else{
						readsWriter.writerFlush();
					}
					System.gc();
					
				}
				
				if(p.getRead().getReadNegativeStrandFlag()){
					strand = '-';
					
					if(!posStrand){
						char methyStatus;
						
						if(p.getBase()==BaseUtilsMore.G){
							methyStatus = 'm';
							NOMeSeqReads cr = new NOMeSeqReads(rawContext.getContig(), rawContext.getLocation().getStart(), sampleContext, ref, methyStatus, p.getQual(), strand, p.getRead().getReadName());
							if(multiThread){
								multiThreadCpgReadsWriter.add(cr);
							}
							else{
								((NOMeSeqReadsWriterImp)readsWriter).add(cr);
							}
							
						}
						else if(p.getBase()==BaseUtilsMore.A){
							methyStatus = 'u';
							NOMeSeqReads cr = new NOMeSeqReads(rawContext.getContig(), rawContext.getLocation().getStart(), sampleContext, ref, methyStatus, p.getQual(), strand, p.getRead().getReadName());
							if(multiThread){
								multiThreadCpgReadsWriter.add(cr);
							}
							else{
								readsWriter.add(cr);
							}
						}

						
					}
				}
				else{
					strand = '+';
					if(posStrand){
						char methyStatus;
						if(p.getBase()==BaseUtilsMore.C){
							methyStatus = 'm';
							NOMeSeqReads cr = new NOMeSeqReads(rawContext.getContig(), rawContext.getLocation().getStart(), sampleContext, ref, methyStatus, p.getQual(), strand, p.getRead().getReadName());
							if(multiThread){
								multiThreadCpgReadsWriter.add(cr);
							}
							else{
								readsWriter.add(cr);
							}
						}
						else if(p.getBase()==BaseUtilsMore.T){
							methyStatus = 'u';
							NOMeSeqReads cr = new NOMeSeqReads(rawContext.getContig(), rawContext.getLocation().getStart(), sampleContext, ref, methyStatus, p.getQual(), strand, p.getRead().getReadName());
							if(multiThread){
								multiThreadCpgReadsWriter.add(cr);
							}
							else{
								readsWriter.add(cr);
							}
						}

					}
				}
				
	
			}

		}

    }
    
    public static BisulfiteArgumentCollection getBAC(){
    	return BAC;
    }
    

  //  public boolean isReduceByInterval() {
   //     return true;
   // }

}
