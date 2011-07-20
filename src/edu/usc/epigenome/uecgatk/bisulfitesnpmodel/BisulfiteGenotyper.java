package edu.usc.epigenome.uecgatk.bisulfitesnpmodel;


import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import org.broad.tribble.vcf.VCFConstants;
import org.broad.tribble.vcf.VCFFilterHeaderLine;
import org.broad.tribble.vcf.VCFHeader;
import org.broad.tribble.vcf.VCFHeaderLine;
import org.broad.tribble.vcf.VCFHeaderLineType;
import org.broad.tribble.vcf.VCFInfoHeaderLine;
import org.broad.tribble.vcf.VCFWriter;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.DownsampleType;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.datasources.reads.SAMDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.helpers.DbSNPHelper;
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
import org.broadinstitute.sting.utils.vcf.VCFUtils;

import org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel.BaseUtilsMore.*;
import org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel.NonRefDependSNPGenotypeLikelihoodsCalculationModel.MethylSNPModel;


//TO DO
//mearue CpG or GpC not by ref seq, but by read status, like 90% of G in the next position of C for CpG
/**
 * A variant caller which unifies the approaches of several disparate callers.  Works for single-sample and
 * multi-sample data.  The user can choose from several different incorporated calculation models.
 */
// todo -- change when UG is generalized to do BAQ as necessary
//@BAQMode(QualityMode = BAQ.QualityMode.DONT_MODIFY, ApplicationTime = BAQ.ApplicationTime.HANDLED_IN_WALKER)
@BAQMode(QualityMode = BAQ.QualityMode.OVERWRITE_QUALS, ApplicationTime = BAQ.ApplicationTime.ON_INPUT)
@Reference(window=@Window(start=-200,stop=200))
@By(DataSource.REFERENCE)
@Downsample(by=DownsampleType.NONE)
public class BisulfiteGenotyper extends LocusWalker<BisulfiteVariantCallContext, BisulfiteGenotyper.BGStatistics> implements TreeReducible<BisulfiteGenotyper.BGStatistics> {


    @ArgumentCollection private BisulfiteArgumentCollection BAC = new BisulfiteArgumentCollection();

    // control the output
    @Output(doc="File to which variants should be written",required=true)
    protected VCFWriter writer = null;

    @Argument(fullName = "verbose_mode", shortName = "verbose", doc = "File to print all of the annotated and detailed debugging output", required = false)
    protected PrintStream verboseWriter = null;
    
    @Argument(fullName="annotation", shortName="A", doc="One or more specific annotations to apply to variant calls", required=false)
    protected List<String> annotationsToUse = new ArrayList<String>();

    @Argument(fullName="group", shortName="G", doc="One or more classes/groups of annotations to apply to variant calls", required=false)
    protected String[] annotationClassesToUse = { "Standard" };

    
    private static boolean autoEstimateC = false;
    private static boolean secondIteration = false;
    
    // the calculation arguments
    private BisulfiteGenotyperEngine BG_engine = null;

    // the annotation engine
    private VariantAnnotatorEngine annotationEngine;
    
    CytosineTypeStatus summary = new CytosineTypeStatus(BAC);

    // enable deletions in the pileup
    public boolean includeReadsWithDeletionAtLoci() { return true; }

    // enable extended events for indels
    public boolean generateExtendedEvents() { return BAC.GLmodel == GenotypeLikelihoodsCalculationModel.Model.DINDEL; }

    /**
     * Inner class for collecting output statistics from the UG
     */
    
    
    public static class BGStatistics {
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
        
        /** The number of Cpg bases called confidently (according to user threshold), either ref or other */
        long nCytosineBasesCalledConfidently = 0;
        
        /** The number of Cpg bases called confidently (according to user threshold), either ref or other */
        long nCpgBasesCalledConfidently = 0;
        
        /** The number of Chg bases called confidently (according to user threshold), either ref or other */
        long nChgBasesCalledConfidently = 0;
        
        /** The number of Chh bases called confidently (according to user threshold), either ref or other */
        long nChhBasesCalledConfidently = 0;
        
        /** The number of Gch bases called confidently (according to user threshold), either ref or other */
        long nGchBasesCalledConfidently = 0;
        
        /** The number of Gcg bases called confidently (according to user threshold), either ref or other */
        long nGcgBasesCalledConfidently = 0;
        
        /** The number of Hcg bases called confidently (according to user threshold), either ref or other */
        long nHcgBasesCalledConfidently = 0;
        
        
        /** The sum of methylation value of Cpg bases called confidently (according to user threshold), either ref or other */
        double sumMethyCytosineBasesCalledConfidently = 0;
        
        /** The sum of methylation value of Cpg bases called confidently (according to user threshold), either ref or other */
        double sumMethyCpgBasesCalledConfidently = 0;
        
        /** The sum of methylation value of Chg bases called confidently (according to user threshold), either ref or other */
        double sumMethyChgBasesCalledConfidently = 0;
        
        /** The sum of methylation value of Chh bases called confidently (according to user threshold), either ref or other */
        double sumMethyChhBasesCalledConfidently = 0;
        
        /** The sum of methylation value of Gch bases called confidently (according to user threshold), either ref or other */
        double sumMethyGchBasesCalledConfidently = 0;
        
        /** The sum of methylation value of Gcg bases called confidently (according to user threshold), either ref or other */
        double sumMethyGcgBasesCalledConfidently = 0;
        
        /** The sum of methylation value of Hcg bases called confidently (according to user threshold), either ref or other */
        double sumMethyHcgBasesCalledConfidently = 0;
  /*      
        HashMap<String, Double[]> otherCytosineListMap = new HashMap<String, Double[]>(); 
        
        if(!BAC.autoEstimateOtherCytosine.isEmpty()){
			String[] tmpArray = BAC.autoEstimateOtherCytosine.split(";");
			for(String tmp : tmpArray){
					Double[] value = new Double[2];
					value[0] = 0.0; //value[0]: number of cytosien; value[1]: sum of methylation level
					value[1] = 0.0;
					otherCytosineListMap.put(tmp, value);
			}
		}
      */  

        double percentCallableOfAll()    { return (100.0 * nBasesCallable) / (nBasesVisited-nExtendedEvents); }
        double percentCalledOfAll()      { return (100.0 * nBasesCalledConfidently) / (nBasesVisited-nExtendedEvents); }
        double percentCalledOfCallable() { return (100.0 * nBasesCalledConfidently) / (nBasesCallable); }
        double percentMethyLevelOfC() { return (sumMethyCytosineBasesCalledConfidently) / (double)(nCytosineBasesCalledConfidently); }
        double percentMethyLevelOfCpg() { return (sumMethyCpgBasesCalledConfidently) / (double)(nCpgBasesCalledConfidently); }
        double percentMethyLevelOfChh() { return (sumMethyChhBasesCalledConfidently) / (double)(nChhBasesCalledConfidently); }
        double percentMethyLevelOfChg() { return (sumMethyChgBasesCalledConfidently) / (double)(nChgBasesCalledConfidently); }
        double percentMethyLevelOfGch() { return (sumMethyGchBasesCalledConfidently) / (double)(nGchBasesCalledConfidently); }
        double percentMethyLevelOfGcg() { return (sumMethyGcgBasesCalledConfidently) / (double)(nGcgBasesCalledConfidently); }
        double percentMethyLevelOfHcg() { return (sumMethyHcgBasesCalledConfidently) / (double)(nHcgBasesCalledConfidently); }
    }

    /**
     * Initialize the samples, output, and genotype calculation model
     *
     **/
    public void initialize() {
        // get all of the unique sample names
        // if we're supposed to assume a single sample, do so
        Set<String> samples = new TreeSet<String>();
        if ( BAC.ASSUME_SINGLE_SAMPLE != null )
            samples.add(BAC.ASSUME_SINGLE_SAMPLE);
        else
            samples = SampleUtils.getSAMFileSamples(getToolkit().getSAMFileHeader());

        //System.err.println("writer-initial: " + writer.toString());
        if(autoEstimateC){
        	if(secondIteration){
        		//annotationEngine = new VariantAnnotatorEngine(getToolkit(), Arrays.asList(annotationClassesToUse), annotationsToUse);
        	}
        	else{
        		annotationEngine = new VariantAnnotatorEngine(getToolkit(), Arrays.asList(annotationClassesToUse), annotationsToUse);
        	}
        }
        else{
        	annotationEngine = new VariantAnnotatorEngine(getToolkit(), Arrays.asList(annotationClassesToUse), annotationsToUse);
        }
  
        BG_engine = new BisulfiteGenotyperEngine(getToolkit(), BAC, logger, annotationEngine, samples);
       
        // initialize the header
        if(autoEstimateC){
        	if(secondIteration){
        		writer.writeHeader(new VCFHeader(getHeaderInfo(), samples));
        		for(String key : summary.cytosineListMap.keySet()){
    				Double[] values = summary.cytosineListMap.get(key);
    				for(Double value : values){
    				//	System.err.println("summary.key: " + key + "\tsummary.value: " + value);
    				}
    			}
    			
    			//System.err.println(summary.chhMethyLevel);
    			//System.err.println(summary.chgMethyLevel);
    			//System.err.println(summary.cpgMethyLevel);
        	}
        	else{
        		//writer.writeHeader(new VCFHeader(getHeaderInfo(), samples));
        	}
        }
        else{
        	writer.writeHeader(new VCFHeader(getHeaderInfo(), samples));
        }
        	
        
    }

    private Set<VCFHeaderLine> getHeaderInfo() {
        Set<VCFHeaderLine> headerInfo = new HashSet<VCFHeaderLine>();

        // all annotation fields from VariantAnnotatorEngine
       // headerInfo.addAll(annotationEngine.getVCFAnnotationDescriptions());

        // annotation (INFO) fields from UnifiedGenotyper
        
        List<ReferenceOrderedDataSource> dataSources = getToolkit().getRodDataSources();
        for ( ReferenceOrderedDataSource source : dataSources ) {
            if ( source.getName().equals(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME) ) {
                headerInfo.add(new VCFInfoHeaderLine(VCFConstants.DBSNP_KEY, 0, VCFHeaderLineType.Flag, "dbSNP Membership"));
            }
            else if ( source.getName().startsWith(VariantAnnotatorEngine.dbPrefix) ) {
                String name = source.getName().substring(VariantAnnotatorEngine.dbPrefix.length());
                headerInfo.add(new VCFInfoHeaderLine(name, 0, VCFHeaderLineType.Flag, name + " Membership"));
            }
        }

        // FORMAT and INFO fields
        headerInfo.addAll(VCFUtils.getSupportedHeaderStrings(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY));

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
    /**
     * Compute at a given locus.
     *
     * @param tracker the meta data tracker
     * @param refContext the reference base
     * @param rawContext contextual information around the locus
     * @return the VariantCallContext object
     */
    public BisulfiteVariantCallContext map(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext) {
    	GenomeLoc thisLoc = refContext.getLocus();
    	int centerCoord = thisLoc.getStart();
    	String thisContig = refContext.getLocus().getContig();
        byte[] contextRef = 
   			this.getToolkit().getReferenceDataSource().getReference().getSubsequenceAt(thisContig, centerCoord-100, centerCoord+100).getBases();
        contextRef = BaseUtilsMore.toUpperCase(contextRef);
        
        CytosineTypeStatus cts = new CytosineTypeStatus(BAC);
        if(secondIteration){
        	cts = summary.clone();
        }
 	
        	BG_engine.setCytosineTypeStatus(cts, contextRef);
      
        
        
    	return BG_engine.calculateLikelihoodsAndGenotypes(tracker, refContext, rawContext);
    }

    public BGStatistics reduceInit() { return new BGStatistics(); }

    public BGStatistics treeReduce(BGStatistics lhs, BGStatistics rhs) {
        lhs.nBasesCallable += rhs.nBasesCallable;
        lhs.nBasesCalledConfidently += rhs.nBasesCalledConfidently;
        lhs.nBasesVisited += rhs.nBasesVisited;
        lhs.nCallsMade += rhs.nCallsMade;
        lhs.nCytosineBasesCalledConfidently += rhs.nCytosineBasesCalledConfidently;
        lhs.nCpgBasesCalledConfidently += rhs.nCpgBasesCalledConfidently;
        lhs.nChhBasesCalledConfidently += rhs.nChhBasesCalledConfidently;
        lhs.nChgBasesCalledConfidently += rhs.nChgBasesCalledConfidently;
        lhs.nGchBasesCalledConfidently += rhs.nGchBasesCalledConfidently;
        lhs.nGcgBasesCalledConfidently += rhs.nGcgBasesCalledConfidently;
        lhs.nHcgBasesCalledConfidently += rhs.nHcgBasesCalledConfidently;
        
        lhs.sumMethyCytosineBasesCalledConfidently += rhs.sumMethyCytosineBasesCalledConfidently;
        lhs.sumMethyCpgBasesCalledConfidently += rhs.sumMethyCpgBasesCalledConfidently;
        lhs.sumMethyChgBasesCalledConfidently += rhs.sumMethyChgBasesCalledConfidently;
        lhs.sumMethyChhBasesCalledConfidently += rhs.sumMethyChhBasesCalledConfidently;
        lhs.sumMethyGchBasesCalledConfidently += rhs.sumMethyGchBasesCalledConfidently;
        lhs.sumMethyGcgBasesCalledConfidently += rhs.sumMethyGcgBasesCalledConfidently;
        lhs.sumMethyHcgBasesCalledConfidently += rhs.sumMethyHcgBasesCalledConfidently;
        
        return lhs;
    }

    public BGStatistics reduce(BisulfiteVariantCallContext value, BGStatistics sum) {
        // we get a point for reaching reduce
        sum.nBasesVisited++;

        // can't call the locus because of no coverage
        if ( value == null )
            return sum;

        // A call was attempted -- the base was potentially callable
        sum.nBasesCallable++;

        // the base was confidently callable
        sum.nBasesCalledConfidently += value.confidentlyCalled ? 1 : 0;
        
        sum.nCytosineBasesCalledConfidently += value.cts.isC ? 1 : 0;
        sum.nCpgBasesCalledConfidently += value.cts.isCpg ? 1 : 0;
        sum.nChgBasesCalledConfidently += value.cts.isChg ? 1 : 0;
        sum.nChhBasesCalledConfidently += value.cts.isChh ? 1 : 0;
        sum.nGchBasesCalledConfidently += value.cts.isGch ? 1 : 0;
        sum.nGcgBasesCalledConfidently += value.cts.isGcg ? 1 : 0;
        sum.nHcgBasesCalledConfidently += value.cts.isHcg ? 1 : 0;

        sum.sumMethyCytosineBasesCalledConfidently += value.cts.isC ? value.cts.cytosineMethyLevel : 0;
        sum.sumMethyCpgBasesCalledConfidently += value.cts.isCpg ? value.cts.cpgMethyLevel : 0;
        sum.sumMethyChgBasesCalledConfidently += value.cts.isChg ? value.cts.chgMethyLevel : 0;
        sum.sumMethyChhBasesCalledConfidently += value.cts.isChh ? value.cts.chhMethyLevel : 0;
        sum.sumMethyGchBasesCalledConfidently += value.cts.isGch ? value.cts.gchMethyLevel : 0;
        sum.sumMethyGcgBasesCalledConfidently += value.cts.isGcg ? value.cts.gcgMethyLevel : 0;
        sum.sumMethyHcgBasesCalledConfidently += value.cts.isHcg ? value.cts.hcgMethyLevel : 0;
        if(sum.sumMethyChhBasesCalledConfidently > 0){
        //	System.err.println(sum.nChhBasesCalledConfidently);
        //    System.err.println(sum.sumMethyChhBasesCalledConfidently);
        //    System.err.println(value.cts.isChh);
        //    System.err.println(value.cts.isCpg);
        }
        
        // can't make a confident variant call here
        if ( value.vc == null )
            return sum;

        try {
            // we are actually making a call
            sum.nCallsMade++;
            if(autoEstimateC){
            	if(secondIteration){
            		writer.add(value.vc, value.refBase);
            	}
            }
            else{
            	writer.add(value.vc, value.refBase);
            }
            
        } catch (IllegalArgumentException e) {
            throw new IllegalArgumentException(e.getMessage() + "; this is often caused by using the --assume_single_sample_reads argument with the wrong sample name");
        }

        return sum;
    }

    public void onTraversalDone(BGStatistics sum) {
        logger.info(String.format("Visited bases                                %d", sum.nBasesVisited));
        logger.info(String.format("Callable bases                               %d", sum.nBasesCallable));
        logger.info(String.format("Confidently called bases                     %d", sum.nBasesCalledConfidently));
        logger.info(String.format("%% callable bases of all loci                 %3.3f", sum.percentCallableOfAll()));
        logger.info(String.format("%% confidently called bases of all loci       %3.3f", sum.percentCalledOfAll()));
        logger.info(String.format("%% confidently called bases of callable loci  %3.3f", sum.percentCalledOfCallable()));
        logger.info(String.format("Actual calls made                            %d", sum.nCallsMade));
        logger.info(String.format("%% Methylation level of Cytosine loci       %3.3f", sum.percentMethyLevelOfC()));
        logger.info(String.format("%% Methylation level of CpG loci       %3.3f", sum.percentMethyLevelOfCpg()));
        logger.info(String.format("%% Methylation level of CHH loci       %3.3f", sum.percentMethyLevelOfChh()));
        logger.info(String.format("%% Methylation level of CHG loci       %3.3f", sum.percentMethyLevelOfChg()));
        logger.info(String.format("%% number of Cytosine loci       %d", sum.nCytosineBasesCalledConfidently));
        logger.info(String.format("%% number of CpG loci       %d", sum.nCpgBasesCalledConfidently));
        logger.info(String.format("%% number of CHH loci       %d", sum.nChhBasesCalledConfidently));
        logger.info(String.format("%% number of CHG loci       %d", sum.nChgBasesCalledConfidently));
        if(BAC.sequencingMode == MethylSNPModel.GM){
        	logger.info(String.format("%% Methylation level of GCH loci       %3.3f", sum.percentMethyLevelOfGch()));
            logger.info(String.format("%% Methylation level of GCG loci       %3.3f", sum.percentMethyLevelOfGcg()));
            logger.info(String.format("%% Methylation level of HCG loci       %3.3f", sum.percentMethyLevelOfHcg()));
        }
        summary.cytosineMethyLevel = Double.isNaN(sum.percentMethyLevelOfC()) ? 0 : sum.percentMethyLevelOfC();
        summary.cpgMethyLevel = Double.isNaN(sum.percentMethyLevelOfCpg()) ? 0 : sum.percentMethyLevelOfCpg();
        summary.chhMethyLevel = Double.isNaN(sum.percentMethyLevelOfChh()) ? 0 : sum.percentMethyLevelOfChh();
        summary.chgMethyLevel = Double.isNaN(sum.percentMethyLevelOfChg()) ? 0 : sum.percentMethyLevelOfChg();
        summary.gchMethyLevel = Double.isNaN(sum.percentMethyLevelOfGch()) ? 0 : sum.percentMethyLevelOfGch();
        summary.gcgMethyLevel = Double.isNaN(sum.percentMethyLevelOfGcg()) ? 0 : sum.percentMethyLevelOfGcg();
        summary.hcgMethyLevel = Double.isNaN(sum.percentMethyLevelOfHcg()) ? 0 : sum.percentMethyLevelOfHcg();
        if(secondIteration){
        	writer.close();
        }

    }
    
    public void setCytosineMethyStatus(CytosineTypeStatus sumCytosine) {
    	
    	summary = sumCytosine.clone();
    	
    	BAC.forceCpg = BAC.autoEstimateCpg ? sumCytosine.cpgMethyLevel : BAC.forceCpg;
    	BAC.forceChg = BAC.autoEstimateChg ? sumCytosine.chgMethyLevel : BAC.forceChg;
    	BAC.forceChh = BAC.autoEstimateChh ? sumCytosine.chhMethyLevel : BAC.forceChh;
    	BAC.forceGch = BAC.autoEstimateGch ? sumCytosine.gchMethyLevel : BAC.forceGch;
    	BAC.forceGcg = BAC.autoEstimateGcg ? sumCytosine.gcgMethyLevel : BAC.forceGcg;
    	BAC.forceHcg = BAC.autoEstimateHcg ? sumCytosine.hcgMethyLevel : BAC.forceHcg;
    	
    	for(String cytosineType : summary.cytosineListMap.keySet()){
			String[] tmpKey = cytosineType.split("-");
			Double[] value = summary.cytosineListMap.get(cytosineType);
			//System.err.println("ctype: " + tmpKey[0]);
			
			
				
				//System.err.println("ctype: " + tmpKey[0]);
				if(tmpKey[0].equalsIgnoreCase("C")){
					value[2] = sumCytosine.cytosineMethyLevel;
				
				}
				else if(tmpKey[0].equalsIgnoreCase("CG")){
					value[2] = BAC.autoEstimateCpg ? sumCytosine.cpgMethyLevel : BAC.forceCpg;
					
				}
				else if(tmpKey[0].equalsIgnoreCase("CHH")){
					value[2] = BAC.autoEstimateChh ? sumCytosine.chhMethyLevel : BAC.forceChh;
				}
				else if(tmpKey[0].equalsIgnoreCase("CHG")){
					value[2] = BAC.autoEstimateChg ? sumCytosine.chgMethyLevel : BAC.forceChg;
				}
				else if(tmpKey[0].equalsIgnoreCase("GCH")){
					value[2] = BAC.autoEstimateGch ? sumCytosine.gchMethyLevel : BAC.forceGch;
					
				}
				else if(tmpKey[0].equalsIgnoreCase("GCG")){
					value[2] = BAC.autoEstimateGcg ? sumCytosine.gcgMethyLevel : BAC.forceGcg;
					
				}
				else if(tmpKey[0].equalsIgnoreCase("HCG")){
					value[2] = BAC.autoEstimateHcg ? sumCytosine.hcgMethyLevel : BAC.forceHcg;
					
				}
				else{
					if(!BAC.autoEstimateOtherCytosine.isEmpty()){
						String[] tmpArray = BAC.autoEstimateOtherCytosine.split(";");
						for(String tmp : tmpArray){
							if(tmp.equalsIgnoreCase(cytosineType)){
								value[2] = sumCytosine.cytosineListMap.get(tmp)[2];
							}
						}
					}
				}
				if(value[2].isNaN())
					value[2] = 0.0;
				summary.cytosineListMap.put(cytosineType, value);	
			
    	}
    	
    	
    }
    
    public void setAutoParameters(boolean autoEstimateC, boolean secondIteration){
    	this.autoEstimateC = autoEstimateC;
    	this.secondIteration = secondIteration;
    }
    
    public CytosineTypeStatus getCytosineMethyStatus() {
    	
		return summary.clone();
    }
    
    public VCFWriter getWriter(){
    	return writer;
    }
   
    public void setWriter(VCFWriter writer){
    	this.writer = writer;
    	//System.err.println("writer-setup: " + this.writer.toString());
    }
    
    public VariantAnnotatorEngine getAnnoEng(){
    	return annotationEngine;
    }
   
    public void setAnnoEng(VariantAnnotatorEngine AnnoEng){
    	this.annotationEngine = AnnoEng;
    	//System.err.println("writer-setup: " + this.writer.toString());
    }
}
