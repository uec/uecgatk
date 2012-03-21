package edu.usc.epigenome.uecgatk.BisSNP;

import java.io.FileNotFoundException;

import java.util.ArrayList;
import java.util.BitSet;

import java.util.Collection;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.StringUtil;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.sting.utils.variantcontext.GenotypesContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContext.Validation;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;

import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;

import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.QualityUtils;

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import edu.usc.epigenome.uecgatk.BisSNP.BadBaseFilterBisulfite;

import edu.usc.epigenome.uecgatk.BisSNP.BisulfiteAlignmentUtils;
import edu.usc.epigenome.uecgatk.BisSNP.BisulfiteEnums.OUTPUT_MODE;
import edu.usc.epigenome.uecgatk.BisSNP.BisulfiteSNPGenotypeLikelihoodsCalculationModel;
import edu.usc.epigenome.uecgatk.BisSNP.BisulfiteVCFConstants;

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

public class BisulfiteGenotyperEngine{

	private BisulfiteArgumentCollection BAC = null;
	
	private BisulfiteVariantCallContext bisulfiteVariantCallContext = null;
	
	private GenomeAnalysisEngine toolkit = null;
	
	private UnifiedGenotyperEngine UGE = null;
	//cytosine pattern and their status
	//private ThreadLocal<HashMap<String,CytosineTypeStatus>> ctss = new ThreadLocal<HashMap<String,CytosineTypeStatus>>();
	//private ThreadLocal<CytosineTypeStatus> ctss = new ThreadLocal<CytosineTypeStatus>();
	
    // the model used for calculating genotypes
    private ThreadLocal<BisulfiteSNPGenotypeLikelihoodsCalculationModel> bglcms = new ThreadLocal<BisulfiteSNPGenotypeLikelihoodsCalculationModel>();

    // the allele frequency likelihoods (allocated once as an optimization)
  //  private ThreadLocal<double[]> log10AlleleFrequencyPosteriors = new ThreadLocal<double[]>();

    // the priors object
 //   private GenotypePriors genotypePriors;

    // samples in input
  //  private Set<String> samples = new TreeSet<String>();

    // the various loggers and writers
 //   private Logger logger = null;

    // fasta reference reader to supplement the edges of the reference sequence for long reads
  //  private IndexedFastaSequenceFile referenceReader;

    // the standard filter to use for calls below the confidence threshold but above the emit threshold
 //   private static final Set<String> filter = new HashSet<String>(1);

 //   private static final int MISMATCH_WINDOW_SIZE = 20;
    
    private static boolean autoEstimateC = false;
    private static boolean secondIteration = false;

	protected double MAX_PHRED = 1000000;
	
	//public static final String LOW_QUAL_FILTER_NAME = "LowQual";
	

	public BisulfiteGenotyperEngine(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext, BisulfiteArgumentCollection BAC, GenomeAnalysisEngine toolkit, boolean autoEstimateC, boolean secondIteration) {
		this.BAC = BAC.clone();
		this.toolkit = toolkit;
		this.autoEstimateC = autoEstimateC;
		this.secondIteration = secondIteration;
		//filter.add(LOW_QUAL_FILTER_NAME);
		this.UGE = initializeUnifiedGenotypeEngine();
		// TODO Auto-generated constructor stub
		calculateLikelihoodsAndGenotypes(tracker, refContext, rawContext);
		
	}
	
	public UnifiedGenotyperEngine initializeUnifiedGenotypeEngine(){
		return new UnifiedGenotyperEngine(toolkit, (UnifiedArgumentCollection)BAC);
	}

	
	/*
	public void setCytosineTypeStatus(CytosineTypeStatus cts){
		HashMap<String,CytosineTypeStatus> cytosineStatuss = new HashMap<String,CytosineTypeStatus>();
		for(String sample : samples){
			cytosineStatuss.put(sample, cts.clone());
		}
		this.ctss.set(cytosineStatuss);
	}
	*/
	public BisulfiteVariantCallContext getBisulfiteVariantCallContext(){
		
		return bisulfiteVariantCallContext;
	}
	
	/**
     * Compute full BisulfiteVariantCallContext at a given locus.
     *
     * @param tracker    the meta data tracker
     * @param refContext the reference base
     * @param rawContext contextual information around the locus
     * @return the BisulfiteVariantCallContext object
     */
    public void calculateLikelihoodsAndGenotypes(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext) {
        Map<String, AlignmentContext> stratifiedContexts = AlignmentContextUtils.splitContextBySampleName(rawContext);
        HashMap<String, BisulfiteContextsGenotypeLikelihoods> BCGLs = new HashMap<String, BisulfiteContextsGenotypeLikelihoods>();
        VariantContext vc = calculateLikelihoods(tracker, refContext, stratifiedContexts, AlignmentContextUtils.ReadOrientation.COMPLETE, null, BCGLs); //get likelihood and methylation pattern information from BisulfiteSNPGenotypeLikelihoodsCalculationModel
        
        if ( vc == null )
            return;
        
        VariantCallContext vcc = UGE.calculateGenotypes(tracker, refContext, rawContext, vc);
        //calculateGenotypes(tracker, refContext, rawContext, stratifiedContexts, BCGLs, vc); //including all Reads group genotypes information into vcc, and provide most probable genotype and alt-allele for all of ReadsGroup 
        
        bisulfiteVariantCallContext = calculateBisulfiteVariantCallContext(BCGLs, vcc, rawContext, refContext); // including bisulfite likelihood information into it
        
        //return vcc;
    }

    

    public VariantContext calculateLikelihoods(RefMetaDataTracker tracker, ReferenceContext refContext, Map<String, AlignmentContext> stratifiedContexts, AlignmentContextUtils.ReadOrientation type, Allele alternateAlleleToUse, Map<String, BisulfiteContextsGenotypeLikelihoods> BCGLs) {
		if ( stratifiedContexts == null ){
			 return null;
		}

        // initialize the data for this thread if that hasn't been done yet
        if ( bglcms.get() == null ) {
            bglcms.set(BisulfiteGenotyperEngine.getGenotypeLikelihoodsCalculationObject(BAC, autoEstimateC, secondIteration));
           // if(refContext.getLocus().getStart() == BAC.testLocus){
            //	System.err.println("initiate\t" + refContext.getLocus().getStart());
           // }
        }

        BisulfiteSNPGenotypeLikelihoodsCalculationModel bglcm = (BisulfiteSNPGenotypeLikelihoodsCalculationModel) bglcms.get();
        //bglcm.initialize();
        
        bglcm.setBsLikelihoods(tracker, refContext, stratifiedContexts, type, alternateAlleleToUse);
        BCGLs = bglcm.getBsContextGenotypeLikelihoods();
       // ctss.set(bglcm.getCytosineTypeStatus());
        
        if (bglcm.getRefAllele() != null)
            return createVariantContextFromLikelihoods(refContext, bglcm.getRefAllele(), BCGLs);
        else
            return null;
          
        
    }
    
    public BisulfiteVariantCallContext calculateBisulfiteVariantCallContext(HashMap<String, BisulfiteContextsGenotypeLikelihoods> BCGLs, VariantCallContext vcc, AlignmentContext rawContext, ReferenceContext refContext){
    	
    	return new BisulfiteVariantCallContext(BCGLs, vcc, rawContext, refContext);
    	
    	
    }
    
    /*
    protected void assignAFPosteriors(double[]likelihoods, double[] log10AFPosteriors){
    	for(int i = 0; i < likelihoods.length; i++){
    		log10AFPosteriors[i] = likelihoods[i];
    	}
    		
    }
	*/
/*
	protected VariantCallContext calculateGenotypes(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext, Map<String, AlignmentContext> stratifiedContexts, Map<String, BisulfiteContextsGenotypeLikelihoods> BCGLs, VariantContext vc) {
		// initialize the data for this thread if that hasn't been done yet
        if ( bglcms.get() == null ) {
            return null;
        }
        BisulfiteVariantCallContext call = null;
        log10AlleleFrequencyPosteriors.set(new double[3]);
        
        for ( String sample : BCGLs.keySet() ) {
        	BisulfiteContextsGenotypeLikelihoods GL = BCGLs.get(sample);
        	
        	assignAFPosteriors(GL.getLikelihoods(),log10AlleleFrequencyPosteriors.get());
        	int bestAFguess = MathUtils.maxElementIndex(log10AlleleFrequencyPosteriors.get());
            int secondAFguess = MathUtils.minElementIndex(log10AlleleFrequencyPosteriors.get());
            for (int i = 0; i < log10AlleleFrequencyPosteriors.get().length; i++){
            	if(i != bestAFguess){	
            		if(log10AlleleFrequencyPosteriors.get()[i] >= log10AlleleFrequencyPosteriors.get()[secondAFguess]){
            			secondAFguess = i;
            		}
            	}
            }
            double[] normalizedPosteriors = MathUtils.normalizeFromLog10(log10AlleleFrequencyPosteriors.get());
            double logRatio;
            if(Double.isInfinite( Math.log10(normalizedPosteriors[secondAFguess]))){
            	logRatio = MAX_PHRED;
            }
            else{
            	logRatio = 10*(Math.log10(normalizedPosteriors[bestAFguess]) - Math.log10(normalizedPosteriors[secondAFguess]));
            }
            
           // if ( BAC.OutputMode != OUTPUT_MODE.EMIT_ALL_SITES && !passesEmitThreshold(logRatio, bestAFguess) ) {
           // 	double sum = 0.0;
           // 	for (int j = 1; j < log10AlleleFrequencyPosteriors.get().length; j++)       	
          //      	sum += normalizedPosteriors[j];
                    
          //     double PofF = Math.min(sum, 1.0);
          //      return estimateReferenceConfidence(stratifiedContexts, BAC.heterozygosity, true, 1.0 - PofF);
           // }
            
            
            GenotypesContext genotypes = GenotypesContext.create();
            HashMap<String, Object> attributes = new HashMap<String, Object>();

            VariantContext dbsnp = null;

            // search for usable record
            for( final VariantContext vc_input : tracker.getValues(BAC.dbsnp, refContext.getLocus()) ) {
                if ( vc_input != null && vc_input.isSNP() ) {               
                    	dbsnp = vc_input;
                    	break;
                }
            }
            String rsID = null;
            if ( dbsnp != null && dbsnp.hasID() ){
            	rsID = dbsnp.getID();
            //	 attributes.put(BisulfiteVCFConstants.ID_KEY, rsID);
            	 attributes.put(VCFConstants.DBSNP_KEY, true);
            }
               

            // if the site was downsampled, record that fact
            if ( rawContext.hasPileupBeenDownsampled() )
                attributes.put(VCFConstants.DOWNSAMPLED_KEY, true);

            
            
            List<Allele> myAlleles = vc.getAlleles();
 
            GenomeLoc loc = refContext.getLocus();

            int endLoc = calculateEndPos(vc.getAlleles(), vc.getReference(), loc);
            
            assignGenotypes(vc,log10AlleleFrequencyPosteriors.get(),bestAFguess, genotypes, -logRatio/10.0, GL, GLs.keySet());
            
            if ( bestAFguess != 0 ) {
                if(genotypes.containsSample(sample)){
                	if(genotypes.get(sample).isHom() && passesCallThreshold(logRatio)){ //only calculate SB score for Heterozygous SNP and Homozygous loci not pass threshold.. 
                		
                	}
                	else{
                		// the overall lod
                        //double overallLog10PofNull = log10AlleleFrequencyPosteriors.get()[0];
                        double overallLog10PofF = MathUtils.log10sumLog10(log10AlleleFrequencyPosteriors.get(), 1);
                        

                        // the forward lod
                        Map<String, BisulfiteBiallelicGenotypeLikelihoods> tmpGLs = new HashMap<String, BisulfiteBiallelicGenotypeLikelihoods>();
                        VariantContext vcForward = calculateLikelihoods(tracker, refContext, stratifiedContexts,AlignmentContextUtils.ReadOrientation.FORWARD, null, tmpGLs);
                        for ( int i = 0; i < log10AlleleFrequencyPosteriors.get().length; i++ )
                        	log10AlleleFrequencyPosteriors.get()[i] = -1.0 * Double.MAX_VALUE;;
                        
                        if(tmpGLs.containsKey(sample)){
                        	assignAFPosteriors(tmpGLs.get(sample).getLikelihoods(),log10AlleleFrequencyPosteriors.get());
                        }
                        
                        //double[] normalizedLog10Posteriors = MathUtils.normalizeFromLog10(log10AlleleFrequencyPosteriors.get(), true);
                        double forwardLog10PofNull = log10AlleleFrequencyPosteriors.get()[0];
                        double forwardLog10PofF = MathUtils.log10sumLog10(log10AlleleFrequencyPosteriors.get(), 1);
                       

                        // the reverse lod
                        tmpGLs = new HashMap<String, BisulfiteBiallelicGenotypeLikelihoods>();
                        VariantContext vcReverse = calculateLikelihoods(tracker, refContext, stratifiedContexts, AlignmentContextUtils.ReadOrientation.REVERSE, null, tmpGLs);
                        
                        for ( int i = 0; i < log10AlleleFrequencyPosteriors.get().length; i++ )
                        	log10AlleleFrequencyPosteriors.get()[i] = -1.0 * Double.MAX_VALUE;;
                        
                        if(tmpGLs.containsKey(sample)){
                        	assignAFPosteriors(tmpGLs.get(sample).getLikelihoods(),log10AlleleFrequencyPosteriors.get());
                        }
                        
                        //normalizedLog10Posteriors = MathUtils.normalizeFromLog10(log10AlleleFrequencyPosteriors.get(), true);
                        double reverseLog10PofNull = log10AlleleFrequencyPosteriors.get()[0];
                        double reverseLog10PofF = MathUtils.log10sumLog10(log10AlleleFrequencyPosteriors.get(), 1);
                        

                        double forwardLod = forwardLog10PofF + reverseLog10PofNull - overallLog10PofF;
                        double reverseLod = reverseLog10PofF + forwardLog10PofNull - overallLog10PofF;
                        

                        // strand score is max bias between forward and reverse strands
                        double strandScore = Math.max(forwardLod, reverseLod);
                        
                      //  System.err.println(vcForward.getStart() + ": " +  forwardLog10PofF + "\t" + reverseLog10PofNull + "\t" + overallLog10PofF + "\t" + forwardLod + "\t" + reverseLog10PofF + "\t" + forwardLog10PofNull + "\t" + reverseLod + "\t" + strandScore);
                        // rescale by a factor of 10
                        strandScore *= 10.0;
                        //logger.debug(String.format("SLOD=%f", strandScore));

                        attributes.put("SB", Double.valueOf(strandScore));
                	}
                }
  
            }

            double cytosineMethyLevel = 0;
            if(vc.getStart() == BAC.testLocus){
            	System.err.println("hah\t" + vc.getStart() + "\tsampleName: " + sample);
            	
            }
           Integer[] cytosineStat = GL.getCytosineStatus();
          //  if(BAC.ASSUME_SINGLE_SAMPLE != null){
            	 if(passesCallThreshold(logRatio)){
            		// System.err.println(sample + "\t" + genotypes.get(sample) + "\t" + genotypes.size());
            		 Genotype genotypeTemp=genotypes.get(sample);
            		// for(Genotype genotypeTemp : genotypes..values()){
            			 if(genotypeTemp.isHomRef()){
                 			 if(genotypeTemp.getAllele(0).getBases()[0]==BaseUtils.C){
                 				
                                cytosineMethyLevel = (double)cytosineStat[1]/(double)(cytosineStat[1] +cytosineStat[3]);
                                attributes.put(BisulfiteVCFConstants.NUMBER_OF_C_KEY, cytosineStat[1]);
                                attributes.put(BisulfiteVCFConstants.NUMBER_OF_T_KEY, cytosineStat[3]);
                                attributes.put(BisulfiteVCFConstants.C_STRAND_KEY, "+");
                                //attributes.put(BisulfiteVCFConstants.CYTOSINE_METHY_VALUE, cytosineMethyLevel);
                                attributes.put(BisulfiteVCFConstants.CYTOSINE_TYPE, getCytosineTypeStatus(false, cytosineMethyLevel, sample));
                 			 }
                 			 else if(genotypeTemp.getAllele(0).getBases()[0]==BaseUtils.G){
                                cytosineMethyLevel = (double)cytosineStat[0]/(double)(cytosineStat[0] + cytosineStat[2]);
                                attributes.put(BisulfiteVCFConstants.NUMBER_OF_C_KEY, cytosineStat[0]);
                                attributes.put(BisulfiteVCFConstants.NUMBER_OF_T_KEY, cytosineStat[2]);
                                attributes.put(BisulfiteVCFConstants.C_STRAND_KEY, "-");
                               // attributes.put(BisulfiteVCFConstants.CYTOSINE_METHY_VALUE, cytosineMethyLevel);
                                attributes.put(BisulfiteVCFConstants.CYTOSINE_TYPE, getCytosineTypeStatus(true, cytosineMethyLevel, sample));
                 			 }
                 		 }
                 		 else if(genotypeTemp.isHomVar()){
                 			if(genotypeTemp.getAllele(1).getBases()[0]==BaseUtils.C){
                                cytosineMethyLevel = (double)cytosineStat[1]/(double)(cytosineStat[1] + cytosineStat[3]);
                                attributes.put(BisulfiteVCFConstants.NUMBER_OF_C_KEY, cytosineStat[1]);
                                attributes.put(BisulfiteVCFConstants.NUMBER_OF_T_KEY, cytosineStat[3]);
                                attributes.put(BisulfiteVCFConstants.C_STRAND_KEY, "+");
                              //  attributes.put(BisulfiteVCFConstants.CYTOSINE_METHY_VALUE, cytosineMethyLevel);
                                attributes.put(BisulfiteVCFConstants.CYTOSINE_TYPE, getCytosineTypeStatus(false, cytosineMethyLevel, sample));
                 			 }
                 			 else if(genotypeTemp.getAllele(1).getBases()[0]==BaseUtils.G){
                                cytosineMethyLevel = (double)cytosineStat[0]/(double)(cytosineStat[0] + cytosineStat[2]);
                                attributes.put(BisulfiteVCFConstants.NUMBER_OF_C_KEY, cytosineStat[0]);
                                attributes.put(BisulfiteVCFConstants.NUMBER_OF_T_KEY, cytosineStat[2]);
                                attributes.put(BisulfiteVCFConstants.C_STRAND_KEY, "-");
                             //   attributes.put(BisulfiteVCFConstants.CYTOSINE_METHY_VALUE, cytosineMethyLevel);
                                attributes.put(BisulfiteVCFConstants.CYTOSINE_TYPE, getCytosineTypeStatus(true, cytosineMethyLevel, sample));
                 			 }
                 		 }
            			 attributes.put(genotypeTemp.getType().toString(), true);
                 	 }
                // }
          //  } 
            	 
            	 

            	 
            	 VariantContextBuilder vcb = new VariantContextBuilder("BG_call", loc.getContig(), loc.getStart(), endLoc,myAlleles);
            	 vcb.attributes(attributes);
            	 vcb.filters(passesCallThreshold(logRatio) ? null : filter);
            	 vcb.genotypes(genotypes);
            	 if ( rsID != null ){
            		 vcb.id(rsID); 
            	 }
            	 
            	 vcb.log10PError(-logRatio/10.0);
            	 vcb.referenceBaseForIndel(refContext.getBase());
            	 VariantContext vcCall = vcb.make();
            	 System.err.println(sample + "\t" + vcCall.getGenotypes().get(sample) + "\t" + vcCall.getGenotypes().size());
            	 if(call == null){
            		 call = new BisulfiteVariantCallContext(vcCall, rawContext,  passesCallThreshold(logRatio), passesEmitThreshold(logRatio));
            		 call.setRefBase(refContext);
            	 }
            	 call.addCytosineTypeStatus(sample,ctss.get().get(sample));
            	 call.setVariantContext(vcCall);
            //calls.put(sample, value)
            //return call;
        }
        
		return call;
	}

    
	protected boolean passesEmitThreshold(double conf, int bestAFguess) {
        return (BAC.OutputMode != OUTPUT_MODE.EMIT_ALL_SITES || bestAFguess != 0) && conf >= Math.min(BAC.STANDARD_CONFIDENCE_FOR_CALLING, BAC.STANDARD_CONFIDENCE_FOR_EMITTING);
    }
	
		*/
	protected static BisulfiteSNPGenotypeLikelihoodsCalculationModel getGenotypeLikelihoodsCalculationObject(BisulfiteArgumentCollection BAC, boolean autoEstimateC, boolean secondIteration) {		
        	return new BisulfiteSNPGenotypeLikelihoodsCalculationModel(BAC, autoEstimateC, secondIteration);  	
    }

	/*
	protected Map<String, AlignmentContext> getFilteredAndStratifiedContexts(BisulfiteArgumentCollection BAC, ReferenceContext refContext, AlignmentContext rawContext) {
		BadBaseFilterBisulfite badReadPileupFilter = new BadBaseFilterBisulfite(refContext, BAC);

        Map<String, AlignmentContext> stratifiedContexts = null;
       
        if ( !rawContext.hasExtendedEventPileup() ) {

            byte ref = refContext.getBase();
            if ( !BaseUtils.isRegularBase(ref) )
                return null;
            
            stratifiedContexts = AlignmentContextUtils.splitContextBySampleName(rawContext.getBasePileup());
            

            if ( !filterPileupBisulfite(stratifiedContexts, badReadPileupFilter) )
                return null;
        }

        return stratifiedContexts;
    }

	 private boolean filterPileupBisulfite(Map<String, AlignmentContext> stratifiedContexts, BadBaseFilterBisulfite badBaseFilter) {
	        int numDeletions = 0, pileupSize = 0;

	        for ( AlignmentContext context : stratifiedContexts.values() ) {
	            ReadBackedPileup pileup = AlignmentContextUtils.stratify(context,AlignmentContextUtils.ReadOrientation.COMPLETE).getBasePileup();
	            for ( PileupElement p : pileup ) {
	                final SAMRecord read = p.getRead();
	                if(read.getDuplicateReadFlag()){ //get rid of duplicate reads
	                	continue;
	                }
	                if ( p.isDeletion() ) {
	                    // if it's a good read, count it
	                    if ( read.getMappingQuality() >= BAC.MIN_MAPPING_QUALTY_SCORE &&
	                         (BAC.USE_BADLY_MATED_READS || !BadMateFilter.hasBadMate(read)) )
	                        numDeletions++;
	                } else {
	                    if ( !(read instanceof GATKSAMRecord) )
	                        throw new ReviewedStingException("The BisulfiteGenotyper currently expects GATKSAMRecords, but instead saw a " + read.getClass());
	                    GATKSAMRecord GATKrecord = (GATKSAMRecord)read;
	                    GATKSAMRecordFilterStorage GATKrecordFilterStor = new GATKSAMRecordFilterStorage(GATKrecord, badBaseFilter);
	                    if ( GATKrecordFilterStor.isGoodBase(p.getOffset()) )
	                        pileupSize++;
	                }
	            }
	        }

	        if ( BAC.OutputMode == OUTPUT_MODE.EMIT_ALL_SITES )
	            return true;

	        // if no coverage?
	        if ( pileupSize == 0 )
	            return false;

	        // too many deletions in the pileup?
	        if ( (BAC.MAX_DELETION_FRACTION >=0 && BAC.MAX_DELETION_FRACTION <=1.0 ) &&
	                (double)numDeletions / (double)(pileupSize + numDeletions) > BAC.MAX_DELETION_FRACTION )
	            return false;

	        return true;
	    }

	
	protected boolean passesCallThreshold(double conf) {
        return conf >= BAC.STANDARD_CONFIDENCE_FOR_CALLING;
    }
	
	protected boolean passesEmitThreshold(double conf) {
        return conf >= BAC.STANDARD_CONFIDENCE_FOR_EMITTING;
    }
	*/
	private int calculateEndPos(List<Allele> alleles, Allele refAllele, GenomeLoc loc) {
        boolean isSNP = true;
        for (Allele a : alleles){
            if (a.getBaseString().length() != 1) {
                isSNP = false;
                break;
            }
        }

        int endLoc = loc.getStart();
        if ( !isSNP )
            endLoc += refAllele.length();

        return endLoc;
    }
	/*
	private String getCytosineTypeStatus(boolean negStrand, double cytosineMethyLevel, String sample){
		int cPos;
		String cTypeStatus = "C";
		if(negStrand){
			cPos = 1;
		}
		else{
			cPos = 0;
		}
	//	boolean cytosineOnly = false;
		for(String cytosineType : ctss.get().get(sample).cytosineListMap.keySet()){
			String[] tmpKey = cytosineType.split("-");
			Double[] value = ctss.get().get(sample).cytosineListMap.get(cytosineType);

			if(Double.compare(value[3], 1.0) == 0){ //in first iteration, it will require higher confidance for calling C, 100 times more than the other type of C; then second iteration, it just 10 times more likelihood than any other type of C
				if(tmpKey[0].equalsIgnoreCase("C")){
					//cTypeStatus = tmpKey[0];
			//		cytosineOnly = true;
				}
				else{
					cTypeStatus = cTypeStatus + "," + tmpKey[0];
				}
				
				 if(Double.isNaN(cytosineMethyLevel))
						 cytosineMethyLevel = 0.0;
				 value[2] = cytosineMethyLevel;
				
				if(tmpKey[0].equalsIgnoreCase("C")){
					ctss.get().get(sample).isC = true;
					ctss.get().get(sample).cytosineMethyLevel = cytosineMethyLevel;
				
				}
				else if(tmpKey[0].equalsIgnoreCase("CG")){
					ctss.get().get(sample).isCpg = true;
					ctss.get().get(sample).cpgMethyLevel = cytosineMethyLevel;
					ctss.get().get(sample).isC = true;
					ctss.get().get(sample).cytosineMethyLevel = cytosineMethyLevel;
				}
				else if(tmpKey[0].equalsIgnoreCase("CH")){
					ctss.get().get(sample).isCph = true;
					ctss.get().get(sample).cphMethyLevel = cytosineMethyLevel;
					ctss.get().get(sample).isC = true;
					ctss.get().get(sample).cytosineMethyLevel = cytosineMethyLevel;
				}
				//else if(tmpKey[0].equalsIgnoreCase("CHH")){
				//	ctss.get().isChh = true;
				//	ctss.get().chhMethyLevel = cytosineMethyLevel;
				//}
				//else if(tmpKey[0].equalsIgnoreCase("CHG")){
				//	ctss.get().isChg = true;
				//	ctss.get().chgMethyLevel = cytosineMethyLevel;
				//}
				if(tmpKey[0].equalsIgnoreCase("GCH")){
					if(BAC.sequencingMode == MethylSNPModel.GM){
						ctss.get().get(sample).isGch = true;
						ctss.get().get(sample).gchMethyLevel = cytosineMethyLevel;
						ctss.get().get(sample).isC = true;
						ctss.get().get(sample).cytosineMethyLevel = cytosineMethyLevel;
					}
					else{
						continue;
					}
					
				}
				else if(tmpKey[0].equalsIgnoreCase("CCH")){
					if(BAC.sequencingMode == MethylSNPModel.GM){
						ctss.get().get(sample).isCch = true;
						ctss.get().get(sample).cchMethyLevel = cytosineMethyLevel;
						ctss.get().get(sample).isC = true;
						ctss.get().get(sample).cytosineMethyLevel = cytosineMethyLevel;
					}
					else{
						continue;
					}
					
				}
				else if(tmpKey[0].equalsIgnoreCase("WCH")){
					if(BAC.sequencingMode == MethylSNPModel.GM){
						ctss.get().get(sample).isWch = true;
						ctss.get().get(sample).wchMethyLevel = cytosineMethyLevel;
						ctss.get().get(sample).isC = true;
						ctss.get().get(sample).cytosineMethyLevel = cytosineMethyLevel;
					}
					else{
						continue;
					}
					
				}
				else if(tmpKey[0].equalsIgnoreCase("GCG")){
					if(BAC.sequencingMode == MethylSNPModel.GM){
						ctss.get().get(sample).isGcg = true;
						ctss.get().get(sample).gcgMethyLevel = cytosineMethyLevel;
						ctss.get().get(sample).isC = true;
						ctss.get().get(sample).cytosineMethyLevel = cytosineMethyLevel;
					}
					else{
						continue;
					}
					
				}
				else if(tmpKey[0].equalsIgnoreCase("CCG")){
					if(BAC.sequencingMode == MethylSNPModel.GM){
						ctss.get().get(sample).isCcg = true;
						ctss.get().get(sample).ccgMethyLevel = cytosineMethyLevel;
						ctss.get().get(sample).isC = true;
						ctss.get().get(sample).cytosineMethyLevel = cytosineMethyLevel;
					}
					else{
						continue;
					}
					
				}
				else if(tmpKey[0].equalsIgnoreCase("WCG")){
					if(BAC.sequencingMode == MethylSNPModel.GM){
						ctss.get().get(sample).isWcg = true;
						ctss.get().get(sample).wcgMethyLevel = cytosineMethyLevel;
						ctss.get().get(sample).isC = true;
						ctss.get().get(sample).cytosineMethyLevel = cytosineMethyLevel;
					}
					else{
						continue;
					}
					
				}
				ctss.get().get(sample).cytosineListMap.put(cytosineType,value);
				
			}
			
 		}
	//	if(cytosineOnly){
			return cTypeStatus;
	//	}
	//	else{
	//		cTypeStatus="";
	//		return cTypeStatus;
	//	}
		
	}
*/
	
	   /**
     * Can be overridden by concrete subclasses
     * @param vc                   variant context with genotype likelihoods
     * @param log10AlleleFrequencyPosteriors    allele frequency results
     * @param AFofMaxLikelihood    allele frequency of max likelihood
     *
     * @return calls
     */
	/*
    public void assignGenotypes(VariantContext vc,
                                                 double[] log10AlleleFrequencyPosteriors,
                                                 int bestAF,
                                                 GenotypesContext calls,
                                                 double logRatio,
                                                 BisulfiteBiallelicGenotypeLikelihoods BGL,
                                                 Set<String> sampleNames) {
        if ( !vc.isVariant() )
            throw new UserException("The VCF record passed in does not contain an ALT allele at " + vc.getChr() + ":" + vc.getStart());

        	GenotypesContext GLs = vc.getGenotypes();
        //	Set<String> sampleNames = GLs.getSampleNames();
        	for(String sampleName : sampleNames){
        		Genotype g = GLs.get(sampleName);
            	
             //   if ( !g.hasLikelihoods() )
             //       continue;
                
                Allele alleleA = BGL.getAlleleA();
                Allele alleleB = BGL.getAlleleB();
                ArrayList<Allele> myAlleles = new ArrayList<Allele>();
                if (bestAF == 0) {
                    myAlleles.add(alleleA);
                    myAlleles.add(alleleA);
                    
                } else if(bestAF == 1) {
                    myAlleles.add(alleleA);
                    myAlleles.add(alleleB);
                   

                }  else {
                    myAlleles.add(alleleB);
                    myAlleles.add(alleleB);
                    
                }
               // System.err.println(sampleName + "\t" + myAlleles.toString() + "\t" + logRatio + "\t" + GLs.getSampleNames());
                //calls.put(sampleName, new Genotype(sampleName, myAlleles, logRatio, null, g.getAttributes(), false));
                calls.add(new Genotype(sampleName, myAlleles, logRatio, null, g.getAttributes(), false));
        	}
        	
        		
    }
	*/
    
    private VariantContext createVariantContextFromLikelihoods(ReferenceContext refContext, Allele refAllele, Map<String, BisulfiteContextsGenotypeLikelihoods> BCGLs) {
        
        List<Allele> noCall = new ArrayList<Allele>();
        noCall.add(Allele.NO_CALL);

        List<Allele> alleles = new ArrayList<Allele>();
        alleles.add(refAllele);
        boolean addedAltAllele = false;

        GenotypesContext genotypes = GenotypesContext.create();;
        for ( BisulfiteContextsGenotypeLikelihoods BCGL : BCGLs.values() ) {
         //   if ( !addedAltAllele ) {
         //       addedAltAllele = true;  
                if(!alleles.contains(BCGL.getAlleleA()))
                	alleles.add(BCGL.getAlleleA());
                if(!alleles.contains(BCGL.getAlleleB()))
                	alleles.add(BCGL.getAlleleB());
          //  }   
                HashMap<String, Object> attributes = new HashMap<String, Object>();
                GenotypeLikelihoods likelihoods = GenotypeLikelihoods.fromLog10Likelihoods(BCGL.getLikelihoods());
                attributes.put(VCFConstants.DEPTH_KEY, BCGL.getDepth());
                attributes.put(VCFConstants.GENOTYPE_POSTERIORS_KEY, likelihoods.getAsString());
                attributes.put(BisulfiteVCFConstants.NUMBER_OF_C_KEY, BCGL.getNumOfCReadsInBisulfiteCStrand());
                attributes.put(BisulfiteVCFConstants.NUMBER_OF_T_KEY, BCGL.getNumOfTReadsInBisulfiteCStrand()); //need to fix it to judge strand
                attributes.put(BisulfiteVCFConstants.BEST_C_PATTERN, BCGL.getBestMatchedCytosinePattern());
                genotypes.add(new Genotype(BCGL.getSample(), alleles, Genotype.NO_LOG10_PERROR, null, attributes, false));

        }

        GenomeLoc loc = refContext.getLocus();
        int endLoc = calculateEndPos(alleles, refAllele, loc);
        
        VariantContextBuilder vcb = new VariantContextBuilder("BG_call", loc.getContig(), loc.getStart(), endLoc,alleles);
        vcb.referenceBaseForIndel(refContext.getBase());
   	 	vcb.log10PError(VariantContext.NO_LOG10_PERROR);
   	 	vcb.genotypes(genotypes);
        return vcb.make();
    }
    
    
    
    
}
