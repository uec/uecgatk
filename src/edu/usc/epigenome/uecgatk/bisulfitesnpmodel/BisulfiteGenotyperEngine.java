package edu.usc.epigenome.uecgatk.bisulfitesnpmodel;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.StringUtil;

import org.apache.log4j.Logger;
import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.GenotypeLikelihoods;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFConstants;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext.StratifiedContextType;
import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.helpers.DbSNPHelper;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.AlleleFrequencyCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.BiallelicGenotypeLikelihoods;
import org.broadinstitute.sting.gatk.walkers.genotyper.DindelGenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.DiploidIndelGenotypePriors;
import org.broadinstitute.sting.gatk.walkers.genotyper.DiploidSNPGenotypePriors;
import org.broadinstitute.sting.gatk.walkers.genotyper.ExactAFCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypePriors;
import org.broadinstitute.sting.gatk.walkers.genotyper.GridSearchAFEstimation;
import org.broadinstitute.sting.gatk.walkers.genotyper.SNPGenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine.OUTPUT_MODE;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedExtendedEventPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.GATKSAMRecordFilter;

import edu.usc.epigenome.uecgatk.bisulfitesnpmodel.BisulfiteAlignmentUtils;
import edu.usc.epigenome.uecgatk.bisulfitesnpmodel.BisulfiteSNPGenotypeLikelihoodsCalculationModel;
import edu.usc.epigenome.uecgatk.bisulfitesnpmodel.BisulfiteVCFConstants;
import edu.usc.epigenome.uecgatk.bisulfitesnpmodel.NonRefDependSNPGenotypeLikelihoodsCalculationModel.MethylSNPModel;

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
	
	//cytosine pattern and their status
	private ThreadLocal<CytosineTypeStatus> ctss = new ThreadLocal<CytosineTypeStatus>();

    // the model used for calculating genotypes
    private ThreadLocal<BisulfiteSNPGenotypeLikelihoodsCalculationModel> bglcms = new ThreadLocal<BisulfiteSNPGenotypeLikelihoodsCalculationModel>();

    // the allele frequency likelihoods (allocated once as an optimization)
    private ThreadLocal<double[]> log10AlleleFrequencyPosteriors = new ThreadLocal<double[]>();

    // the priors object
 //   private GenotypePriors genotypePriors;

    // samples in input
    private Set<String> samples = new TreeSet<String>();

    // the various loggers and writers
    private Logger logger = null;

    // fasta reference reader to supplement the edges of the reference sequence for long reads
    private IndexedFastaSequenceFile referenceReader;

    // the standard filter to use for calls below the confidence threshold but above the emit threshold
    private static final Set<String> filter = new HashSet<String>(1);

    private static final int MISMATCH_WINDOW_SIZE = 20;
    
    private static boolean autoEstimateC = false;
    private static boolean secondIteration = false;

	protected double MAX_PHRED = 1000000;
	
	public static final String LOW_QUAL_FILTER_NAME = "LowQual";

	public BisulfiteGenotyperEngine(GenomeAnalysisEngine toolkit,
			BisulfiteArgumentCollection BAC, Logger logger,
			Set<String> samples) {
		this.samples = new TreeSet<String>(samples);
		initialize(toolkit, BAC, logger, samples.size());
		// TODO Auto-generated constructor stub
	}
	

	protected void initialize(GenomeAnalysisEngine toolkit, BisulfiteArgumentCollection BAC, Logger logger, int numSamples) {
        this.BAC = BAC.clone();
        this.logger = logger;
     //   genotypePriors = BisulfiteGenotyperEngine.createGenotypePriors(BAC);
        filter.add(LOW_QUAL_FILTER_NAME);

        try {
            referenceReader = new CachingIndexedFastaSequenceFile(toolkit.getArguments().referenceFile);
        }
        catch(FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile(toolkit.getArguments().referenceFile,ex);
        }
    }
	
	public void setCytosineTypeStatus(CytosineTypeStatus cts){
		this.ctss.set(cts.clone());
	}
	
	/**
     * Compute full BisulfiteVariantCallContext at a given locus.
     *
     * @param tracker    the meta data tracker
     * @param refContext the reference base
     * @param rawContext contextual information around the locus
     * @return the VariantCallContext object
     */
    public BisulfiteVariantCallContext calculateLikelihoodsAndGenotypes(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext) {
        Map<String, StratifiedAlignmentContext> stratifiedContexts = getFilteredAndStratifiedContexts(BAC, refContext, rawContext);
        Map<String, BisulfiteBiallelicGenotypeLikelihoods> GLs = new HashMap<String, BisulfiteBiallelicGenotypeLikelihoods>();
        VariantContext vc = calculateLikelihoods(tracker, refContext, stratifiedContexts, StratifiedAlignmentContext.StratifiedContextType.COMPLETE, null, GLs);
        
        if ( vc == null )
            return null;

        BisulfiteVariantCallContext vcc = calculateGenotypes(tracker, refContext, rawContext, stratifiedContexts, GLs, vc);
        
        return vcc;
    }
	//
	//protected static GenotypePriors createGenotypePriors(BisulfiteArgumentCollection BAC) {
   //     return new BisulfiteDiploidSNPGenotypePriors();
   // }
	


    protected VariantContext calculateLikelihoods(RefMetaDataTracker tracker, ReferenceContext refContext, Map<String, StratifiedAlignmentContext> stratifiedContexts, StratifiedAlignmentContext.StratifiedContextType type, Allele alternateAlleleToUse, Map<String, BisulfiteBiallelicGenotypeLikelihoods> GLs) {
		if ( stratifiedContexts == null ){
			 return null;
		}

        // initialize the data for this thread if that hasn't been done yet
        if ( bglcms.get() == null ) {
            bglcms.set(BisulfiteGenotyperEngine.getGenotypeLikelihoodsCalculationObject(logger, BAC));
            
        }

        BisulfiteSNPGenotypeLikelihoodsCalculationModel bglcm = (BisulfiteSNPGenotypeLikelihoodsCalculationModel) bglcms.get();
        bglcm.initialize(ctss.get(), BAC, autoEstimateC, secondIteration);
        
        Allele refAllele = bglcm.getBsLikelihoods(tracker, refContext, stratifiedContexts, type, GLs, alternateAlleleToUse);
       
        if (refAllele != null)
            return createVariantContextFromLikelihoods(refContext, refAllele, GLs);
        else
            return null;
          
        
    }
    
    protected void assignAFPosteriors(double[]likelihoods, double[] log10AFPosteriors){
    	for(int i = 0; i < likelihoods.length; i++){
    		log10AFPosteriors[i] = likelihoods[i];
    	}
    		
    }
	

	protected BisulfiteVariantCallContext calculateGenotypes(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext, Map<String, StratifiedAlignmentContext> stratifiedContexts, Map<String, BisulfiteBiallelicGenotypeLikelihoods> GLs, VariantContext vc) {
		// initialize the data for this thread if that hasn't been done yet
        if ( bglcms.get() == null ) {
            return null;
        }
        log10AlleleFrequencyPosteriors.set(new double[3]);
        for ( String sample : GLs.keySet() ) {
        	BisulfiteBiallelicGenotypeLikelihoods GL = GLs.get(sample);
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
            
            if ( BAC.OutputMode != OUTPUT_MODE.EMIT_ALL_SITES && !passesEmitThreshold(logRatio, bestAFguess) ) {
            	double sum = 0.0;
            	for (int j = 1; j < log10AlleleFrequencyPosteriors.get().length; j++)       	
                	sum += normalizedPosteriors[j];
                    
               double PofF = Math.min(sum, 1.0);
                return estimateReferenceConfidence(stratifiedContexts, BAC.heterozygosity, true, 1.0 - PofF);
            }
            
            HashMap<String, Genotype> genotypes = new HashMap<String, Genotype>();
            
            HashMap<String, Object> attributes = new HashMap<String, Object>();

            String rsID = DbSNPHelper.rsIDOfFirstRealSNP(tracker.getReferenceMetaData(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME));
            if ( rsID != null ){
            	 attributes.put(VariantContext.ID_KEY, rsID);
            	 attributes.put(VCFConstants.DBSNP_KEY, true);
            }
               

            // if the site was downsampled, record that fact
            if ( rawContext.hasPileupBeenDownsampled() )
                attributes.put(VCFConstants.DOWNSAMPLED_KEY, true);

            
            
            Set<Allele> myAlleles = vc.getAlleles();
 
            GenomeLoc loc = refContext.getLocus();

            int endLoc = calculateEndPos(vc.getAlleles(), vc.getReference(), loc);
            
            assignGenotypes(vc,log10AlleleFrequencyPosteriors.get(),bestAFguess, genotypes, logRatio/10.0, GL);
            
            if ( bestAFguess != 0 ) {
                if(genotypes.containsKey(sample)){
                	if(genotypes.get(sample).isHom() && passesCallThreshold(logRatio)){ //only calculate SB score for Heterozygous SNP and Homozygous loci not pass threshold.. 
                		
                	}
                	else{
                		// the overall lod
                        //double overallLog10PofNull = log10AlleleFrequencyPosteriors.get()[0];
                        double overallLog10PofF = MathUtils.log10sumLog10(log10AlleleFrequencyPosteriors.get(), 1);
                        

                        // the forward lod
                        Map<String, BisulfiteBiallelicGenotypeLikelihoods> tmpGLs = new HashMap<String, BisulfiteBiallelicGenotypeLikelihoods>();
                        VariantContext vcForward = calculateLikelihoods(tracker, refContext, stratifiedContexts, StratifiedAlignmentContext.StratifiedContextType.FORWARD, null, tmpGLs);
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
                        VariantContext vcReverse = calculateLikelihoods(tracker, refContext, stratifiedContexts, StratifiedAlignmentContext.StratifiedContextType.REVERSE, null, tmpGLs);
                        
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
            
           Integer[] cytosineStat = GL.getCytosineStatus();
            if(BAC.ASSUME_SINGLE_SAMPLE != null){
            	 if(passesCallThreshold(logRatio)){
            		 
            		 for(Genotype genotypeTemp : genotypes.values()){
            			 if(genotypeTemp.isHomRef()){
                 			 if(genotypeTemp.getAllele(0).getBases()[0]==BaseUtils.C){
                 				
                                cytosineMethyLevel = (double)cytosineStat[1]/(double)(cytosineStat[1] +cytosineStat[3]);
                                attributes.put(BisulfiteVCFConstants.NUMBER_OF_C_KEY, cytosineStat[1]);
                                attributes.put(BisulfiteVCFConstants.NUMBER_OF_T_KEY, cytosineStat[3]);
                                attributes.put(BisulfiteVCFConstants.C_IN_NEG_STRAND_KEY, false);
                                //attributes.put(BisulfiteVCFConstants.CYTOSINE_METHY_VALUE, cytosineMethyLevel);
                                attributes.put(BisulfiteVCFConstants.CYTOSINE_TYPE, getCytosineTypeStatus(false, cytosineMethyLevel));
                 			 }
                 			 else if(genotypeTemp.getAllele(0).getBases()[0]==BaseUtils.G){
                                cytosineMethyLevel = (double)cytosineStat[0]/(double)(cytosineStat[0] + cytosineStat[2]);
                                attributes.put(BisulfiteVCFConstants.NUMBER_OF_C_KEY, cytosineStat[0]);
                                attributes.put(BisulfiteVCFConstants.NUMBER_OF_T_KEY, cytosineStat[2]);
                                attributes.put(BisulfiteVCFConstants.C_IN_NEG_STRAND_KEY, true);
                               // attributes.put(BisulfiteVCFConstants.CYTOSINE_METHY_VALUE, cytosineMethyLevel);
                                attributes.put(BisulfiteVCFConstants.CYTOSINE_TYPE, getCytosineTypeStatus(true, cytosineMethyLevel));
                 			 }
                 		 }
                 		 else if(genotypeTemp.isHomVar()){
                 			if(genotypeTemp.getAllele(1).getBases()[0]==BaseUtils.C){
                                cytosineMethyLevel = (double)cytosineStat[1]/(double)(cytosineStat[1] + cytosineStat[3]);
                                attributes.put(BisulfiteVCFConstants.NUMBER_OF_C_KEY, cytosineStat[1]);
                                attributes.put(BisulfiteVCFConstants.NUMBER_OF_T_KEY, cytosineStat[3]);
                                attributes.put(BisulfiteVCFConstants.C_IN_NEG_STRAND_KEY, false);
                              //  attributes.put(BisulfiteVCFConstants.CYTOSINE_METHY_VALUE, cytosineMethyLevel);
                                attributes.put(BisulfiteVCFConstants.CYTOSINE_TYPE, getCytosineTypeStatus(false, cytosineMethyLevel));
                 			 }
                 			 else if(genotypeTemp.getAllele(1).getBases()[0]==BaseUtils.G){
                                cytosineMethyLevel = (double)cytosineStat[0]/(double)(cytosineStat[0] + cytosineStat[2]);
                                attributes.put(BisulfiteVCFConstants.NUMBER_OF_C_KEY, cytosineStat[0]);
                                attributes.put(BisulfiteVCFConstants.NUMBER_OF_T_KEY, cytosineStat[2]);
                                attributes.put(BisulfiteVCFConstants.C_IN_NEG_STRAND_KEY, true);
                             //   attributes.put(BisulfiteVCFConstants.CYTOSINE_METHY_VALUE, cytosineMethyLevel);
                                attributes.put(BisulfiteVCFConstants.CYTOSINE_TYPE, getCytosineTypeStatus(true, cytosineMethyLevel));
                 			 }
                 		 }
            			 attributes.put(genotypeTemp.getType().toString(), true);
                 	 }
                 }
            } 
            VariantContext vcCall = new VariantContext("BG_call", loc.getContig(), loc.getStart(), endLoc,
                    myAlleles, genotypes, logRatio/10.0, passesCallThreshold(logRatio) ? null : filter, attributes);
            BisulfiteVariantCallContext call = new BisulfiteVariantCallContext(vcCall, rawContext,  passesCallThreshold(logRatio), ctss.get(), passesEmitThreshold(logRatio));
            call.setRefBase(refContext.getBase());
            return call;
        }
		return null;
	}
	
	protected boolean passesEmitThreshold(double conf, int bestAFguess) {
        return (BAC.OutputMode != OUTPUT_MODE.EMIT_ALL_SITES || bestAFguess != 0) && conf >= Math.min(BAC.STANDARD_CONFIDENCE_FOR_CALLING, BAC.STANDARD_CONFIDENCE_FOR_EMITTING);
    }
	
	
	protected static BisulfiteSNPGenotypeLikelihoodsCalculationModel getGenotypeLikelihoodsCalculationObject(Logger logger, BisulfiteArgumentCollection BAC) {		
        	return new BisulfiteSNPGenotypeLikelihoodsCalculationModel(BAC, logger);  	
    }

	protected Map<String, StratifiedAlignmentContext> getFilteredAndStratifiedContexts(BisulfiteArgumentCollection BAC, ReferenceContext refContext, AlignmentContext rawContext) {
		BadBaseFilterBisulfite badReadPileupFilter = new BadBaseFilterBisulfite(refContext, BAC);

        Map<String, StratifiedAlignmentContext> stratifiedContexts = null;
       
        if ( !rawContext.hasExtendedEventPileup() ) {

            byte ref = refContext.getBase();
            if ( !BaseUtils.isRegularBase(ref) )
                return null;
            
            stratifiedContexts = StratifiedAlignmentContext.splitContextBySampleName(rawContext.getBasePileup(), BAC.ASSUME_SINGLE_SAMPLE);
            

            if ( !filterPileupBisulfite(stratifiedContexts, badReadPileupFilter) )
                return null;
        }

        return stratifiedContexts;
    }

	 private boolean filterPileupBisulfite(Map<String, StratifiedAlignmentContext> stratifiedContexts, BisulfiteGenotyperEngine.BadBaseFilterBisulfite badBaseFilter) {
	        int numDeletions = 0, pileupSize = 0;

	        for ( StratifiedAlignmentContext context : stratifiedContexts.values() ) {
	            ReadBackedPileup pileup = context.getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getBasePileup();
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
	                    GATKrecord.setGoodBases(badBaseFilter, true);
	                    if ( GATKrecord.isGoodBase(p.getOffset()) )
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

	 //copy from GATK, since it is not public class there
	protected class BadBaseFilterBisulfite implements GATKSAMRecordFilter {
        private ReferenceContext refContext;
        private final BisulfiteArgumentCollection BAC;

        public BadBaseFilterBisulfite(ReferenceContext refContext, BisulfiteArgumentCollection BAC) {
            this.refContext = refContext;
            this.BAC = BAC;
        }

        @Override
        public BitSet getGoodBases(final GATKSAMRecord record) {
            BitSet bitset = new BitSet(record.getReadLength());

            // if the mapping quality is too low or the mate is bad, we can just zero out the whole read and continue
            if ( record.getMappingQuality() < BAC.MIN_MAPPING_QUALTY_SCORE ||
                 (!BAC.USE_BADLY_MATED_READS && BadMateFilter.hasBadMate(record)) ) {
            	//System.out.println("bad mates");
            	return bitset;
            }
            
            byte[] quals = record.getBaseQualities();
            for (int i = 0; i < quals.length; i++) {
                if ( quals[i] >= BAC.MIN_BASE_QUALTY_SCORE )
                    bitset.set(i);
            }

            // if a read is too long for the reference context, extend the context (being sure not to extend past the end of the chromosome)
            if ( record.getAlignmentEnd() > refContext.getWindow().getStop() ) {
                GenomeLoc window = refContext.getGenomeLocParser().createGenomeLoc(refContext.getLocus().getContig(), refContext.getWindow().getStart(), Math.min(record.getAlignmentEnd(), referenceReader.getSequenceDictionary().getSequence(refContext.getLocus().getContig()).getSequenceLength()));
                byte[] bases = referenceReader.getSubsequenceAt(window.getContig(), window.getStart(), window.getStop()).getBases();
                StringUtil.toUpperCase(bases);
                refContext = new ReferenceContext(refContext.getGenomeLocParser(),refContext.getLocus(), window, bases);
            }            

            BitSet mismatches;
            if(BAC.sequencingMode == MethylSNPModel.BM || BAC.sequencingMode == MethylSNPModel.GM){
            	mismatches = BisulfiteAlignmentUtils.mismatchesInRefWindow(record, refContext, BAC.MAX_MISMATCHES, MISMATCH_WINDOW_SIZE, BAC.sequencingMode, BAC.pairedEndMode);
            }
            else{
            	mismatches = AlignmentUtils.mismatchesInRefWindow(record, refContext, BAC.MAX_MISMATCHES, MISMATCH_WINDOW_SIZE);
            }
            if ( mismatches != null )
                bitset.and(mismatches);

            return bitset;
        }
    }
	
	//copy from GATK since private class issue
	private BisulfiteVariantCallContext estimateReferenceConfidence(Map<String, StratifiedAlignmentContext> contexts, double theta, boolean ignoreCoveredSamples, double initialPofRef) {
        if ( contexts == null )
            return null;

        double P_of_ref = initialPofRef;

        // for each sample that we haven't examined yet
        for ( String sample : samples ) {
            boolean isCovered = contexts.containsKey(sample);
            if ( ignoreCoveredSamples && isCovered )
                continue;


            int depth = 0;

            if (isCovered) {
                AlignmentContext context =  contexts.get(sample).getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE);

                if (context.hasBasePileup())
                    depth = context.getBasePileup().size();
                else if (context.hasExtendedEventPileup())
                    depth = context.getExtendedEventPileup().size();
            }

            P_of_ref *= 1.0 - (theta / 2.0) * MathUtils.binomialProbability(0, depth, 0.5);
        }

        return new BisulfiteVariantCallContext(QualityUtils.phredScaleErrorRate(1.0 - P_of_ref) >= BAC.STANDARD_CONFIDENCE_FOR_CALLING, ctss.get(),QualityUtils.phredScaleErrorRate(1.0 - P_of_ref) >= BAC.STANDARD_CONFIDENCE_FOR_EMITTING);
    }
	
	protected boolean passesCallThreshold(double conf) {
        return conf >= BAC.STANDARD_CONFIDENCE_FOR_CALLING;
    }
	
	protected boolean passesEmitThreshold(double conf) {
        return conf >= BAC.STANDARD_CONFIDENCE_FOR_EMITTING;
    }
	
	private int calculateEndPos(Set<Allele> alleles, Allele refAllele, GenomeLoc loc) {
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
	
	private String getCytosineTypeStatus(boolean negStrand, double cytosineMethyLevel){
		int cPos;
		String cTypeStatus = "C";
		if(negStrand){
			cPos = 1;
		}
		else{
			cPos = 0;
		}
		boolean first = true;
		for(String cytosineType : ctss.get().cytosineListMap.keySet()){
			String[] tmpKey = cytosineType.split("-");
			Double[] value = ctss.get().cytosineListMap.get(cytosineType);

			if(Double.compare(value[3], 1.0) == 0){ //in first iteration, it will require higher confidance for calling C, 100 times more than the other type of C; then second iteration, it just 10 times more likelihood than any other type of C
				if(first){
					cTypeStatus = tmpKey[0];
					first = false;
				}
				else{
					cTypeStatus = cTypeStatus + "," + tmpKey[0];
				}
				
				 if(Double.isNaN(cytosineMethyLevel))
						 cytosineMethyLevel = 0.0;
				 value[2] = cytosineMethyLevel;
				
				if(tmpKey[0].equalsIgnoreCase("C")){
					ctss.get().isC = true;
					ctss.get().cytosineMethyLevel = cytosineMethyLevel;
				
				}
				else if(tmpKey[0].equalsIgnoreCase("CG")){
					ctss.get().isCpg = true;
					ctss.get().cpgMethyLevel = cytosineMethyLevel;
				}
				else if(tmpKey[0].equalsIgnoreCase("CH")){
					ctss.get().isCph = true;
					ctss.get().cphMethyLevel = cytosineMethyLevel;
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
						ctss.get().isGch = true;
						ctss.get().gchMethyLevel = cytosineMethyLevel;
					}
					else{
						continue;
					}
					
				}
				else if(tmpKey[0].equalsIgnoreCase("HCH")){
					if(BAC.sequencingMode == MethylSNPModel.GM){
						ctss.get().isHch = true;
						ctss.get().hchMethyLevel = cytosineMethyLevel;
					}
					else{
						continue;
					}
					
				}
				else if(tmpKey[0].equalsIgnoreCase("WCH")){
					if(BAC.sequencingMode == MethylSNPModel.GM){
						ctss.get().isWch = true;
						ctss.get().wchMethyLevel = cytosineMethyLevel;
					}
					else{
						continue;
					}
					
				}
				else if(tmpKey[0].equalsIgnoreCase("GCG")){
					if(BAC.sequencingMode == MethylSNPModel.GM){
						ctss.get().isGcg = true;
						ctss.get().gcgMethyLevel = cytosineMethyLevel;
					}
					else{
						continue;
					}
					
				}
				else if(tmpKey[0].equalsIgnoreCase("HCG")){
					if(BAC.sequencingMode == MethylSNPModel.GM){
						ctss.get().isHcg = true;
						ctss.get().hcgMethyLevel = cytosineMethyLevel;
					}
					else{
						continue;
					}
					
				}
				else if(tmpKey[0].equalsIgnoreCase("WCG")){
					if(BAC.sequencingMode == MethylSNPModel.GM){
						ctss.get().isWcg = true;
						ctss.get().wcgMethyLevel = cytosineMethyLevel;
					}
					else{
						continue;
					}
					
				}
				ctss.get().cytosineListMap.put(cytosineType,value);
				
			}
			
 		}
		return cTypeStatus;
	}

	   /**
     * Can be overridden by concrete subclasses
     * @param vc                   variant context with genotype likelihoods
     * @param log10AlleleFrequencyPosteriors    allele frequency results
     * @param AFofMaxLikelihood    allele frequency of max likelihood
     *
     * @return calls
     */
    public void assignGenotypes(VariantContext vc,
                                                 double[] log10AlleleFrequencyPosteriors,
                                                 int bestAF,
                                                 HashMap<String, Genotype> calls,
                                                 double logRatio,
                                                 BiallelicGenotypeLikelihoods BGL) {
        if ( !vc.isVariant() )
            throw new UserException("The VCF record passed in does not contain an ALT allele at " + vc.getChr() + ":" + vc.getStart());

        	Map<String, Genotype> GLs = vc.getGenotypes();
              for ( Map.Entry<String, Genotype> sample : GLs.entrySet() ) {
                if ( !sample.getValue().hasLikelihoods() )
                    continue;
                Genotype g = sample.getValue();
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
                calls.put(sample.getKey(), new Genotype(sample.getKey(), myAlleles, logRatio, null, g.getAttributes(), false));
              }
    }
	
    private VariantContext createVariantContextFromLikelihoods(ReferenceContext refContext, Allele refAllele, Map<String, BisulfiteBiallelicGenotypeLikelihoods> GLs) {
        
        List<Allele> noCall = new ArrayList<Allele>();
        noCall.add(Allele.NO_CALL);

        Set<Allele> alleles = new HashSet<Allele>();
        alleles.add(refAllele);
        boolean addedAltAllele = false;

        HashMap<String, Genotype> genotypes = new HashMap<String, Genotype>();
        for ( BisulfiteBiallelicGenotypeLikelihoods GL : GLs.values() ) {
            if ( !addedAltAllele ) {
                addedAltAllele = true;              
                alleles.add(GL.getAlleleA());
                alleles.add(GL.getAlleleB());
                
                HashMap<String, Object> attributes = new HashMap<String, Object>();
                GenotypeLikelihoods likelihoods = new GenotypeLikelihoods(GL.getLikelihoods());
                attributes.put(VCFConstants.DEPTH_KEY, GL.getDepth());
                attributes.put(VCFConstants.GENOTYPE_LIKELIHOODS_KEY, likelihoods.getAsString());
                genotypes.put(GL.getSample(), new Genotype(GL.getSample(), noCall, Genotype.NO_NEG_LOG_10PERROR, null, attributes, false));
            }

        }

        GenomeLoc loc = refContext.getLocus();
        int endLoc = calculateEndPos(alleles, refAllele, loc);
        
        return new VariantContext("BG_call",
                loc.getContig(),
                loc.getStart(),
                endLoc,
                alleles,
                genotypes,
                VariantContext.NO_NEG_LOG_10PERROR,
                null,
                null);
    }
    
    public void setAutoParameters(boolean autoEstimateC, boolean secondIteration){
    	this.autoEstimateC = autoEstimateC;
    	this.secondIteration = secondIteration;
    }
    
    public enum OUTPUT_MODE {
        EMIT_VARIANTS_ONLY, //only confident variants
        EMIT_ALL_CONFIDENT_SITES,
        EMIT_ALL_SITES,
        EMIT_ALL_CPG, //only confident cpgs
        EMIT_ALL_CYTOSINES, //only confident cytosines
        EMIT_HET_SNPS_ONLY, //only confident heterozygous snps
        DEFAULT_FOR_TCGA //output two vcf files: 1. confident variants vcf file (not contain cytosine info), 2. confident cpgs vcf file (only contain homozygous cpg info), 3. in future, output cpg_read files
    }
    
}
