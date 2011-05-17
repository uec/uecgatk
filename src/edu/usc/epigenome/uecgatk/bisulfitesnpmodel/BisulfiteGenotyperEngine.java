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

import org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel.BisulfiteAlignmentUtils;
import org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel.BisulfiteSNPGenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel.BisulfiteVCFConstants;
import org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel.NonRefDependSNPGenotypeLikelihoodsCalculationModel.MethylSNPModel;

public class BisulfiteGenotyperEngine{

	private BisulfiteArgumentCollection BAC = null;
	
	private CytosineTypeStatus cts = null;
	
	// the annotation engine
    private VariantAnnotatorEngine annotationEngine;
	
    
    
    // the model used for calculating genotypes
    private ThreadLocal<GenotypeLikelihoodsCalculationModel> glcm = new ThreadLocal<GenotypeLikelihoodsCalculationModel>();

    private Allele refAllele = null;
    

    // because the allele frequency priors are constant for a given i, we cache the results to avoid having to recompute everything
    //private double[] log10AlleleFrequencyPriors;

    // the allele frequency likelihoods (allocated once as an optimization)
    private ThreadLocal<double[]> log10AlleleFrequencyPosteriors = new ThreadLocal<double[]>();

    // the priors object
    private GenotypePriors genotypePriors;

    // samples in input
    private Set<String> samples = new TreeSet<String>();

    // the various loggers and writers
    private Logger logger = null;

    // fasta reference reader to supplement the edges of the reference sequence for long reads
    private IndexedFastaSequenceFile referenceReader;



    // the standard filter to use for calls below the confidence threshold but above the emit threshold
    private static final Set<String> filter = new HashSet<String>(1);

    private static final int MISMATCH_WINDOW_SIZE = 20;
	

	public byte[] CONTEXTREF = null;
	public Integer[] CYTOSINE_STATUS = null;
	
	protected double MAX_PHRED = 1000000;
	//public static double phredLiklihoodConfidance;
	
	

	public BisulfiteGenotyperEngine(GenomeAnalysisEngine toolkit,
			BisulfiteArgumentCollection BAC, Logger logger,
			VariantAnnotatorEngine engine,
			Set<String> samples) {
		//super(toolkit, BAC, logger, verboseWriter, engine, samples);
		this.samples = new TreeSet<String>(samples);
		initialize(toolkit, BAC, logger, engine, samples.size());
		// TODO Auto-generated constructor stub
	}
	

	protected void initialize(GenomeAnalysisEngine toolkit, BisulfiteArgumentCollection BAC, Logger logger, VariantAnnotatorEngine engine, int numSamples) {
        // note that, because we cap the base quality by the mapping quality, minMQ cannot be less than minBQ
        this.BAC = BAC.clone();
        
        
        this.logger = logger;
        this.annotationEngine = engine;
        
        genotypePriors = BisulfiteGenotyperEngine.createGenotypePriors(BAC);

        try {
            referenceReader = new CachingIndexedFastaSequenceFile(toolkit.getArguments().referenceFile);
        }
        catch(FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile(toolkit.getArguments().referenceFile,ex);
        }
    }
	
	public void setCytosineTypeStatus(CytosineTypeStatus cts, byte[] contextRef){
		this.cts = cts;
		this.CONTEXTREF = contextRef;
	}
	
	//@Override
	/**
     * Compute full calls at a given locus.
     *
     * @param tracker    the meta data tracker
     * @param refContext the reference base
     * @param rawContext contextual information around the locus
     * @return the VariantCallContext object
     */
    public BisulfiteVariantCallContext calculateLikelihoodsAndGenotypes(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext) {
        Map<String, StratifiedAlignmentContext> stratifiedContexts = getFilteredAndStratifiedContexts(BAC, refContext, rawContext);
        Map<String, BiallelicGenotypeLikelihoods> GLs = new HashMap<String, BiallelicGenotypeLikelihoods>();
        VariantContext vc = calculateLikelihoods(tracker, refContext, stratifiedContexts, StratifiedAlignmentContext.StratifiedContextType.COMPLETE, null, GLs);
        if ( vc == null )
            return null;
            

        BisulfiteVariantCallContext vcc = calculateGenotypes(tracker, refContext, rawContext, stratifiedContexts, GLs, vc);
        //System.out.println(vcc.vc.getGenotypesSortedByName());
       
       // if ( vcc == null ){
        	
       // 	System.out.println("vcc null--position: " + refContext.getLocus().getStart() + " refBase: " + refContext.getBase());
        	
       // }
       // if(!vcc.confidentlyCalled)
        //	System.out.println("vcc no confident--position: " + refContext.getLocus().getStart() + " refBase: " + refContext.getBase());
        
        return vcc;
    }
	
	//can not override unifiedGenotyperEngine, as it is static method
	protected static GenotypePriors createGenotypePriors(BisulfiteArgumentCollection BAC) {
        return new BisulfiteDiploidSNPGenotypePriors();
    }
	


    protected VariantContext calculateLikelihoods(RefMetaDataTracker tracker, ReferenceContext refContext, Map<String, StratifiedAlignmentContext> stratifiedContexts, StratifiedAlignmentContext.StratifiedContextType type, Allele alternateAlleleToUse, Map<String, BiallelicGenotypeLikelihoods> GLs) {
		if ( stratifiedContexts == null ){
			//System.out.println("no stratifiedContexts now");
			 return null;
		}

        // initialize the data for this thread if that hasn't been done yet
        if ( glcm.get() == null ) {
            glcm.set(BisulfiteGenotyperEngine.getGenotypeLikelihoodsCalculationObject(logger, BAC));
            
        }

        
        BisulfiteSNPGenotypeLikelihoodsCalculationModel bglcm = (BisulfiteSNPGenotypeLikelihoodsCalculationModel) glcm.get();
        bglcm.initialize(cts, BAC, CONTEXTREF);
        refAllele = bglcm.getLikelihoods(tracker, refContext, stratifiedContexts, type, genotypePriors, GLs, alternateAlleleToUse);
        
        CYTOSINE_STATUS = bglcm.getCytosineStatus();
        //CYTOSINE_TYPE_STATUS = bglcm.getCytosineTypeStatus();
        if (refAllele != null)
            return createVariantContextFromLikelihoods(refContext, refAllele, GLs);
        else
            return null;
          
        
    }
    
    protected void assignAFPosteriors(double[]likelihoods, double[] log10AFPosteriors){
    	for(int i = 0; i < likelihoods.length; i++){
    		//System.err.println(likelihoods[i]);
    		
    		log10AFPosteriors[i] = likelihoods[i];
    	}
    		
    }
	

	protected BisulfiteVariantCallContext calculateGenotypes(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext, Map<String, StratifiedAlignmentContext> stratifiedContexts, Map<String, BiallelicGenotypeLikelihoods> GLs, VariantContext vc) {
		// initialize the data for this thread if that hasn't been done yet
        if ( glcm.get() == null ) {
            return null;
        }
        log10AlleleFrequencyPosteriors.set(new double[3]);
        for ( BiallelicGenotypeLikelihoods GL : GLs.values() ) {
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
                // technically, at this point our confidence in a reference call isn't accurately estimated
                //  because it didn't take into account samples with no data, so let's get a better estimate
            	double sum = 0.0;
            	for (int j = 1; j < log10AlleleFrequencyPosteriors.get().length; j++)       	
                	sum += normalizedPosteriors[j];
                    
               double PofF = Math.min(sum, 1.0);
                return estimateReferenceConfidence(stratifiedContexts, genotypePriors.getHeterozygosity(), true, 1.0 - PofF);
            }
            
            HashMap<String, Genotype> genotypes = new HashMap<String, Genotype>();
            
            HashMap<String, Object> attributes = new HashMap<String, Object>();

            String rsID = DbSNPHelper.rsIDOfFirstRealSNP(tracker.getReferenceMetaData(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME));
            if ( rsID != null )
                attributes.put(VariantContext.ID_KEY, rsID);

            // if the site was downsampled, record that fact
            if ( rawContext.hasPileupBeenDownsampled() )
                attributes.put(VCFConstants.DOWNSAMPLED_KEY, true);

            Set<Allele> myAlleles = vc.getAlleles();
 
            GenomeLoc loc = refContext.getLocus();

            int endLoc = calculateEndPos(vc.getAlleles(), vc.getReference(), loc);
            
            assignGenotypes(vc,log10AlleleFrequencyPosteriors.get(),bestAFguess, genotypes, logRatio/10.0, GL);
            
            // strip out the alternate allele if it's a ref call
           // if ( bestAFguess == 0 && BAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.DISCOVERY ) {
           //     myAlleles = new HashSet<Allele>(1);
           //     myAlleles.add(refAllele);
           // }
            
            
            
            double cytosineMethyLevel = 0;
            //System.err.println("myAlleles: " + genotypes.values().toString());
            //genotypes.put(GL.getSample(), new Genotype(GL.getSample(), myAlleles, logRatio/10.0, null, attributes, false));
            
            if(BAC.ASSUME_SINGLE_SAMPLE != null){
            	 if(passesCallThreshold(logRatio)){
            		 
            		 for(Genotype genotypeTemp : genotypes.values()){
            			 //System.err.println("genotype: " + genotypeTemp.getGenotypeString() + "\tloc" + loc.getStart());
            			 //System.err.println("genotype: " + (char)genotypeTemp.getAllele(0).getBases()[0]);
            			 //System.err.println("genotype: " + genotypeTemp.getType());
            			 if(genotypeTemp.isHomRef()){
                 			 if(genotypeTemp.getAllele(0).getBases()[0]==BaseUtils.C){
                 				attributes.put(BisulfiteVCFConstants.NUMBER_OF_C_KEY, CYTOSINE_STATUS[1]);
                                attributes.put(BisulfiteVCFConstants.NUMBER_OF_T_KEY, CYTOSINE_STATUS[3]);
                                attributes.put(BisulfiteVCFConstants.C_IN_NEG_STRAND_KEY, false);
                                cytosineMethyLevel = (double)CYTOSINE_STATUS[1]/(double)(CYTOSINE_STATUS[1] + CYTOSINE_STATUS[3]);
                                attributes.put(BisulfiteVCFConstants.CYTOSINE_TYPE, getCytosineTypeStatus(false, cytosineMethyLevel));
                               
                                
                 			 }
                 			 else if(genotypeTemp.getAllele(0).getBases()[0]==BaseUtils.G){
                 				attributes.put(BisulfiteVCFConstants.NUMBER_OF_C_KEY, CYTOSINE_STATUS[0]);
                                attributes.put(BisulfiteVCFConstants.NUMBER_OF_T_KEY, CYTOSINE_STATUS[2]);
                                attributes.put(BisulfiteVCFConstants.C_IN_NEG_STRAND_KEY, true);
                                cytosineMethyLevel = (double)CYTOSINE_STATUS[0]/(double)(CYTOSINE_STATUS[0] + CYTOSINE_STATUS[2]);
                                attributes.put(BisulfiteVCFConstants.CYTOSINE_TYPE, getCytosineTypeStatus(true, cytosineMethyLevel));
                                
                 			 }
                 		 }
                 		 else if(genotypeTemp.isHomVar()){
                 			if(genotypeTemp.getAllele(1).getBases()[0]==BaseUtils.C){
                 				attributes.put(BisulfiteVCFConstants.NUMBER_OF_C_KEY, CYTOSINE_STATUS[1]);
                                attributes.put(BisulfiteVCFConstants.NUMBER_OF_T_KEY, CYTOSINE_STATUS[3]);
                                attributes.put(BisulfiteVCFConstants.C_IN_NEG_STRAND_KEY, false);
                                cytosineMethyLevel = (double)CYTOSINE_STATUS[1]/(double)(CYTOSINE_STATUS[1] + CYTOSINE_STATUS[3]);
                                attributes.put(BisulfiteVCFConstants.CYTOSINE_TYPE, getCytosineTypeStatus(false, cytosineMethyLevel));
                                
                 			 }
                 			 else if(genotypeTemp.getAllele(1).getBases()[0]==BaseUtils.G){
                 				attributes.put(BisulfiteVCFConstants.NUMBER_OF_C_KEY, CYTOSINE_STATUS[0]);
                                attributes.put(BisulfiteVCFConstants.NUMBER_OF_T_KEY, CYTOSINE_STATUS[2]);
                                attributes.put(BisulfiteVCFConstants.C_IN_NEG_STRAND_KEY, true);
                                cytosineMethyLevel = (double)CYTOSINE_STATUS[0]/(double)(CYTOSINE_STATUS[0] + CYTOSINE_STATUS[2]);
                                attributes.put(BisulfiteVCFConstants.CYTOSINE_TYPE, getCytosineTypeStatus(true, cytosineMethyLevel));
                                
                 			 }
                 		 }
            			 attributes.put(genotypeTemp.getType().toString(), true);
                 	 }
            		 if(CYTOSINE_STATUS[0] + CYTOSINE_STATUS[1] > 0 && loc.getStart() == BAC.testLocus){
            			// System.err.println("CYTOSINE_STATUS[0]: " + CYTOSINE_STATUS[0] + "\tCYTOSINE_STATUS[1]: " + CYTOSINE_STATUS[1] + "\tCYTOSINE_STATUS[2]: " + CYTOSINE_STATUS[2] + "\tCYTOSINE_STATUS[3]: " + CYTOSINE_STATUS[3]);
                		// System.err.println("cytosineMethyLevel: " + cytosineMethyLevel + "\tcts: " + cts.chgMethyLevel + "\t" + cts.chhMethyLevel + "\t" + cts.cpgMethyLevel + "\t" + logRatio);
                		 
            		 }
            		 
                 }
            }
           
           // genotypes.put(GL.getSample(), new Genotype(GL.getSample(), myAlleles, logRatio/10.0, null, attributes, false));
            
            VariantContext vcCall = new VariantContext("BG_call", loc.getContig(), loc.getStart(), endLoc,
                    myAlleles, genotypes, logRatio/10.0, passesCallThreshold(logRatio) ? null : filter, attributes);

            if ( annotationEngine != null ) {
                // first off, we want to use the *unfiltered* context for the annotations
                ReadBackedPileup pileup = null;
                if (rawContext.hasExtendedEventPileup())
                    pileup = rawContext.getExtendedEventPileup();
                else if (rawContext.hasBasePileup())
                    pileup = rawContext.getBasePileup();
                stratifiedContexts = StratifiedAlignmentContext.splitContextBySampleName(pileup, BAC.ASSUME_SINGLE_SAMPLE);

                Collection<VariantContext> variantContexts = annotationEngine.annotateContext(tracker, refContext, stratifiedContexts, vcCall);
                vcCall = variantContexts.iterator().next(); // we know the collection will always have exactly 1 element.
            }
            //if(vcCall != null){
            //	 System.out.println(vcCall.getChr() + "\t" + vcCall.getStart() + "\t" + vcCall.getEnd());
            //}
            
            
           
           
            BisulfiteVariantCallContext call = new BisulfiteVariantCallContext(vcCall, passesCallThreshold(logRatio), cts);
            call.setRefBase(refContext.getBase());
            return call;
        }
		return null;
	}
	
	protected boolean passesEmitThreshold(double conf, int bestAFguess) {
        return (BAC.OutputMode == OUTPUT_MODE.EMIT_ALL_CONFIDENT_SITES || bestAFguess != 0) && conf >= BAC.STANDARD_CONFIDENCE_FOR_CALLING;
    }
	
	
	protected static GenotypeLikelihoodsCalculationModel getGenotypeLikelihoodsCalculationObject(Logger logger, BisulfiteArgumentCollection BAC) {		
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

	        // now, test for bad pileups

	        // in all_bases mode, it doesn't matter
	        if ( BAC.OutputMode == OUTPUT_MODE.EMIT_ALL_SITES )
	            return true;

	        // is there no coverage?
	        if ( pileupSize == 0 )
	            return false;

	        // are there too many deletions in the pileup?
	        if ( (BAC.MAX_DELETION_FRACTION >=0 && BAC.MAX_DELETION_FRACTION <=1.0 ) &&
	                (double)numDeletions / (double)(pileupSize + numDeletions) > BAC.MAX_DELETION_FRACTION )
	            return false;

	        return true;
	    }

	 
	protected class BadBaseFilterBisulfite implements GATKSAMRecordFilter {
        private ReferenceContext refContext;
        private final BisulfiteArgumentCollection BAC;

        public BadBaseFilterBisulfite(ReferenceContext refContext, BisulfiteArgumentCollection BAC) {
            this.refContext = refContext;
            this.BAC = BAC;
        }

        @Override
        public BitSet getGoodBases(final GATKSAMRecord record) {
            // all bits are set to false by default
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
            //System.out.println(BAC.bisulfiteSpace);
            if(BAC.sequencingMode == MethylSNPModel.BM || BAC.sequencingMode == MethylSNPModel.GM){
            	mismatches = BisulfiteAlignmentUtils.mismatchesInRefWindow(record, refContext, BAC.MAX_MISMATCHES, MISMATCH_WINDOW_SIZE, BAC.sequencingMode);
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

        return new BisulfiteVariantCallContext(QualityUtils.phredScaleErrorRate(1.0 - P_of_ref) >= BAC.STANDARD_CONFIDENCE_FOR_CALLING, cts);
    }
	
	protected boolean passesCallThreshold(double conf) {
        return conf >= BAC.STANDARD_CONFIDENCE_FOR_CALLING;
    }
	
	private int calculateEndPos(Set<Allele> alleles, Allele refAllele, GenomeLoc loc) {
        // TODO - temp fix until we can deal with extended events properly
        // for indels, stop location is one more than ref allele length
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
		for(String cytosineType : cts.cytosineListMap.keySet()){
			String[] tmpKey = cytosineType.split("-");
			Double[] value = cts.cytosineListMap.get(cytosineType);
			//System.err.println("ctype: " + tmpKey[0]);
			//cTypeStatus = cTypeStatus + "," + tmpKey[0];
			if(value[cPos] >= Math.log10(BAC.cTypeThreshold)){
				cTypeStatus = cTypeStatus + "," + tmpKey[0];
				 if(Double.isNaN(cytosineMethyLevel))
						 cytosineMethyLevel = 0.0;
				//System.err.println("cTypeStatus: " + cTypeStatus + "\tLikelihood: " + value[cPos]);
				if(tmpKey[0].equalsIgnoreCase("CG")){
					
					cts.isCpg = true;
					cts.cpgMethyLevel = cytosineMethyLevel;
				}
				else if(tmpKey[0].equalsIgnoreCase("CHH")){
					cts.isChh = true;
					cts.chhMethyLevel = cytosineMethyLevel;
				}
				else if(tmpKey[0].equalsIgnoreCase("CHG")){
					cts.isChg = true;
					cts.chgMethyLevel = cytosineMethyLevel;
				}
				if(tmpKey[0].equalsIgnoreCase("GCH")){
					if(BAC.sequencingMode == MethylSNPModel.GM){
						cts.isGch = true;
						cts.gchMethyLevel = cytosineMethyLevel;
					}
					else{
						continue;
					}
					
				}
				else if(tmpKey[0].equalsIgnoreCase("GCG")){
					if(BAC.sequencingMode == MethylSNPModel.GM){
						cts.isGcg = true;
						cts.gcgMethyLevel = cytosineMethyLevel;
					}
					else{
						continue;
					}
					
				}
				else if(tmpKey[0].equalsIgnoreCase("HCG")){
					if(BAC.sequencingMode == MethylSNPModel.GM){
						cts.isHcg = true;
						cts.hcgMethyLevel = cytosineMethyLevel;
					}
					else{
						continue;
					}
					
				}
				
				
			}
			//if(cytosineMethyLevel > 0)
				//System.err.println("methy: " + value[2] + "\tliklihood_pos: " + value[0]  + "\tliklihood_neg: " + value[1]);
 		}
		
		return cTypeStatus;
	}

//    private static boolean isValidDeletionFraction(double d) {
 //       return ( d >= 0.0 && d <= 1.0 );
 //   }

	//protected static AlleleFrequencyCalculationModel getBisulfiteAlleleFrequencyCalculationObject(int N, Logger logger, PrintStream verboseWriter, BisulfiteArgumentCollection BAC) {
     //   AlleleFrequencyCalculationModel afcm;
      //  afcm = new BisulfiteExactAFCalculationModel(N, logger, verboseWriter);
       // return afcm;
    //}
	
	//COPY FROM gatk, NEED TO BE REPLACED
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

       // Allele refAllele = vc.getReference();
      //  Allele altAllele = vc.getAlternateAllele(0);
       
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


                
               // System.err.println("bestAF: " + bestAF + "\t" + log10AlleleFrequencyPosteriors[0] + "\t" + log10AlleleFrequencyPosteriors[1] + "\t" + log10AlleleFrequencyPosteriors[2] + "\t" + logRatio + "\t" + alleleA.getBaseString() + "\t" + alleleB.getBaseString());
                calls.put(sample.getKey(), new Genotype(sample.getKey(), myAlleles, logRatio, null, g.getAttributes(), false));
              }
        // return calls;
    }
	
    private VariantContext createVariantContextFromLikelihoods(ReferenceContext refContext, Allele refAllele, Map<String, BiallelicGenotypeLikelihoods> GLs) {
        // no-call everyone for now
        List<Allele> noCall = new ArrayList<Allele>();
        noCall.add(Allele.NO_CALL);

        Set<Allele> alleles = new HashSet<Allele>();
        alleles.add(refAllele);
        boolean addedAltAllele = false;

        HashMap<String, Genotype> genotypes = new HashMap<String, Genotype>();
        for ( BiallelicGenotypeLikelihoods GL : GLs.values() ) {
            if ( !addedAltAllele ) {
                addedAltAllele = true;
                alleles.add(GL.getAlleleA());
                alleles.add(GL.getAlleleB());
            }

            HashMap<String, Object> attributes = new HashMap<String, Object>();
            GenotypeLikelihoods likelihoods = new GenotypeLikelihoods(GL.getLikelihoods());
            attributes.put(VCFConstants.DEPTH_KEY, GL.getDepth());
            attributes.put(VCFConstants.GENOTYPE_LIKELIHOODS_KEY, likelihoods.getAsString());

            genotypes.put(GL.getSample(), new Genotype(GL.getSample(), noCall, Genotype.NO_NEG_LOG_10PERROR, null, attributes, false));
            //System.err.println("refAllele: " + refAllele.getBaseString() + "\tGL.getAlleleA(): " + GL.getAlleleA() + "\tGL.getAlleleB(): " + GL.getAlleleB());
        }

        GenomeLoc loc = refContext.getLocus();
        int endLoc = calculateEndPos(alleles, refAllele, loc);
       // for(Allele a : alleles){
       // 	 System.err.println(a.getBaseString());
        //}
       
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
    
}
