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

public class BisulfiteGenotyperEngine extends UnifiedGenotyperEngine {

	public boolean bisulfiteSpace = false;
	public static byte[] CONTEXTSEQ = null;
	protected double MAX_PHRED = 1000000;
	//public static double phredLiklihoodConfidance;
	
	public BisulfiteGenotyperEngine(GenomeAnalysisEngine toolkit,
			UnifiedArgumentCollection UAC) {
		super(toolkit, UAC);
		

		// TODO Auto-generated constructor stub
	}

	public BisulfiteGenotyperEngine(GenomeAnalysisEngine toolkit,
			UnifiedArgumentCollection UAC, Logger logger,
			PrintStream verboseWriter, VariantAnnotatorEngine engine,
			Set<String> samples) {
		super(toolkit, UAC, logger, verboseWriter, engine, samples);

		// TODO Auto-generated constructor stub
	}
	
	@Override
	protected void initialize(GenomeAnalysisEngine toolkit, UnifiedArgumentCollection UAC, Logger logger, PrintStream verboseWriter, VariantAnnotatorEngine engine, int numSamples) {
        // note that, because we cap the base quality by the mapping quality, minMQ cannot be less than minBQ
        this.UAC = UAC.clone();
        this.UAC.MIN_MAPPING_QUALTY_SCORE = Math.max(UAC.MIN_MAPPING_QUALTY_SCORE, UAC.MIN_BASE_QUALTY_SCORE);
        this.bisulfiteSpace = UAC.bisulfiteSpace;
        this.logger = logger;
        this.verboseWriter = verboseWriter;
        this.annotationEngine = engine;

        N = 2 * numSamples;
        log10AlleleFrequencyPriors = new double[N+1];
        computeAlleleFrequencyPriors(N);
        genotypePriors = BisulfiteGenotyperEngine.createGenotypePriors(UAC);

        filter.add(LOW_QUAL_FILTER_NAME);

        try {
            referenceReader = new CachingIndexedFastaSequenceFile(toolkit.getArguments().referenceFile);
        }
        catch(FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile(toolkit.getArguments().referenceFile,ex);
        }
    }
	
	public static void setContextSeq(byte[] contextSeq){
		CONTEXTSEQ = contextSeq;
	}
	
	@Override
	/**
     * Compute full calls at a given locus.
     *
     * @param tracker    the meta data tracker
     * @param refContext the reference base
     * @param rawContext contextual information around the locus
     * @return the VariantCallContext object
     */
    public VariantCallContext calculateLikelihoodsAndGenotypes(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext) {
        Map<String, StratifiedAlignmentContext> stratifiedContexts = getFilteredAndStratifiedContexts(UAC, refContext, rawContext);
        VariantContext vc = calculateLikelihoods(tracker, refContext, stratifiedContexts, StratifiedAlignmentContext.StratifiedContextType.COMPLETE, null);
        if ( vc == null ){
        	
        //	System.out.println("vc null--position: " + refContext.getLocus().getStart() + " refBase: " + refContext.getBase());
        	
        	return null;
        }
            

        VariantCallContext vcc = calculateGenotypes(tracker, refContext, rawContext, stratifiedContexts, vc);
        //System.out.println(vcc.vc.getGenotypesSortedByName());
        vcc.vc = GLsToPLs(vcc.vc);
       // if ( vcc == null ){
        	
       // 	System.out.println("vcc null--position: " + refContext.getLocus().getStart() + " refBase: " + refContext.getBase());
        	
       // }
       // if(!vcc.confidentlyCalled)
        //	System.out.println("vcc no confident--position: " + refContext.getLocus().getStart() + " refBase: " + refContext.getBase());
        
        return vcc;
    }
	
	//can not override unifiedGenotyperEngine, as it is static method
	protected static GenotypePriors createGenotypePriors(UnifiedArgumentCollection UAC) {
        GenotypePriors priors;
        switch ( UAC.GLmodel ) {
        	case SNP:
            // use flat priors for GLs
        		priors = new DiploidSNPGenotypePriors();
        		break;
        	case DINDEL:
            // create flat priors for Indels, actual priors will depend on event length to be genotyped
        		priors = new DiploidIndelGenotypePriors();
        		break;
        	case BSSNP:
            // create flat priors for Indels, actual priors will depend on event length to be genotyped
        		priors = new BisulfiteDiploidSNPGenotypePriors();
        		break;
            default: throw new IllegalArgumentException("Unexpected GenotypeCalculationModel " + UAC.GLmodel);
        }

        return priors;
    }
	

	@Override
    protected VariantContext calculateLikelihoods(RefMetaDataTracker tracker, ReferenceContext refContext, Map<String, StratifiedAlignmentContext> stratifiedContexts, StratifiedAlignmentContext.StratifiedContextType type, Allele alternateAlleleToUse) {
		if ( stratifiedContexts == null ){
			//System.out.println("no stratifiedContexts now");
			 return null;
		}

        // initialize the data for this thread if that hasn't been done yet
        if ( glcm.get() == null ) {
            glcm.set(BisulfiteGenotyperEngine.getGenotypeLikelihoodsCalculationObject(logger, UAC));
            
        }

        Map<String, BiallelicGenotypeLikelihoods> GLs = new HashMap<String, BiallelicGenotypeLikelihoods>();
        BisulfiteSNPGenotypeLikelihoodsCalculationModel bglcm = (BisulfiteSNPGenotypeLikelihoodsCalculationModel) glcm.get();
        bglcm.initialize(CONTEXTSEQ);
        Allele refAllele = bglcm.getLikelihoods(tracker, refContext, stratifiedContexts, type, genotypePriors, GLs, alternateAlleleToUse);
        
        
        
        if (refAllele != null){
        	//System.out.println("there is refAllele now");
        	return createVariantContextFromLikelihoods(refContext, refAllele, GLs);
        }   
        else{
        	//System.out.println("no refAllele now");
			 return null;
        }
    }
	
	@Override
	protected VariantCallContext calculateGenotypes(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc) {
		// initialize the data for this thread if that hasn't been done yet
        if ( afcm.get() == null ) {
            log10AlleleFrequencyPosteriors.set(new double[N+1]);
            afcm.set(getAlleleFrequencyCalculationObject(N, logger, verboseWriter, UAC));
        }

        // estimate our confidence in a reference call and return
        if ( vc.getNSamples() == 0 )
            return estimateReferenceConfidence(stratifiedContexts, genotypePriors.getHeterozygosity(), false, 1.0);

        // 'zero' out the AFs (so that we don't have to worry if not all samples have reads at this position)
        clearAFarray(log10AlleleFrequencyPosteriors.get());
        afcm.get().getLog10PNonRef(tracker, refContext, vc.getGenotypes(), log10AlleleFrequencyPriors, log10AlleleFrequencyPosteriors.get());

        // find the most likely frequency
        int bestAFguess = MathUtils.maxElementIndex(log10AlleleFrequencyPosteriors.get());
        int secondAFguess = MathUtils.minElementIndex(log10AlleleFrequencyPosteriors.get());
        for (int i = 0; i <= N; i++){
        	if(i != bestAFguess){	
        		if(log10AlleleFrequencyPosteriors.get()[i] >= log10AlleleFrequencyPosteriors.get()[secondAFguess]){
        			secondAFguess = i;
        		}
        	}
        }
        // calculate p(f>0)
        double[] normalizedPosteriors = MathUtils.normalizeFromLog10(log10AlleleFrequencyPosteriors.get());
        double logRatio;
        if(Double.isInfinite( Math.log10(normalizedPosteriors[secondAFguess]))){
        	logRatio = MAX_PHRED;
        }
        else{
        	logRatio = 10*(Math.log10(normalizedPosteriors[bestAFguess]) - Math.log10(normalizedPosteriors[secondAFguess]));
        }
       
        //for (int i = 0; i <= N; i++){
        //	System.out.println("AFposterior: " + log10AlleleFrequencyPosteriors.get()[i]);
        //	System.out.println("normalizeAFposterior: " + normalizedPosteriors[i]);
        //}
        //System.out.println("normalizeBestAFposterior: " + sum);
        //for (int i = 1; i <= N; i++)       	
        	//sum += normalizedPosteriors[i];
            
        //double PofF = Math.min(sum, 1.0); // deal with precision errors

       // double phredScaledConfidence;
       
           // phredScaledConfidence = QualityUtils.phredScaleErrorRate(1.0 - PofF);
          //  if ( Double.isInfinite(phredScaledConfidence) ) {
          //  	phredScaledConfidence = MAX_PHRED;
            	//   sum = 0.0;
                //for (int i = 1; i <= N; i++) {
                 //   if ( log10AlleleFrequencyPosteriors.get()[i] == AlleleFrequencyCalculationModel.VALUE_NOT_CALCULATED )
                  //      break;
                   // sum += log10AlleleFrequencyPosteriors.get()[i];
                //}
               // phredScaledConfidence = (MathUtils.compareDoubles(sum, 0.0) == 0 ? 0 : -10.0 * sum);
        //    }
            
       // System.out.println("phredScaledConfidence: " + phredScaledConfidence);
        // return a null call if we don't pass the confidence cutoff or the most likely allele frequency is zero
        if ( UAC.OutputMode != OUTPUT_MODE.EMIT_ALL_SITES && !passesEmitThreshold(logRatio, bestAFguess) ) {
            // technically, at this point our confidence in a reference call isn't accurately estimated
            //  because it didn't take into account samples with no data, so let's get a better estimate
        	double sum = 0.0;
        	for (int j = 1; j <= N; j++)       	
            	sum += normalizedPosteriors[j];
                
           double PofF = Math.min(sum, 1.0);
            return estimateReferenceConfidence(stratifiedContexts, genotypePriors.getHeterozygosity(), true, 1.0 - PofF);
        }

        // create the genotypes
        Map<String, Genotype> genotypes = afcm.get().assignGenotypes(vc, log10AlleleFrequencyPosteriors.get(), bestAFguess);

        // print out stats if we have a writer
        if ( verboseWriter != null )
            printVerboseData(refContext.getLocus().toString(), vc, logRatio, logRatio, normalizedPosteriors);

        // *** note that calculating strand bias involves overwriting data structures, so we do that last
        HashMap<String, Object> attributes = new HashMap<String, Object>();

        String rsID = DbSNPHelper.rsIDOfFirstRealSNP(tracker.getReferenceMetaData(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME));
        if ( rsID != null )
            attributes.put(VariantContext.ID_KEY, rsID);

        // if the site was downsampled, record that fact
        if ( rawContext.hasPileupBeenDownsampled() )
            attributes.put(VCFConstants.DOWNSAMPLED_KEY, true);

/*
        if ( !UAC.NO_SLOD && bestAFguess != 0 ) {
            final boolean DEBUG_SLOD = false;

            // the overall lod
            //double overallLog10PofNull = log10AlleleFrequencyPosteriors.get()[0];
            double overallLog10PofF = MathUtils.log10sumLog10(log10AlleleFrequencyPosteriors.get(), 1);
            if ( DEBUG_SLOD ) System.out.println("overallLog10PofF=" + overallLog10PofF);

            // the forward lod
            VariantContext vcForward = calculateLikelihoods(tracker, refContext, stratifiedContexts, StratifiedAlignmentContext.StratifiedContextType.FORWARD, vc.getAlternateAllele(0));
            clearAFarray(log10AlleleFrequencyPosteriors.get());
            afcm.get().getLog10PNonRef(tracker, refContext, vcForward.getGenotypes(), log10AlleleFrequencyPriors, log10AlleleFrequencyPosteriors.get());
            //double[] normalizedLog10Posteriors = MathUtils.normalizeFromLog10(log10AlleleFrequencyPosteriors.get(), true);
            double forwardLog10PofNull = log10AlleleFrequencyPosteriors.get()[0];
            double forwardLog10PofF = MathUtils.log10sumLog10(log10AlleleFrequencyPosteriors.get(), 1);
            if ( DEBUG_SLOD ) System.out.println("forwardLog10PofNull=" + forwardLog10PofNull + ", forwardLog10PofF=" + forwardLog10PofF);

            // the reverse lod
            VariantContext vcReverse = calculateLikelihoods(tracker, refContext, stratifiedContexts, StratifiedAlignmentContext.StratifiedContextType.REVERSE, vc.getAlternateAllele(0));
            clearAFarray(log10AlleleFrequencyPosteriors.get());
            afcm.get().getLog10PNonRef(tracker, refContext, vcReverse.getGenotypes(), log10AlleleFrequencyPriors, log10AlleleFrequencyPosteriors.get());
            //normalizedLog10Posteriors = MathUtils.normalizeFromLog10(log10AlleleFrequencyPosteriors.get(), true);
            double reverseLog10PofNull = log10AlleleFrequencyPosteriors.get()[0];
            double reverseLog10PofF = MathUtils.log10sumLog10(log10AlleleFrequencyPosteriors.get(), 1);
            if ( DEBUG_SLOD ) System.out.println("reverseLog10PofNull=" + reverseLog10PofNull + ", reverseLog10PofF=" + reverseLog10PofF);

            double forwardLod = forwardLog10PofF + reverseLog10PofNull - overallLog10PofF;
            double reverseLod = reverseLog10PofF + forwardLog10PofNull - overallLog10PofF;
            if ( DEBUG_SLOD ) System.out.println("forward lod=" + forwardLod + ", reverse lod=" + reverseLod);

            // strand score is max bias between forward and reverse strands
            double strandScore = Math.max(forwardLod, reverseLod);
            // rescale by a factor of 10
            strandScore *= 10.0;
            //logger.debug(String.format("SLOD=%f", strandScore));

            attributes.put("SB", Double.valueOf(strandScore));
        }
*/
        GenomeLoc loc = refContext.getLocus();

        int endLoc = calculateEndPos(vc.getAlleles(), vc.getReference(), loc);

        Set<Allele> myAlleles = vc.getAlleles();
        // strip out the alternate allele if it's a ref call
        if ( bestAFguess == 0 && UAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.DISCOVERY ) {
            myAlleles = new HashSet<Allele>(1);
            myAlleles.add(vc.getReference());
        }
        VariantContext vcCall = new VariantContext("UG_call", loc.getContig(), loc.getStart(), endLoc,
                myAlleles, genotypes, logRatio/10.0, passesCallThreshold(logRatio) ? null : filter, attributes);

        if ( annotationEngine != null ) {
            // first off, we want to use the *unfiltered* context for the annotations
            ReadBackedPileup pileup = null;
            if (rawContext.hasExtendedEventPileup())
                pileup = rawContext.getExtendedEventPileup();
            else if (rawContext.hasBasePileup())
                pileup = rawContext.getBasePileup();
            stratifiedContexts = StratifiedAlignmentContext.splitContextBySampleName(pileup, UAC.ASSUME_SINGLE_SAMPLE);

            Collection<VariantContext> variantContexts = annotationEngine.annotateContext(tracker, refContext, stratifiedContexts, vcCall);
            vcCall = variantContexts.iterator().next(); // we know the collection will always have exactly 1 element.
        }
        //if(vcCall != null){
        //	 System.out.println(vcCall.getChr() + "\t" + vcCall.getStart() + "\t" + vcCall.getEnd());
        //}
       
        VariantCallContext call = new VariantCallContext(vcCall, passesCallThreshold(logRatio));
        call.setRefBase(refContext.getBase());
        return call;
	}
	
	@Override
	protected boolean passesEmitThreshold(double conf, int bestAFguess) {
        return (UAC.OutputMode == OUTPUT_MODE.EMIT_ALL_CONFIDENT_SITES || bestAFguess != 0) && conf >= UAC.STANDARD_CONFIDENCE_FOR_CALLING;
    }
	
	
	protected static GenotypeLikelihoodsCalculationModel getGenotypeLikelihoodsCalculationObject(Logger logger, UnifiedArgumentCollection UAC) {
		GenotypeLikelihoodsCalculationModel glcm;
        switch ( UAC.GLmodel ) {
            case SNP:
                glcm = new SNPGenotypeLikelihoodsCalculationModel(UAC, logger);
                break;
           case DINDEL:
                glcm = new DindelGenotypeLikelihoodsCalculationModel(UAC, logger);
                break;
           case BSSNP:
               glcm = new BisulfiteSNPGenotypeLikelihoodsCalculationModel(UAC, logger);
               break;
            default: throw new IllegalArgumentException("Unexpected GenotypeCalculationModel " + UAC.GLmodel);
        }

        return glcm;
    }

	
	
	@Override
	protected Map<String, StratifiedAlignmentContext> getFilteredAndStratifiedContexts(UnifiedArgumentCollection UAC, ReferenceContext refContext, AlignmentContext rawContext) {
		BadBaseFilterBisulfite badReadPileupFilter = new BadBaseFilterBisulfite(refContext, UAC);

        Map<String, StratifiedAlignmentContext> stratifiedContexts = null;
       
        if ( UAC.GLmodel == GenotypeLikelihoodsCalculationModel.Model.DINDEL && rawContext.hasExtendedEventPileup() ) {

            ReadBackedExtendedEventPileup rawPileup = rawContext.getExtendedEventPileup();

            // filter the context based on min mapping quality
            ReadBackedExtendedEventPileup pileup = rawPileup.getMappingFilteredPileup(UAC.MIN_MAPPING_QUALTY_SCORE);

            // don't call when there is no coverage
            if ( pileup.size() == 0 && UAC.OutputMode != OUTPUT_MODE.EMIT_ALL_SITES )
                return null;

            // stratify the AlignmentContext and cut by sample
            stratifiedContexts = StratifiedAlignmentContext.splitContextBySampleName(pileup, UAC.ASSUME_SINGLE_SAMPLE);

        } else if ( (UAC.GLmodel == GenotypeLikelihoodsCalculationModel.Model.SNP || UAC.GLmodel == GenotypeLikelihoodsCalculationModel.Model.BSSNP) && !rawContext.hasExtendedEventPileup() ) {

            byte ref = refContext.getBase();
            if ( !BaseUtils.isRegularBase(ref) )
                return null;

            // stratify the AlignmentContext and cut by sample
            stratifiedContexts = StratifiedAlignmentContext.splitContextBySampleName(rawContext.getBasePileup(), UAC.ASSUME_SINGLE_SAMPLE);
            //System.out.println("here!");
            //System.out.println(rawContext.getBasePileup().size());
            // filter the reads (and test for bad pileups)
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
	                    if ( read.getMappingQuality() >= UAC.MIN_MAPPING_QUALTY_SCORE &&
	                         (UAC.USE_BADLY_MATED_READS || !BadMateFilter.hasBadMate(read)) )
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
	        if ( UAC.OutputMode == OUTPUT_MODE.EMIT_ALL_SITES )
	            return true;

	        // is there no coverage?
	        if ( pileupSize == 0 )
	            return false;

	        // are there too many deletions in the pileup?
	        if ( isValidDeletionFraction(UAC.MAX_DELETION_FRACTION) &&
	                (double)numDeletions / (double)(pileupSize + numDeletions) > UAC.MAX_DELETION_FRACTION )
	            return false;

	        return true;
	    }

	 
	protected class BadBaseFilterBisulfite implements GATKSAMRecordFilter {
        private ReferenceContext refContext;
        private final UnifiedArgumentCollection UAC;

        public BadBaseFilterBisulfite(ReferenceContext refContext, UnifiedArgumentCollection UAC) {
            this.refContext = refContext;
            this.UAC = UAC;
        }

        @Override
        public BitSet getGoodBases(final GATKSAMRecord record) {
            // all bits are set to false by default
            BitSet bitset = new BitSet(record.getReadLength());

            // if the mapping quality is too low or the mate is bad, we can just zero out the whole read and continue
            if ( record.getMappingQuality() < UAC.MIN_MAPPING_QUALTY_SCORE ||
                 (!UAC.USE_BADLY_MATED_READS && BadMateFilter.hasBadMate(record)) ) {
            	//System.out.println("bad mates");
            	return bitset;
            }

            byte[] quals = record.getBaseQualities();
            for (int i = 0; i < quals.length; i++) {
                if ( quals[i] >= UAC.MIN_BASE_QUALTY_SCORE )
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
            //System.out.println(UAC.bisulfiteSpace);
            if(UAC.bisulfiteSpace){
            	mismatches = BisulfiteAlignmentUtils.mismatchesInRefWindow(record, refContext, UAC.MAX_MISMATCHES, MISMATCH_WINDOW_SIZE, UAC.bisulfiteSpace);
            }
            else{
            	mismatches = AlignmentUtils.mismatchesInRefWindow(record, refContext, UAC.MAX_MISMATCHES, MISMATCH_WINDOW_SIZE);
            }
            if ( mismatches != null )
                bitset.and(mismatches);

            return bitset;
        }
    }
	
	@Override
	protected void computeAlleleFrequencyPriors(int N) {
        // calculate the allele frequency priors for 1-N
        double sum = 0.0;
        double heterozygosity;

        if (UAC.GLmodel == GenotypeLikelihoodsCalculationModel.Model.DINDEL)
            heterozygosity = UAC.INDEL_HETEROZYGOSITY;
        else
            heterozygosity = UAC.heterozygosity;
        if(N==2 && bisulfiteSpace){
        	for (int i = 0; i <= N; i++) {
                log10AlleleFrequencyPriors[i] = 0;
            }
        }
        else{
        	 for (int i = 1; i <= N; i++) {
                 double value = heterozygosity / (double)i;
                 log10AlleleFrequencyPriors[i] = Math.log10(value);
                 sum += value;
             }

             // null frequency for AF=0 is (1 - sum(all other frequencies))
             log10AlleleFrequencyPriors[0] = Math.log10(1.0 - sum);
        }
    }
	

	//protected static AlleleFrequencyCalculationModel getBisulfiteAlleleFrequencyCalculationObject(int N, Logger logger, PrintStream verboseWriter, UnifiedArgumentCollection UAC) {
     //   AlleleFrequencyCalculationModel afcm;
      //  afcm = new BisulfiteExactAFCalculationModel(N, logger, verboseWriter);
       // return afcm;
    //}
	
}
