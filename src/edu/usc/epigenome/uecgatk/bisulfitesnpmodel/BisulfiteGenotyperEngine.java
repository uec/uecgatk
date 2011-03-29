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

import edu.usc.epigenome.uecgatk.bisulfitesnpmodel.BisulfiteAlignmentUtils;
import edu.usc.epigenome.uecgatk.bisulfitesnpmodel.BisulfiteSNPGenotypeLikelihoodsCalculationModel;

public class BisulfiteGenotyperEngine extends UnifiedGenotyperEngine {

	public boolean bisulfiteSpace = false;
	
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

        Allele refAllele = glcm.get().getLikelihoods(tracker, refContext, stratifiedContexts, type, genotypePriors, GLs, alternateAlleleToUse);
        
        
        if (refAllele != null){
        	//System.out.println("there is refAllele now");
        	return createVariantContextFromLikelihoods(refContext, refAllele, GLs);
        }   
        else{
        	//System.out.println("no refAllele now");
			 return null;
        }
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
