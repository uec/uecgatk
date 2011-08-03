package edu.usc.epigenome.uecgatk.bisulfitesnpmodel;

import static java.lang.Math.log10;
import static java.lang.Math.pow;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMUtils;

import org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel.BisulfiteDiploidSNPGenotypePriors;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.genotyper.DiploidSNPGenotypeLikelihoods;
import org.broadinstitute.sting.gatk.walkers.genotyper.DiploidSNPGenotypePriors;
import org.broadinstitute.sting.gatk.walkers.genotyper.PerFragmentPileupElement;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.genotype.DiploidGenotype;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;


public class BisulfiteDiploidSNPGenotypeLikelihoods implements Cloneable  {
	protected BisulfiteDiploidSNPGenotypePriors priors = null;
    // TODO: don't calculate this each time through

	protected BisulfiteArgumentCollection BAC;
    protected double BISULFITE_CONVERSION_RATE;
    protected double OVER_CONVERSION_RATE;
 //   protected static double CPG_METHYLATION_RATE = 0;
 //   protected static double CPH_METHYLATION_RATE = 0;
    protected String DETERMINED_CYTOSINE_TYPE_POS = "";
    protected double DETERMINED_CYTOSINE_TYPE_C_LIKELIHOOD_POS = Double.NEGATIVE_INFINITY;
    protected double DETERMINED_CYTOSINE_TYPE_C_METHY_POS = 0;
    //protected static boolean DETERMINED_CYTOSINE_TYPE_NEGSTRAND = false;
    protected String DETERMINED_CYTOSINE_TYPE_NEG = "";
    protected double DETERMINED_CYTOSINE_TYPE_C_LIKELIHOOD_NEG = Double.NEGATIVE_INFINITY;
    protected double DETERMINED_CYTOSINE_TYPE_C_METHY_NEG = 0;
    protected double ABSOLUTE_CYTOSINE_TYPE_C_LIKELIHOOD = 0.05;
    
 //   protected static double UNDETERMINED_CYTOSINE_TYPE_C_METHY = 0;
    
    public static final double HUMAN_HETEROZYGOSITY = 1e-3;
 //   public static final double CEU_HETEROZYGOSITY = 1e-3;
 //   public static final double YRI_HETEROZYGOSITY = 1.0 / 850;
    public static final double PROB_OF_REFERENCE_ERROR = 1e-6;
    protected final static double log10_3 = log10(3.0);
    protected double[] log10Likelihoods = null;
    protected double[] log10Posteriors = null;


    // TODO: don't calculate this each time through
    protected double log10_PCR_error_3;
    protected double log10_1_minus_PCR_error;
    public boolean VERBOSE = false;
    /*
	public BisulfiteDiploidSNPGenotypeLikelihoods() {
		// TODO Auto-generated constructor stub
		super.priors = new DiploidSNPGenotypePriors();
        log10_PCR_error_3 = log10(DEFAULT_PCR_ERROR_RATE) - log10_3;
        this.priors = new BisulfiteDiploidSNPGenotypePriors();
        setToZeroBs();
	}

	public BisulfiteDiploidSNPGenotypeLikelihoods(
			DiploidSNPGenotypePriors priors, double PCR_error_rate) {
		super(priors, PCR_error_rate);
		// TODO Auto-generated constructor stub

        this.priors = new BisulfiteDiploidSNPGenotypePriors();
        setToZeroBs();
	}
	
	public BisulfiteDiploidSNPGenotypeLikelihoods(
			BisulfiteDiploidSNPGenotypePriors priors, double PCR_error_rate) {
		this.priors = priors;

        setToZeroBs();
	}
	*/
	public BisulfiteDiploidSNPGenotypeLikelihoods(
			RefMetaDataTracker tracker, ReferenceContext ref, BisulfiteDiploidSNPGenotypePriors priors, BisulfiteArgumentCollection BAC, double[] cytosineMethyStatus) {
		this.priors = priors;
		this.BAC = BAC;
		this.BISULFITE_CONVERSION_RATE = BAC.bsRate;
		this.OVER_CONVERSION_RATE = BAC.overRate;
		log10_PCR_error_3 = log10(BAC.PCR_error) - log10_3;
        log10_1_minus_PCR_error = log10(1.0 - BAC.PCR_error);
        this.DETERMINED_CYTOSINE_TYPE_C_METHY_POS = cytosineMethyStatus[0];
        this.DETERMINED_CYTOSINE_TYPE_C_METHY_NEG = cytosineMethyStatus[1];
		setToZeroBs();
//		this.CPG_METHYLATION_RATE = cpgMethyRate;
//		this.CPH_METHYLATION_RATE = cphMethyRate;
//		DETERMINED_CYTOSINE_TYPE_C_METHY = CPH_METHYLATION_RATE;
		
	}
	
	public void setPriors(RefMetaDataTracker tracker, ReferenceContext ref, double PCR_error_rate, double novelDbsnpHet, double validateDbsnpHet){
		
		this.priors.setPriors(tracker, ref, HUMAN_HETEROZYGOSITY, PROB_OF_REFERENCE_ERROR, novelDbsnpHet, validateDbsnpHet);
        setToZeroBs();
	}
	
	public void setPriorsBasedOnContextRef(RefMetaDataTracker tracker, ReferenceContext ref, double PCR_error_rate, double bisulfiteConversionRate, double novelDbsnpHet, double validateDbsnpHet, CytosineTypeStatus cts, byte[] contextRef){
		double refMethyStatus = 0;
		for(String cytosineType : cts.cytosineListMap.keySet()){
			String[] tmpKey = cytosineType.split("-");
			if(tmpKey[0].equalsIgnoreCase("C")){
				continue;
			}
			Integer cytosinePos = Integer.parseInt(tmpKey[1]) - 1;
			int i = 0;
			for(; i < tmpKey[0].length(); i++){
				//if(i == cytosinePos)
				//	continue;
				int elementOffset = i - cytosinePos + 100;
				if(VERBOSE){
					System.err.println("elementOffset: " + elementOffset + "\tcytosinePos: " + cytosinePos + "\ti: " + i + "\tcontextRef[elementOffset]" + contextRef[elementOffset] + "\ttmpKey[0].getBytes()[i]" + tmpKey[0].getBytes()[i]);
					//System.err.println( "value[0]: " + cts.cytosineListMap.get(cytosineType)[0] + "\tvalue[1]: " + cts.cytosineListMap.get(cytosineType)[1]);
				}
				if(elementOffset < 0 || elementOffset > contextRef.length-1)
					continue;
				if(!BaseUtilsMore.iupacCodeEqual(tmpKey[0].getBytes()[i], contextRef[elementOffset], false)){
					break;
				}
			}
			if(i >= tmpKey[0].length()){
				refMethyStatus = cts.cytosineListMap.get(cytosineType)[2];
				break;
			}
			i = 0;
			for(; i < tmpKey[0].length(); i++){
			//	if(i == cytosinePos)
			//		continue;
				int elementOffset = 100 - (i - cytosinePos);
				if(VERBOSE){
					System.err.println("elementOffset: " + elementOffset + "\tcytosinePos: " + cytosinePos + "\ti: " + i + "\tcontextRef[elementOffset]" + contextRef[elementOffset] + "\ttmpKey[0].getBytes()[i]" + tmpKey[0].getBytes()[i]);
					//System.err.println( "value[0]: " + cts.cytosineListMap.get(cytosineType)[0] + "\tvalue[1]: " + cts.cytosineListMap.get(cytosineType)[1]);
				}
				if(elementOffset < 0 || elementOffset > contextRef.length-1)
					continue;
				if(!BaseUtilsMore.iupacCodeEqual(BaseUtilsMore.iupacCodeComplement(tmpKey[0].getBytes()[i]), contextRef[elementOffset], false)){
					break;
				}
			}
			if(i >= tmpKey[0].length()){
				refMethyStatus = cts.cytosineListMap.get(cytosineType)[2];
				break;
			}	
			if(VERBOSE){
				System.err.println("cytosineType: " + cytosineType + "\trefMethyStatus: " + refMethyStatus + "\ti: " + i);
				//System.err.println( "value[0]: " + cts.cytosineListMap.get(cytosineType)[0] + "\tvalue[1]: " + cts.cytosineListMap.get(cytosineType)[1]);
			}
		}
		if(VERBOSE){
			System.err.println("refMethyStatus: " + refMethyStatus);
			//System.err.println( "value[0]: " + cts.cytosineListMap.get(cytosineType)[0] + "\tvalue[1]: " + cts.cytosineListMap.get(cytosineType)[1]);
		}
		this.priors.setPriors(tracker, ref, HUMAN_HETEROZYGOSITY, PROB_OF_REFERENCE_ERROR, BISULFITE_CONVERSION_RATE, OVER_CONVERSION_RATE, refMethyStatus, novelDbsnpHet, validateDbsnpHet, cts);
        setToZeroBs();
	}
	
	public void checkCytosineStatus(ReadBackedPileup pileup, CytosineTypeStatus cts, double threshold, boolean autoEstimateC, boolean secondIteration){
		//boolean useDefaulMethyLevel = true;
		//double defaultMethy = 0;
	//	HashMap<String, Double> cTypeLikelihood = new HashMap<String, Double>();
		//HashMap<String, Double> cTypeReverseLikelihood = new HashMap<String, Double>();
		for(String cytosineType : cts.cytosineListMap.keySet()){
			String[] tmpKey = cytosineType.split("-");
			//if(tmpKey[0].equalsIgnoreCase("C"))
			//	continue;
			if(VERBOSE){
				System.err.println(tmpKey[0]);
				//System.err.println( "value[0]: " + cts.cytosineListMap.get(cytosineType)[0] + "\tvalue[1]: " + cts.cytosineListMap.get(cytosineType)[1]);
			}
			
			Integer cytosinePos = Integer.parseInt(tmpKey[1]) - 1;
			double[] adjacentCytosineSeqLikelihood = new double[tmpKey[0].length()];
			double[] adjacentCytosineSeqLikelihoodReverseStrand = new double[tmpKey[0].length()];
			double[] adjacentNotCytosineSeqLikelihood = new double[tmpKey[0].length()];
			double[] adjacentNotCytosineSeqLikelihoodReverseStrand = new double[tmpKey[0].length()];
			ReadBackedPileup pileupPositiveStrand = pileup.getPositiveStrandPileup();
	        for(PileupElement p : pileupPositiveStrand ){
	        	SAMRecord samRecord = new SAMRecord(null);
	        	int cytosineOffset = p.getOffset();
	        	if(BAC.pairedEndMode){
	        		try {
						samRecord = (SAMRecord) p.getRead().clone();
					} catch (CloneNotSupportedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
		        	
		        	boolean Paired = samRecord.getReadPairedFlag();
		        	//boolean firstOfPair = samRecord.getFirstOfPairFlag();
		        	boolean secondOfPair = samRecord.getSecondOfPairFlag();
		        	//int cytosineOffset = p.getOffset();
		        	if (samRecord.getNotPrimaryAlignmentFlag())
					{
						continue;
					}
					
					// Inverted dups, count only one end
					if (samRecord.getAlignmentStart() == samRecord.getMateAlignmentStart() && samRecord.getReadNegativeStrandFlag() == samRecord.getMateNegativeStrandFlag())
					{
						if (samRecord.getSecondOfPairFlag()) continue;
						//System.err.printf("Inverted dup %d%s (%s)\n", samRecord.getAlignmentStart(), samRecord.getReadNegativeStrandFlag()?"-":"+", PicardUtils.getReadString(samRecord, true));
					}
		        	if (Paired  && !BAC.allowBadMates && !samRecord.getProperPairFlag())
					{
						continue;
					}
		        	
		        	if(secondOfPair){
		        		if(VERBOSE){
		        			System.err.println("NegativeStrandFlag: " + samRecord.getReadNegativeStrandFlag() + "\t" + "MateNegativeStrandFlag: " + samRecord.getMateNegativeStrandFlag() + "\tbase: " + samRecord.getReadBases()[cytosineOffset] + "\t" + "baseQ: " + samRecord.getBaseQualities()[cytosineOffset] + "\tbase-1: " + samRecord.getReadBases()[cytosineOffset-1] + "\t" + "baseQ-1: " + samRecord.getBaseQualities()[cytosineOffset-1] + "\tbase+1: " + samRecord.getReadBases()[cytosineOffset+1] + "\t" + "baseQ+1: " + samRecord.getBaseQualities()[cytosineOffset+1]);
		        			System.err.println("getReadString: " + samRecord.getReadString());
		        		}
		        		//samRecord.setReadBases(BaseUtils.simpleReverseComplement(samRecord.getReadBases()));
			        	samRecord.setReadNegativeStrandFlag(!samRecord.getReadNegativeStrandFlag());
			        	//cytosineOffset = samRecord.getReadLength() - 1 - cytosineOffset;
			        	//samRecord.setBaseQualities(BaseUtilsMore.simpleReverse(samRecord.getBaseQualities()));
			        	samRecord.setSecondOfPairFlag(!secondOfPair);
		        		if(VERBOSE){
		        			System.err.println("secondOfPair: " + secondOfPair + "\t" + "getMateAlignmentStart: " + samRecord.getMateAlignmentStart() + "\t" + "MateNegativeStrandFlag: " + samRecord.getMateNegativeStrandFlag());
			        		System.err.println("getReadString: " + samRecord.getReadString() + "\t" + "getAlignmentStart: " + samRecord.getAlignmentStart() + "\t" + "getUnclippedEnd: " + samRecord.getUnclippedEnd() + "\t" +"NegativeStrandFlag: " + samRecord.getReadNegativeStrandFlag() + "\tcytosineOffset: " + cytosineOffset + "\tbase: " + samRecord.getReadBases()[cytosineOffset] + "\t" + "baseQ: " + samRecord.getBaseQualities()[cytosineOffset] + "\tbase-1: " + samRecord.getReadBases()[cytosineOffset-1] + "\t" + "baseQ-1: " + samRecord.getBaseQualities()[cytosineOffset-1] + "\tbase+1: " + samRecord.getReadBases()[cytosineOffset+1] + "\t" + "baseQ+1: " + samRecord.getBaseQualities()[cytosineOffset+1]);
			        		
		        		}
		        			        		
		        	}
	        	}
	        	else{
	        		samRecord = p.getRead();
	        	}
				
	        	
				
				for(int i = 0; i < tmpKey[0].length(); i++){
					//if(i == cytosinePos)
					//	continue;
					int elementOffset = i - cytosinePos + cytosineOffset;
					if(elementOffset < 0 || elementOffset > samRecord.getReadLength()-1)
						continue;
					//System.err.println("1: " + elementOffset + "\t" + cytosinePos + "\t" + cytosineOffset + "\t" + samRecord.getAlignmentStart() );
					PileupElement tmpP = new PileupElement(samRecord,elementOffset);
					if ( !usableBase(tmpP, true) ){
			        	continue;
			        }
					byte qual = tmpP.getQual();
					if ( qual > SAMUtils.MAX_PHRED_SCORE )
			            throw new UserException.MalformedBAM(tmpP.getRead(), String.format("the maximum allowed quality score is %d, but a quality of %d was observed in read %s.  Perhaps your BAM incorrectly encodes the quality scores in Sanger format; see http://en.wikipedia.org/wiki/FASTQ_format for more details", SAMUtils.MAX_PHRED_SCORE, qual, tmpP.getRead().getReadName()));
			        qual = (byte)Math.min((int)tmpP.getQual(), tmpP.getMappingQual());
			        byte observedBase = tmpP.getBase();
			        byte qualityScore = qual;
			        byte queryByte = (byte)tmpKey[0].charAt(i);
			        double error = pow(10, -qualityScore/10.0);
			        if ( observedBase == 0 )
			            continue;
			        adjacentCytosineSeqLikelihood[i] += BaseUtilsMore.iupacCodeEqual(queryByte, observedBase, samRecord.getReadNegativeStrandFlag()) ? log10(1.0 - error) : log10(error/3.0);
			        adjacentNotCytosineSeqLikelihood[i] += BaseUtilsMore.iupacCodeEqual(queryByte, observedBase, samRecord.getReadNegativeStrandFlag()) ? log10(error/3.0) : log10(1.0 - error);
			        if(VERBOSE){
			        	//System.err.println("forwardStrand-1: queryByte: " + (char)queryByte + "\tobservedBase: " + (char)observedBase + "\t" + BaseUtilsMore.iupacCodeEqual(queryByte, observedBase, false) + "\telementOffset: " + elementOffset + "\tcytosinePos: " + cytosinePos + "\tcytosineOffset: " + cytosineOffset + "\ti: " + i);
			        }
				}
				
				//int cytosineOffsetOpposite = samRecord.getReadLength() - p.getOffset();
				for(int i = 0; i < tmpKey[0].length(); i++){
					//if(i == cytosinePos)
					//	continue;
					int elementOffset = cytosineOffset - (i - cytosinePos);
					if(elementOffset < 0 || elementOffset > samRecord.getReadLength()-1)
						continue;
					//System.err.println("2: " + elementOffset + "\t" + cytosinePos + "\t" + cytosineOffset + "\t" + samRecord.getAlignmentStart() );
					PileupElement tmpP = new PileupElement(samRecord,elementOffset);
					if ( !usableBase(tmpP, true) ){
			        	continue;
			        }
					byte qual = tmpP.getQual();
					if ( qual > SAMUtils.MAX_PHRED_SCORE )
			            throw new UserException.MalformedBAM(tmpP.getRead(), String.format("the maximum allowed quality score is %d, but a quality of %d was observed in read %s.  Perhaps your BAM incorrectly encodes the quality scores in Sanger format; see http://en.wikipedia.org/wiki/FASTQ_format for more details", SAMUtils.MAX_PHRED_SCORE, qual, tmpP.getRead().getReadName()));
			        qual = (byte)Math.min((int)tmpP.getQual(), tmpP.getMappingQual());
			        byte observedBase = tmpP.getBase();
			        byte qualityScore = qual;
			       // byte queryByte = (byte)tmpKey[0].charAt(tmpKey[0].length() - 1 - i);
			        byte queryByte = (byte)tmpKey[0].charAt(i);
			        queryByte = BaseUtilsMore.iupacCodeComplement(queryByte);
			        double error = pow(10, -qualityScore/10.0);
			        if ( observedBase == 0 )
			            continue;
			        adjacentCytosineSeqLikelihoodReverseStrand[i] += BaseUtilsMore.iupacCodeEqual(queryByte, observedBase, samRecord.getReadNegativeStrandFlag()) ? log10(1.0 - error) : log10(error/3.0);
			        adjacentNotCytosineSeqLikelihoodReverseStrand[i] += BaseUtilsMore.iupacCodeEqual(queryByte, observedBase, samRecord.getReadNegativeStrandFlag()) ? log10(error/3.0) : log10(1.0 - error);
			        if(VERBOSE){
			        	//System.err.println("reverseStrand-1: queryByte: " + (char)queryByte + "\tobservedBase: " + (char)observedBase + "\t" + BaseUtilsMore.iupacCodeEqual(queryByte, observedBase, false) + "\telementOffset: " + elementOffset + "\tcytosinePos: " + cytosinePos + "\tcytosineOffset: " + cytosineOffset + "\ti: " + i + "\tadjacentCytosineSeqLikelihoodReverseStrand[i]" + adjacentCytosineSeqLikelihoodReverseStrand[i]);
			        }
				}           
	        }
	        
	        ReadBackedPileup pileupNegativeStrand = pileup.getNegativeStrandPileup();
	        for(PileupElement p : pileupNegativeStrand ){
	        	SAMRecord samRecord = new SAMRecord(null);
	        	int cytosineOffset = p.getOffset();
	        	if(BAC.pairedEndMode){
	        		try {
						samRecord = (SAMRecord) p.getRead().clone();
					} catch (CloneNotSupportedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
		        	boolean secondOfPair = samRecord.getSecondOfPairFlag();
		        	//int cytosineOffset = p.getOffset();
		        	boolean Paired = samRecord.getReadPairedFlag();
		        	if (samRecord.getNotPrimaryAlignmentFlag())
					{
						continue;
					}
					
					// Inverted dups, count only one end
					if (samRecord.getAlignmentStart() == samRecord.getMateAlignmentStart() && samRecord.getReadNegativeStrandFlag() == samRecord.getMateNegativeStrandFlag())
					{
						if (samRecord.getSecondOfPairFlag()) continue;
						//System.err.printf("Inverted dup %d%s (%s)\n", samRecord.getAlignmentStart(), samRecord.getReadNegativeStrandFlag()?"-":"+", PicardUtils.getReadString(samRecord, true));
					}
		        	if (Paired  && !BAC.allowBadMates && !samRecord.getProperPairFlag())
					{
						continue;
					}
		        	if(secondOfPair){
		        		if(VERBOSE){
		        			System.err.println("NegativeStrandFlag: " + samRecord.getReadNegativeStrandFlag() + "\t" + "MateNegativeStrandFlag: " + samRecord.getMateNegativeStrandFlag() + "\tbase: " + samRecord.getReadBases()[cytosineOffset] + "\t" + "baseQ: " + samRecord.getBaseQualities()[cytosineOffset] + "\tbase-1: " + samRecord.getReadBases()[cytosineOffset-1] + "\t" + "baseQ-1: " + samRecord.getBaseQualities()[cytosineOffset-1] + "\tbase+1: " + samRecord.getReadBases()[cytosineOffset+1] + "\t" + "baseQ+1: " + samRecord.getBaseQualities()[cytosineOffset+1]);
		        			System.err.println("getReadString: " + samRecord.getReadString());
		        		}
		        		//samRecord.setReadBases(BaseUtils.simpleReverseComplement(samRecord.getReadBases()));
			        	samRecord.setReadNegativeStrandFlag(!samRecord.getReadNegativeStrandFlag());
			        	//cytosineOffset = samRecord.getReadLength() - 1 - cytosineOffset;
			        	//samRecord.setBaseQualities(BaseUtilsMore.simpleReverse(samRecord.getBaseQualities()));
			        	samRecord.setSecondOfPairFlag(!secondOfPair);
		        		if(VERBOSE){
		        			System.err.println("secondOfPair: " + secondOfPair + "\t" + "getMateAlignmentStart: " + samRecord.getMateAlignmentStart() + "\t" + "MateNegativeStrandFlag: " + samRecord.getMateNegativeStrandFlag());
			        		System.err.println("getReadString: " + samRecord.getReadString() + "\t" + "getAlignmentStart: " + samRecord.getAlignmentStart() + "\t" + "getUnclippedEnd: " + samRecord.getUnclippedEnd() + "\t" +"NegativeStrandFlag: " + samRecord.getReadNegativeStrandFlag() + "\tcytosineOffset: " + cytosineOffset + "\tbase: " + samRecord.getReadBases()[cytosineOffset] + "\t" + "baseQ: " + samRecord.getBaseQualities()[cytosineOffset] + "\tbase-1: " + samRecord.getReadBases()[cytosineOffset-1] + "\t" + "baseQ-1: " + samRecord.getBaseQualities()[cytosineOffset-1] + "\tbase+1: " + samRecord.getReadBases()[cytosineOffset+1] + "\t" + "baseQ+1: " + samRecord.getBaseQualities()[cytosineOffset+1]);

		        		}
		        			        		
		        	}
	        	}
	        	else{
	        		samRecord = p.getRead();
	        	}
				
				for(int i = 0; i < tmpKey[0].length(); i++){
				//	if(i == cytosinePos)
				//		continue;
					int elementOffset = i - cytosinePos + cytosineOffset;
					if(elementOffset < 0 || elementOffset > samRecord.getReadLength()-1)
						continue;
					//System.err.println("3: " + elementOffset + "\t" + cytosinePos + "\t" + cytosineOffset + "\t" + samRecord.getAlignmentStart() );
					PileupElement tmpP = new PileupElement(samRecord,elementOffset);
					if ( !usableBase(tmpP, true) ){
			        	continue;
			        }
					byte qual = tmpP.getQual();
					if ( qual > SAMUtils.MAX_PHRED_SCORE )
			            throw new UserException.MalformedBAM(tmpP.getRead(), String.format("the maximum allowed quality score is %d, but a quality of %d was observed in read %s.  Perhaps your BAM incorrectly encodes the quality scores in Sanger format; see http://en.wikipedia.org/wiki/FASTQ_format for more details", SAMUtils.MAX_PHRED_SCORE, qual, tmpP.getRead().getReadName()));
			        qual = (byte)Math.min((int)tmpP.getQual(), tmpP.getMappingQual());
			        byte observedBase = tmpP.getBase();
			        byte qualityScore = qual;
			        byte queryByte = (byte)tmpKey[0].charAt(i);
			        //queryByte = BaseUtilsMore.iupacCodeComplement(queryByte);
			        double error = pow(10, -qualityScore/10.0);
			        if ( observedBase == 0 )
			            continue;
			        adjacentCytosineSeqLikelihood[i] += BaseUtilsMore.iupacCodeEqual(queryByte, observedBase, samRecord.getReadNegativeStrandFlag()) ? log10(1.0 - error) : log10(error/3.0);
			        adjacentNotCytosineSeqLikelihood[i] += BaseUtilsMore.iupacCodeEqual(queryByte, observedBase, samRecord.getReadNegativeStrandFlag()) ? log10(error/3.0) : log10(1.0 - error);
			        if(VERBOSE){
			        	//System.err.println("forwardStrand-2: queryByte: " + (char)queryByte + "\tobservedBase: " + (char)observedBase + "\t" + BaseUtilsMore.iupacCodeEqual(queryByte, observedBase, true) + "\telementOffset: " + elementOffset + "\tcytosinePos: " + cytosinePos + "\tcytosineOffset: " + cytosineOffset + "\ti: " + i);
			        }
				}
				
				//int cytosineOffsetOpposite = samRecord.getReadLength() - p.getOffset();
				for(int i = 0; i < tmpKey[0].length(); i++){
				//	if(i == cytosinePos)
				//		continue;
					int elementOffset = cytosineOffset - (i - cytosinePos);
					if(elementOffset < 0 || elementOffset > samRecord.getReadLength()-1)
						continue;
					//System.err.println("4: " + elementOffset + "\t" + cytosinePos + "\t" + cytosineOffset + "\t" + samRecord.getAlignmentStart() );
					PileupElement tmpP = new PileupElement(samRecord,elementOffset);
					//System.err.println("4: " + tmpP.toString() );
					if ( !usableBase(tmpP, true) ){
			        	continue;
			        }
					byte qual = tmpP.getQual();
					if ( qual > SAMUtils.MAX_PHRED_SCORE )
			            throw new UserException.MalformedBAM(tmpP.getRead(), String.format("the maximum allowed quality score is %d, but a quality of %d was observed in read %s.  Perhaps your BAM incorrectly encodes the quality scores in Sanger format; see http://en.wikipedia.org/wiki/FASTQ_format for more details", SAMUtils.MAX_PHRED_SCORE, qual, tmpP.getRead().getReadName()));
			        qual = (byte)Math.min((int)tmpP.getQual(), tmpP.getMappingQual());
			        byte observedBase = tmpP.getBase();
			        byte qualityScore = qual;
			        byte queryByte = (byte)tmpKey[0].charAt(i);
			        queryByte = BaseUtilsMore.iupacCodeComplement(queryByte);
			        double error = pow(10, -qualityScore/10.0);
			        if ( observedBase == 0 )
			            continue;
			        adjacentCytosineSeqLikelihoodReverseStrand[i] += BaseUtilsMore.iupacCodeEqual(queryByte, observedBase, samRecord.getReadNegativeStrandFlag()) ? log10(1.0 - error) : log10(error/3.0);
			        adjacentNotCytosineSeqLikelihoodReverseStrand[i] += BaseUtilsMore.iupacCodeEqual(queryByte, observedBase, samRecord.getReadNegativeStrandFlag()) ? log10(error/3.0) : log10(1.0 - error);
			        if(VERBOSE){
			        	//System.err.println("reverseStrand-2: queryByte: " + (char)queryByte + "\tobservedBase: " + (char)observedBase + "\t" + BaseUtilsMore.iupacCodeEqual(queryByte, observedBase, true) + "\telementOffset: " + elementOffset + "\tcytosinePos: " + cytosinePos + "\tcytosineOffset: " + cytosineOffset + "\ti: " + i + "\tadjacentCytosineSeqLikelihoodReverseStrand[i]" + adjacentCytosineSeqLikelihoodReverseStrand[i]);
			        }
				}
				
	        }
	     
	        double sum = 0;
	        double sumReverse = 0;
	        double sumNot = 0;
	        double sumNotReverse = 0;
	        
	        for(int i = 0; i < adjacentCytosineSeqLikelihood.length; i++)
	        	sum += adjacentCytosineSeqLikelihood[i];
	        for(int i = 0; i < adjacentCytosineSeqLikelihoodReverseStrand.length; i++)
	        	sumReverse += adjacentCytosineSeqLikelihoodReverseStrand[i];
	        
	        for(int i = 0; i < adjacentNotCytosineSeqLikelihood.length; i++)
	        	sumNot += adjacentNotCytosineSeqLikelihood[i];
	        for(int i = 0; i < adjacentNotCytosineSeqLikelihoodReverseStrand.length; i++)
	        	sumNotReverse += adjacentNotCytosineSeqLikelihoodReverseStrand[i];
	        
	      
	        Double[] value = cts.cytosineListMap.get(cytosineType);
	        //value[0] = Double.NEGATIVE_INFINITY;
	        //value[1] = Double.NEGATIVE_INFINITY;
	        if(tmpKey[0].equalsIgnoreCase("C")){
	        	if(sum - sumNot >= ((autoEstimateC && !secondIteration) ? threshold : threshold/10.0) && sum >  ((autoEstimateC && !secondIteration) ? log10(1 - ABSOLUTE_CYTOSINE_TYPE_C_LIKELIHOOD) : log10(ABSOLUTE_CYTOSINE_TYPE_C_LIKELIHOOD))){
	        		//DETERMINED_CYTOSINE_TYPE_NEGSTRAND = false;
	        		

	    	       // 
	        		if(DETERMINED_CYTOSINE_TYPE_POS.isEmpty())
	        			DETERMINED_CYTOSINE_TYPE_POS = cytosineType;

	    	       // useDefaulMethyLevel = false;
	    	        value[0] = sum - sumNot;
	    	        value[3] = 1.0;
	    	      //  System.err.println("pos");
	        	}
	        	if(sumReverse - sumNotReverse >= ((autoEstimateC && !secondIteration) ? threshold : threshold/10.0) && sumReverse > ((autoEstimateC && !secondIteration) ? log10(1 - ABSOLUTE_CYTOSINE_TYPE_C_LIKELIHOOD) : log10(ABSOLUTE_CYTOSINE_TYPE_C_LIKELIHOOD))){
	        	
	    	        //cts.cytosineListMap.put(cytosineType, value);
	        		if(DETERMINED_CYTOSINE_TYPE_NEG.isEmpty())
	        			DETERMINED_CYTOSINE_TYPE_NEG = cytosineType;


	    	        value[1] = sumReverse - sumNotReverse;
	    	        value[3] = 1.0;
	    	      //  System.err.println("neg");
	        	}
	        }
	        else{
	        	if(sum - sumNot >= ((autoEstimateC && !secondIteration) ? threshold : threshold/10.0) && sum >  ((autoEstimateC && !secondIteration) ? log10(1 - ABSOLUTE_CYTOSINE_TYPE_C_LIKELIHOOD) : log10(ABSOLUTE_CYTOSINE_TYPE_C_LIKELIHOOD))){
	        		//DETERMINED_CYTOSINE_TYPE_NEGSTRAND = false;
	        		if(sum - sumNot > DETERMINED_CYTOSINE_TYPE_C_LIKELIHOOD_POS ){
	        			DETERMINED_CYTOSINE_TYPE_C_LIKELIHOOD_POS = sum - sumNot;
	        			DETERMINED_CYTOSINE_TYPE_POS = cytosineType;
	        			DETERMINED_CYTOSINE_TYPE_C_METHY_POS = value[2];
	        		}
	    	        
	    	       // 
	    	        
	    	        
	    	       // useDefaulMethyLevel = false;
	    	        value[0] = sum - sumNot;
	    	        
	    	      
	    	        
	    	        
	    	      //  System.err.println("pos");
	        	}
	        	if(sumReverse - sumNotReverse >= ((autoEstimateC && !secondIteration) ? threshold : threshold/10.0) && sumReverse > ((autoEstimateC && !secondIteration) ? log10(1 - ABSOLUTE_CYTOSINE_TYPE_C_LIKELIHOOD) : log10(ABSOLUTE_CYTOSINE_TYPE_C_LIKELIHOOD))){
	        		//DETERMINED_CYTOSINE_TYPE_NEGSTRAND = true;
	        		//Double[] value = cts.cytosineListMap.get(cytosineType);
	        		
	    	       
	    	        
	    	        if(sumReverse - sumNotReverse > DETERMINED_CYTOSINE_TYPE_C_LIKELIHOOD_NEG ){
	    	        	 DETERMINED_CYTOSINE_TYPE_C_LIKELIHOOD_NEG = sumReverse - sumNotReverse;
	 	    	        //cts.cytosineListMap.put(cytosineType, value);
	 	    	        DETERMINED_CYTOSINE_TYPE_NEG = cytosineType;
	 	    	        DETERMINED_CYTOSINE_TYPE_C_METHY_NEG = value[2];
	        		}
	    	       // useDefaulMethyLevel = false;
	    	        value[1] = sumReverse - sumNotReverse;
	    	        
	    	      //  System.err.println("neg");
	        	}
	        } 
	        cts.cytosineListMap.put(cytosineType, value);
	        
	        	 if(VERBOSE){
	 	        	System.err.println("cytosineType: " + cytosineType + "\tsum: " + sum + "\tsumNot: " + sumNot + "\tsumReverse" + sumReverse + "\tsumNotReverse" + sumNotReverse + "\tthreshold: " + threshold);
	 	        	System.err.println("DETERMINED_CYTOSINE_TYPE_POS: " + DETERMINED_CYTOSINE_TYPE_POS + "\tDETERMINED_CYTOSINE_TYPE_C_METHY_POS: " + DETERMINED_CYTOSINE_TYPE_C_METHY_POS  + "\tDETERMINED_CYTOSINE_TYPE_C_LIKELIHOOD_POS: " + DETERMINED_CYTOSINE_TYPE_C_LIKELIHOOD_POS + "\tvalue[0]:" + value[0]);
	 				 System.err.println("DETERMINED_CYTOSINE_TYPE_NEG: " + DETERMINED_CYTOSINE_TYPE_NEG + "\tDETERMINED_CYTOSINE_TYPE_C_METHY_NEG: " + DETERMINED_CYTOSINE_TYPE_C_METHY_NEG + "\tDETERMINED_CYTOSINE_TYPE_C_LIKELIHOOD_NEG: " + DETERMINED_CYTOSINE_TYPE_C_LIKELIHOOD_NEG + "\tvalue[1]:" + value[1] + "\tvalue[2]:" + value[2] + "\tvalue[3]:" + value[3]);
	 	        }	
	        	
	    }
		
	    if(!DETERMINED_CYTOSINE_TYPE_POS.isEmpty()){
	    	Double[] value = cts.cytosineListMap.get(DETERMINED_CYTOSINE_TYPE_POS);
	    	value[3] = 1.0;
	    	cts.cytosineListMap.put(DETERMINED_CYTOSINE_TYPE_POS, value);
	    }
	    else if(!DETERMINED_CYTOSINE_TYPE_NEG.isEmpty()){
	    	Double[] value = cts.cytosineListMap.get(DETERMINED_CYTOSINE_TYPE_NEG);
	    	value[3] = 1.0;
	    	cts.cytosineListMap.put(DETERMINED_CYTOSINE_TYPE_NEG, value);
	    }
	    
		//if(useDefaulMethyLevel && (!autoEstimateC || secondIteration)){
		//	DETERMINED_CYTOSINE_TYPE_C_METHY_NEG = defaultMethy;
		//	DETERMINED_CYTOSINE_TYPE_C_METHY_POS = defaultMethy;
		//}
		if(VERBOSE){
			// System.out.println("DETERMINED_CYTOSINE_TYPE_POS: " + DETERMINED_CYTOSINE_TYPE_POS + "\tDETERMINED_CYTOSINE_TYPE_C_METHY_POS: " + DETERMINED_CYTOSINE_TYPE_C_METHY_POS);
			 //System.out.println("DETERMINED_CYTOSINE_TYPE_NEG: " + DETERMINED_CYTOSINE_TYPE_NEG + "\tDETERMINED_CYTOSINE_TYPE_C_METHY_NEG: " + DETERMINED_CYTOSINE_TYPE_C_METHY_NEG);
		}
		
		
		//return cts;
	}
	
	/*
	public int add(ReadBackedPileup pileup, boolean ignoreBadBases, boolean capBaseQualsAtMappingQual, byte refNextBase, byte refPreBase) {
        int n = 0;
        
        ReadBackedPileup pileupPositiveStrand = pileup.getPositiveStrandPileup();
        for(PileupElement p : pileupPositiveStrand )
        	n += add(p, ignoreBadBases, capBaseQualsAtMappingQual, false, refNextBase, refPreBase);
        
        ReadBackedPileup pileupNegativeStrand = pileup.getNegativeStrandPileup();
        for(PileupElement p : pileupNegativeStrand )
        	n += add(p, ignoreBadBases, capBaseQualsAtMappingQual, true, refNextBase, refPreBase);
 
        return n;
    }
	*/
	//@Override
	public int add(ReadBackedPileup pileup, boolean ignoreBadBases, boolean capBaseQualsAtMappingQual) {
        int n = 0;
        
        ReadBackedPileup pileupPositiveStrand = pileup.getPositiveStrandPileup();
        for(PileupElement p : pileupPositiveStrand )
        	n += add(p, ignoreBadBases, capBaseQualsAtMappingQual, false);
        
        ReadBackedPileup pileupNegativeStrand = pileup.getNegativeStrandPileup();
        for(PileupElement p : pileupNegativeStrand )
        	n += add(p, ignoreBadBases, capBaseQualsAtMappingQual, true);
 
        return n;
    }
	
	public int add(PileupElement p, boolean ignoreBadBases, boolean capBaseQualsAtMappingQual, boolean negStrand) {
        // TODO-- Right now we assume that there are at most 2 reads per fragment.  This assumption is fine
        // TODO--   given the current state of next-gen sequencing, but may need to be fixed in the future.
        // TODO--   However, when that happens, we'll need to be a lot smarter about the caching we do here.
        byte observedBase = 0, qualityScore = 0;
        if ( !usableBase(p, ignoreBadBases) ){
        	//System.out.println("bad base: " + p.getBase());
        	return 0;
        }
            
        SAMRecord samRecord = p.getRead();
        if(BAC.pairedEndMode){
        	try {
    			samRecord = (SAMRecord) p.getRead().clone();
    		} catch (CloneNotSupportedException e) {
    			// TODO Auto-generated catch block
    			e.printStackTrace();
    		}
        	
        	boolean Paired = samRecord.getReadPairedFlag();
        	boolean secondOfPair = samRecord.getSecondOfPairFlag();
        	if (samRecord.getNotPrimaryAlignmentFlag())
    		{
    			return 0;
    		}
    		
    		// Inverted dups, count only one end
    		if (samRecord.getAlignmentStart() == samRecord.getMateAlignmentStart() && samRecord.getReadNegativeStrandFlag() == samRecord.getMateNegativeStrandFlag())
    		{
    			if (samRecord.getSecondOfPairFlag()) return 0;
    			//System.err.printf("Inverted dup %d%s (%s)\n", samRecord.getAlignmentStart(), samRecord.getReadNegativeStrandFlag()?"-":"+", PicardUtils.getReadString(samRecord, true));
    		}
        	if (Paired  && !BAC.allowBadMates && !samRecord.getProperPairFlag())
    		{
        		return 0;
    		}
        	if(secondOfPair){
           	 //observedBase = BaseUtils.simpleComplement(p.getBase()); 
           	 negStrand = !negStrand;
        	}
        }
		
    	
    	
        byte qual = p.getQual();
        if ( qual > SAMUtils.MAX_PHRED_SCORE )
            throw new UserException.MalformedBAM(p.getRead(), String.format("the maximum allowed quality score is %d, but a quality of %d was observed in read %s.  Perhaps your BAM incorrectly encodes the quality scores in Sanger format; see http://en.wikipedia.org/wiki/FASTQ_format for more details", SAMUtils.MAX_PHRED_SCORE, qual, p.getRead().getReadName()));
        if ( capBaseQualsAtMappingQual )
            qual = (byte)Math.min((int)p.getQual(), p.getMappingQual());
        
        observedBase = p.getBase(); 
        
        qualityScore = qual;
         
        // abort early if we didn't see any good bases
        if ( observedBase == 0 )
            return 0;


        BisulfiteDiploidSNPGenotypeLikelihoods gl;
        gl = getCalculateGenotypeLikelihoods(observedBase, qualityScore, negStrand);


        // for bad bases, there are no likelihoods
        if ( gl == null )
            return 0;

        double[] likelihoods = gl.getLikelihoods();

        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            double likelihood = likelihoods[g.ordinal()];
            
            //if ( VERBOSE ) {
            //    System.out.printf("  L(%c | G=%s, Q=%d, S=%s) = %f / %f%n",
            //            observedBase, g, qualityScore, pow(10,likelihood) * 100, likelihood);
            //}

            log10Likelihoods[g.ordinal()] += likelihood;
            log10Posteriors[g.ordinal()] += likelihood;
        }

        return 1;
    }
/*
	public int add(PileupElement p, boolean ignoreBadBases, boolean capBaseQualsAtMappingQual, boolean negStrand, byte refNextBase, byte refPreBase) {
        // TODO-- Right now we assume that there are at most 2 reads per fragment.  This assumption is fine
        // TODO--   given the current state of next-gen sequencing, but may need to be fixed in the future.
        // TODO--   However, when that happens, we'll need to be a lot smarter about the caching we do here.
        byte observedBase = 0, qualityScore = 0;
        if ( !usableBase(p, ignoreBadBases) ){
        	//System.out.println("bad base: " + p.getBase());
        	return 0;
        }
            	
        byte qual = p.getQual();
        if ( qual > SAMUtils.MAX_PHRED_SCORE )
            throw new UserException.MalformedBam(p.getRead(), String.format("the maximum allowed quality score is %d, but a quality of %d was observed in read %s.  Perhaps your BAM incorrectly encodes the quality scores in Sanger format; see http://en.wikipedia.org/wiki/FASTQ_format for more details", SAMUtils.MAX_PHRED_SCORE, qual, p.getRead().getReadName()));
        if ( capBaseQualsAtMappingQual )
            qual = (byte)Math.min((int)p.getQual(), p.getMappingQual());
         
        observedBase = p.getBase();
        qualityScore = qual;
         
        // abort early if we didn't see any good bases
        if ( observedBase == 0 )
            return 0;


        BisulfiteDiploidSNPGenotypeLikelihoods gl;
        gl = getCalculateGenotypeLikelihoods(observedBase, qualityScore, negStrand, refNextBase, refPreBase);


        // for bad bases, there are no likelihoods
        if ( gl == null )
            return 0;

        double[] likelihoods = gl.getLikelihoods();

        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            double likelihood = likelihoods[g.ordinal()];
            
            //if ( VERBOSE ) {
            //    System.out.printf("  L(%c | G=%s, Q=%d, S=%s) = %f / %f%n",
            //            observedBase, g, qualityScore, pow(10,likelihood) * 100, likelihood);
            //}

            log10Likelihoods[g.ordinal()] += likelihood;
            log10Posteriors[g.ordinal()] += likelihood;
        }

        return 1;
    }
	
*/	
	/*
	 * for reads: AAGGCCTTaaggcctt (Q=30)
	 * likelihood of all of the genotype would be (err = pow(10,Q/-10))
	 * AA: (1-err)^2 * (err/3)^10 * (1-err+err*conv/3)^2 * (err/3 * (1-conv))^2
	 * AG: [1/2 * (1-err) + 1/2 * err/3]^4 *  [1/2 * err/3 + 1/2 * (1-err)]^2 * (err/3)^8 * [1/2 * (1-err) + 1/2 * (err/3) * conv]
	 * AC: [1/2 * (1-err) + 1/2 * err/3 * (1-conv)]^2 * (err/3)^2 * [1/2 * (err/3) * (1-conv) + 1/2 * (1-err)* (1-conv)]^2 * [1/2 * err/3 + 1/2 * (err/3 + conv)]^2 * [1/2 * (1-err) + 1/2 * (err/3 + conv)]^2 * (err/3)^2 * [1/2 * err/3 + 1/2 * (1-err)]^2 * (err/3)^2 
	 * AT: 
	 * 
	 * 
	*/
	/*
	protected BisulfiteDiploidSNPGenotypeLikelihoods getCalculateGenotypeLikelihoods(byte observedBase, byte qualityScore, boolean negStrand, byte refNextBase, byte refPreBase){

		try {
			BisulfiteDiploidSNPGenotypeLikelihoods gl = (BisulfiteDiploidSNPGenotypeLikelihoods)this.clone();
			gl.setToZeroBs();
			double error = pow(10, -qualityScore/10.0);
			double pOfBase1 = 0.5, pOfBase2 = 1.0 - pOfBase1;
			for ( DiploidGenotype g : DiploidGenotype.values() ) {
				double likelihood = 0.0;
				if(!negStrand){
						switch(observedBase){
							case BaseUtils.C:
								if(refNextBase == BaseUtils.G){
									likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) * ( (1.0-CPG_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPG_METHYLATION_RATE ) : pOfBase1 * (error/3.0) * ( (1.0-CPG_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPG_METHYLATION_RATE );
									likelihood += observedBase == g.base2 ? pOfBase2 * (1.0-error) *  ( (1.0-CPG_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPG_METHYLATION_RATE ) : pOfBase2 * (error/3.0) * ( (1.0-CPG_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPG_METHYLATION_RATE );
									if ( VERBOSE ) {
									//	System.out.println("flag1: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
									}
								}
								else{
									likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) * ((1.0-CPH_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPH_METHYLATION_RATE ) : pOfBase1 * (error/3.0) * ((1.0-CPH_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPH_METHYLATION_RATE );
									likelihood += observedBase == g.base2 ? pOfBase2 * (1.0-error) * ((1.0-CPH_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPH_METHYLATION_RATE ) : pOfBase2 * (error/3.0) * ((1.0-CPH_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPH_METHYLATION_RATE );
									if ( VERBOSE ) {
									//	System.out.println("flag2: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
									}
								}
								break;
							case BaseUtils.T:
								if(refNextBase == BaseUtils.G){
									likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) : (g.base1 == BaseUtils.C ? pOfBase1 * (error/3.0 + (1.0-CPG_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE) : pOfBase1 * (error/3.0));
									likelihood += observedBase == g.base2 ? pOfBase1 * (1.0-error) : (g.base2 == BaseUtils.C ? pOfBase2 * (error/3.0 + (1.0-CPG_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE) : pOfBase2 * (error/3.0));
									if ( VERBOSE ) {
									//	System.out.println("flag3: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
									}
								}
								else{
									likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) : (g.base1 == BaseUtils.C ? pOfBase1 * (error/3.0 + (1.0-CPH_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE) : pOfBase1 * (error/3.0));
									likelihood += observedBase == g.base2 ? pOfBase1 * (1.0-error) : (g.base2 == BaseUtils.C ? pOfBase2 * (error/3.0 + (1.0-CPH_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE) : pOfBase2 * (error/3.0));
									if ( VERBOSE ) {
									//	System.out.println("flag4: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
									}
								}
								
								break;
							default: 
								likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) : pOfBase1 * (error/3.0);
								likelihood += observedBase == g.base2 ? pOfBase2 * (1.0-error) : pOfBase2 * (error/3.0);
								if ( VERBOSE ) {
									//System.out.println("flag5: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
								}
									
						}
				}
				else{
						switch(observedBase){
						case BaseUtils.G:
							if(refPreBase == BaseUtils.C){
								likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) * ( (1.0-CPG_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPG_METHYLATION_RATE ) : pOfBase1 * (error/3.0) * ( (1.0-CPG_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPG_METHYLATION_RATE );
								likelihood += observedBase == g.base2 ? pOfBase2 * (1.0-error) *  ( (1.0-CPG_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPG_METHYLATION_RATE ) : pOfBase2 * (error/3.0) * ( (1.0-CPG_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPG_METHYLATION_RATE );
								if ( VERBOSE ) {
									//System.out.println("flag6: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
								}
							}
							else{
								likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) * ((1.0-CPH_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPH_METHYLATION_RATE) : pOfBase1 * (error/3.0) * ((1.0-CPH_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPH_METHYLATION_RATE);
								likelihood += observedBase == g.base2 ? pOfBase2 * (1.0-error) * ((1.0-CPH_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPH_METHYLATION_RATE) : pOfBase2 * (error/3.0) * ((1.0-CPH_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPH_METHYLATION_RATE);
								if ( VERBOSE ) {
									//System.out.println("flag7: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
								}
							}
							break;
						case BaseUtils.A:
							if(refPreBase == BaseUtils.C){
								likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) : (g.base1 == BaseUtils.G ? pOfBase1 * (error/3.0 + (1.0-CPG_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE) : pOfBase1 * (error/3.0));
								likelihood += observedBase == g.base2 ? pOfBase1 * (1.0-error) : (g.base2 == BaseUtils.G ? pOfBase2 * (error/3.0 + (1.0-CPG_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE) : pOfBase2 * (error/3.0));
								if ( VERBOSE ) {
									//System.out.println("flag8: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
								}
							}
							else{
								likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) : (g.base1 == BaseUtils.G ? pOfBase1 * (error/3.0 + (1.0-CPH_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE) : pOfBase1 * (error/3.0));
								likelihood += observedBase == g.base2 ? pOfBase1 * (1.0-error) : (g.base2 == BaseUtils.G ? pOfBase2 * (error/3.0 + (1.0-CPH_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE) : pOfBase2 * (error/3.0));
								if ( VERBOSE ) {
									//System.out.println("flag9: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
								}
							}
							
							break;
						default: 
							likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) : pOfBase1 * (error/3.0);
							likelihood += observedBase == g.base2 ? pOfBase2 * (1.0-error) : pOfBase2 * (error/3.0);
							if ( VERBOSE ) {
								//System.out.println("flag9: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
							}
						}
				}
				 double logLikelihood = log10(likelihood);
	             gl.log10Likelihoods[g.ordinal()] += logLikelihood;
	             gl.log10Posteriors[g.ordinal()] += logLikelihood;
			}
			
			if ( VERBOSE ) {
                for ( DiploidGenotype g : DiploidGenotype.values() ) { System.out.printf("%s\t", g); }
                System.out.println();
                for ( DiploidGenotype g : DiploidGenotype.values() ) { System.out.printf("%.2f\t", gl.log10Likelihoods[g.ordinal()]); }
                System.out.println();
                for ( DiploidGenotype g : DiploidGenotype.values() ) { System.out.printf("%.2f\t", gl.getPriors()[g.ordinal()]); }
                System.out.println();
            }

            return gl;
            
		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			throw new RuntimeException(e);
		}
		
		
		
	}
	*/

	protected BisulfiteDiploidSNPGenotypeLikelihoods getCalculateGenotypeLikelihoods(byte observedBase, byte qualityScore, boolean negStrand){

		try {
			BisulfiteDiploidSNPGenotypeLikelihoods gl = (BisulfiteDiploidSNPGenotypeLikelihoods)this.clone();
			gl.setToZeroBs();
			double error = pow(10, -qualityScore/10.0);
			double pOfBase1 = 0.5, pOfBase2 = 1.0 - pOfBase1;
			for ( DiploidGenotype g : DiploidGenotype.values() ) {
				double likelihood = 0.0;
				if(!negStrand){
						switch(observedBase){
							case BaseUtils.C:
						//		if(!DETERMINED_CYTOSINE_TYPE_NEGSTRAND){
									likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) * ( (1.0-DETERMINED_CYTOSINE_TYPE_C_METHY_POS) * (1.0-BISULFITE_CONVERSION_RATE) + DETERMINED_CYTOSINE_TYPE_C_METHY_POS * (1.0-OVER_CONVERSION_RATE)) : pOfBase1 * (error/3.0);
									likelihood += observedBase == g.base2 ? pOfBase2 * (1.0-error) *  ( (1.0-DETERMINED_CYTOSINE_TYPE_C_METHY_POS) * (1.0-BISULFITE_CONVERSION_RATE) + DETERMINED_CYTOSINE_TYPE_C_METHY_POS * (1.0-OVER_CONVERSION_RATE) ) : pOfBase2 * (error/3.0);
							//	}
						//		else{
							//		likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) * ( (1.0-UNDETERMINED_CYTOSINE_TYPE_C_METHY) * (1.0-BISULFITE_CONVERSION_RATE) + UNDETERMINED_CYTOSINE_TYPE_C_METHY ) : pOfBase1 * (error/3.0) * ( (1.0-UNDETERMINED_CYTOSINE_TYPE_C_METHY) * (1.0-BISULFITE_CONVERSION_RATE) + UNDETERMINED_CYTOSINE_TYPE_C_METHY );
								//	likelihood += observedBase == g.base2 ? pOfBase2 * (1.0-error) *  ( (1.0-UNDETERMINED_CYTOSINE_TYPE_C_METHY) * (1.0-BISULFITE_CONVERSION_RATE) + UNDETERMINED_CYTOSINE_TYPE_C_METHY ) : pOfBase2 * (error/3.0) * ( (1.0-UNDETERMINED_CYTOSINE_TYPE_C_METHY) * (1.0-BISULFITE_CONVERSION_RATE) + UNDETERMINED_CYTOSINE_TYPE_C_METHY );
								//}
									
									if ( VERBOSE ) {
									//	System.out.println("flag1: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
									}
								
								break;
							case BaseUtils.T:
							
									likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) : (g.base1 == BaseUtils.C ? pOfBase1 * (error/3.0 + (1.0-DETERMINED_CYTOSINE_TYPE_C_METHY_POS) * BISULFITE_CONVERSION_RATE + DETERMINED_CYTOSINE_TYPE_C_METHY_POS * OVER_CONVERSION_RATE) : pOfBase1 * (error/3.0));
									likelihood += observedBase == g.base2 ? pOfBase1 * (1.0-error) : (g.base2 == BaseUtils.C ? pOfBase2 * (error/3.0 + (1.0-DETERMINED_CYTOSINE_TYPE_C_METHY_POS) * BISULFITE_CONVERSION_RATE + DETERMINED_CYTOSINE_TYPE_C_METHY_POS * OVER_CONVERSION_RATE) : pOfBase2 * (error/3.0));

								if ( VERBOSE ) {
									//	System.out.println("flag3: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
									}
								
								break;
							default: 
								likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) : pOfBase1 * (error/3.0);
								likelihood += observedBase == g.base2 ? pOfBase2 * (1.0-error) : pOfBase2 * (error/3.0);
								if ( VERBOSE ) {
									//System.out.println("flag5: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
								}
									
						}
				}
				else{
						switch(observedBase){
						case BaseUtils.G:
						//	if(DETERMINED_CYTOSINE_TYPE_NEGSTRAND){
								likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) * ( (1.0-DETERMINED_CYTOSINE_TYPE_C_METHY_NEG) * (1.0-BISULFITE_CONVERSION_RATE) + DETERMINED_CYTOSINE_TYPE_C_METHY_NEG * (1.0-OVER_CONVERSION_RATE) ) : pOfBase1 * (error/3.0);
								likelihood += observedBase == g.base2 ? pOfBase2 * (1.0-error) *  ( (1.0-DETERMINED_CYTOSINE_TYPE_C_METHY_NEG) * (1.0-BISULFITE_CONVERSION_RATE) + DETERMINED_CYTOSINE_TYPE_C_METHY_NEG * (1.0-OVER_CONVERSION_RATE) ) : pOfBase2 * (error/3.0);
								if ( VERBOSE ) {
									//System.out.println("flag6: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
								}
					//		}
					//		else{
					//			likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) * ((1.0-UNDETERMINED_CYTOSINE_TYPE_C_METHY) * (1.0-BISULFITE_CONVERSION_RATE) + UNDETERMINED_CYTOSINE_TYPE_C_METHY) : pOfBase1 * (error/3.0) * ((1.0-UNDETERMINED_CYTOSINE_TYPE_C_METHY) * (1.0-BISULFITE_CONVERSION_RATE) + CPH_METHYLATION_RATE);
					//			likelihood += observedBase == g.base2 ? pOfBase2 * (1.0-error) * ((1.0-UNDETERMINED_CYTOSINE_TYPE_C_METHY) * (1.0-BISULFITE_CONVERSION_RATE) + UNDETERMINED_CYTOSINE_TYPE_C_METHY) : pOfBase2 * (error/3.0) * ((1.0-UNDETERMINED_CYTOSINE_TYPE_C_METHY) * (1.0-BISULFITE_CONVERSION_RATE) + CPH_METHYLATION_RATE);
					//			if ( VERBOSE ) {
									//System.out.println("flag7: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
					//			}
							//}
							break;
						case BaseUtils.A:
							
								likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) : (g.base1 == BaseUtils.G ? pOfBase1 * (error/3.0 + (1.0-DETERMINED_CYTOSINE_TYPE_C_METHY_NEG) * BISULFITE_CONVERSION_RATE + DETERMINED_CYTOSINE_TYPE_C_METHY_NEG * OVER_CONVERSION_RATE) : pOfBase1 * (error/3.0));
								likelihood += observedBase == g.base2 ? pOfBase1 * (1.0-error) : (g.base2 == BaseUtils.G ? pOfBase2 * (error/3.0 + (1.0-DETERMINED_CYTOSINE_TYPE_C_METHY_NEG) * BISULFITE_CONVERSION_RATE + DETERMINED_CYTOSINE_TYPE_C_METHY_NEG * OVER_CONVERSION_RATE) : pOfBase2 * (error/3.0));
								if ( VERBOSE ) {
									//System.out.println("flag9: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
								}
							
							
							break;
						default: 
							likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) : pOfBase1 * (error/3.0);
							likelihood += observedBase == g.base2 ? pOfBase2 * (1.0-error) : pOfBase2 * (error/3.0);
							if ( VERBOSE ) {
								//System.out.println("flag9: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
							}
						}
				}
				 double logLikelihood = log10(likelihood);
	             gl.log10Likelihoods[g.ordinal()] += logLikelihood;
	             gl.log10Posteriors[g.ordinal()] += logLikelihood;
			}
			
			if ( VERBOSE ) {
                for ( DiploidGenotype g : DiploidGenotype.values() ) { System.out.printf("%s\t", g); }
                System.out.println();
                for ( DiploidGenotype g : DiploidGenotype.values() ) { System.out.printf("%.2f\t", gl.log10Likelihoods[g.ordinal()]); }
                System.out.println();
                for ( DiploidGenotype g : DiploidGenotype.values() ) { System.out.printf("%.2f\t", gl.getPriors()[g.ordinal()]); }
                System.out.println();
            }

            return gl;
            
		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			throw new RuntimeException(e);
		}
		
		
		
	}
	
	public double[] getLikelihoods() {
        return log10Likelihoods;
    }
	
	public double[] getPosteriors() {
        return log10Posteriors;
    }
	
	protected static boolean usableBase(PileupElement p, boolean ignoreBadBases) {
        // ignore deletions, Q0 bases, and filtered bases
        if ( p.isDeletion() ||
                p.getQual() == 0 ||
                (p.getRead() instanceof GATKSAMRecord &&
                 !((GATKSAMRecord)p.getRead()).isGoodBase(p.getOffset())) )
            return false;

        return ( !ignoreBadBases || !badBase(p.getBase()) );
    }

    /**
     * Returns true when the observedBase is considered bad and shouldn't be processed by this object.  A base
     * is considered bad if:
     *
     *   Criterion 1: observed base isn't a A,C,T,G or lower case equivalent
     *
     * @param observedBase observed base
     * @return true if the base is a bad base
     */
    protected static boolean badBase(byte observedBase) {
        return BaseUtils.simpleBaseToBaseIndex(observedBase) == -1;
    }
	
    //@Override
    public double[] getPriors() {
        return this.priors.getPriors();
    }
    
    public void setPriors(BisulfiteDiploidSNPGenotypePriors priors) {
        this.priors = priors;
        //log10Posteriors = new double[DiploidGenotype.values().length];
        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            int i = g.ordinal();
            log10Posteriors[i] = priors.getPriors()[i] + log10Likelihoods[i];
        }
    }
 // likelihoods are all zeros

    protected void setToZeroBs() {
    	//System.err.println("test!!!");
    	log10Likelihoods = new double[DiploidGenotype.values().length];
    	for ( DiploidGenotype g : DiploidGenotype.values() ) {
            int i = g.ordinal();
            log10Likelihoods[i] = log10(1-PROB_OF_REFERENCE_ERROR);
        }

       // if(this.priors.getPriors() == null){
        //	System.err.println("priors null!!!");
       // }
       // System.err.println("log10Likelihoods null!!!");
        log10Posteriors = this.priors.getPriors().clone();     // posteriors are all the priors
    }

   // @Override
    protected Object clone() throws CloneNotSupportedException {
    	BisulfiteDiploidSNPGenotypeLikelihoods c = (BisulfiteDiploidSNPGenotypeLikelihoods)super.clone();
        c.priors = priors;
        c.log10Likelihoods = log10Likelihoods.clone();
        c.log10Posteriors = log10Posteriors.clone();
        return c;
    }
    
    /**
     * Sets the priors
     * @param priors priors
     */
	
}
