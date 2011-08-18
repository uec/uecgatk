package edu.usc.epigenome.uecgatk.bisulfitesnpmodel;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import net.sf.samtools.SAMRecord;

import org.apache.log4j.Logger;
import org.broad.tribble.Feature;
import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel.BisulfiteDiploidSNPGenotypeLikelihoods;
import org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel.BisulfiteDiploidSNPGenotypePriors;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext.StratifiedContextType;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.helpers.DbSNPHelper;
import org.broadinstitute.sting.gatk.walkers.genotyper.BiallelicGenotypeLikelihoods;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypePriors;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.genotype.DiploidGenotype;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;



public class BisulfiteSNPGenotypeLikelihoodsCalculationModel extends
		NonRefDependSNPGenotypeLikelihoodsCalculationModel {

	private BisulfiteArgumentCollection BAC;
	
	protected Byte bestAllele = null;
	protected Byte alternateAllele = null;
	//protected final boolean useAlleleFromVCF;
	protected long testLoc;
	//protected boolean isCGI = false;
	//public byte[] CONTEXTREF = null;
	protected static Integer numCNegStrand = 0;
	protected static Integer numTNegStrand = 0;
	protected static Integer numCPosStrand = 0;
	protected static Integer numTPosStrand = 0;
	//protected static String CYTOSINE_TYPE_LIST = null;
	//protected String DETERMINED_CYTOSINE_TYPE = "C";
	private CytosineTypeStatus cts = null;
	//protected static Boolean cytosineStrand = false;
	//protected static String cytosineWindowContext = "C";
	private static boolean autoEstimateC = false;
    private static boolean secondIteration = false;
    private double FALT_METHY_STATUS = 0.5;
	

	
	
	public BisulfiteSNPGenotypeLikelihoodsCalculationModel(
			UnifiedArgumentCollection UAC, Logger logger) {
		super(UAC, logger);
		//useAlleleFromVCF = UAC.GenotypingMode == GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES;
		
		// TODO Auto-generated constructor stub
	}
	
	@Override
	public void initialize(CytosineTypeStatus cts, BisulfiteArgumentCollection BAC, boolean autoEstimateC, boolean secondIteration){
		
		this.BAC = BAC;
		this.cts = cts;
		this.testLoc = BAC.testLocus;
		this.autoEstimateC = autoEstimateC;
		this.secondIteration = secondIteration;
		if(BAC.sequencingMode == MethylSNPModel.NM){
			FALT_METHY_STATUS = 0.0;
		}
			
	//	CYTOSINE_TYPE_LIST = BAC.cytosineType;

	}
	
	
	@Override
	public Allele getLikelihoods(RefMetaDataTracker tracker,
			ReferenceContext ref,
			Map<String, StratifiedAlignmentContext> contexts,
			StratifiedContextType contextType, GenotypePriors priors,
			Map<String, BiallelicGenotypeLikelihoods> GLs,
			Allele alternateAlleleToUse) {
		if ( !(priors instanceof BisulfiteDiploidSNPGenotypePriors) )
            throw new StingException("Only Bisulfite diploid-based SNP priors are supported in the BSSNP GL model");

		byte refBase = ref.getBase();
		
		
        Allele refAllele = Allele.create(refBase, true);
       
        numCNegStrand = 0;
        numTNegStrand = 0;
        numCPosStrand = 0;
        numTPosStrand = 0;
        

        for ( Map.Entry<String, StratifiedAlignmentContext> sample : contexts.entrySet() ) {
            ReadBackedPileup pileup = sample.getValue().getContext(contextType).getBasePileup();
            for ( PileupElement p : pileup ) {
            	SAMRecord samRecord = p.getRead();
            	int offset = p.getOffset();
            	if(offset < 0)//is deletion
            		continue;
            	if(BAC.pairedEndMode){
            		try {
    					samRecord = (SAMRecord) p.getRead().clone();
    				} catch (CloneNotSupportedException e) {
    					// TODO Auto-generated catch block
    					e.printStackTrace();
    				}
                	boolean Paired = samRecord.getReadPairedFlag();	
    	        	boolean secondOfPair = samRecord.getSecondOfPairFlag();
    	        	
    	        	//byte base = p.getBase();
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
    	        	
    	        	if((pileup.getLocation().getStart()) == testLoc){
            			System.err.println("NegativeStrandFlag: " + samRecord.getReadNegativeStrandFlag() + "\t" + "MateNegativeStrandFlag: " + samRecord.getMateNegativeStrandFlag() + "\tbase: " + samRecord.getReadBases()[offset] + "\t" + "baseQ: " + samRecord.getBaseQualities()[offset]);
            			System.err.println("getReadString: " + samRecord.getReadString() + "\tsecond: " + secondOfPair);
            		}
    	        	if(secondOfPair){	        		
    		        	samRecord.setReadNegativeStrandFlag(!samRecord.getReadNegativeStrandFlag());        		
    	        	}
    	        	if((pileup.getLocation().getStart()) == testLoc){
            			System.err.println("proper paired: " + samRecord.getProperPairFlag() + "\t" + "getMateAlignmentStart: " + samRecord.getMateAlignmentStart() + "\t" + "MateNegativeStrandFlag: " + samRecord.getMateNegativeStrandFlag());
    	        		System.err.println("getReadString: " + samRecord.getReadString() + "\t" + "getAlignmentStart: " + samRecord.getAlignmentStart() + "\t" + "getUnclippedEnd: " + samRecord.getUnclippedEnd() + "\t" +"NegativeStrandFlag: " + samRecord.getReadNegativeStrandFlag() + "\tcytosineOffset: " + offset + "\tbase: " + samRecord.getReadBases()[offset] + "\t" + "baseQ: " + samRecord.getBaseQualities()[offset]);
    	        		System.err.println("getBase: " + p.getBase() + "\tp.getRead(): " + p.getRead().getReadBases()[offset]);
            		}
            	}
				
	        	boolean negStrand = samRecord.getReadNegativeStrandFlag();
				int alignmentS = samRecord.getAlignmentStart();
				int	onRefCoord = (negStrand) ? samRecord.getUnclippedEnd() : alignmentS;
				
				if(((GATKSAMRecord)p.getRead()).isGoodBase(offset)){
					if(negStrand){
						if(p.getBase()==BaseUtils.G){
							numCNegStrand++;
						}
						else if(p.getBase()==BaseUtils.A){
							numTNegStrand++;
						}
						else{
							
						}
						
					}
					else{
						if(p.getBase()==BaseUtils.C){
							numCPosStrand++;
						}
						else if(p.getBase()==BaseUtils.T){
							numTPosStrand++;
						}
						else{
							
						}
					}
				}
				
									
				if((pileup.getLocation().getStart()) == testLoc){
					System.err.println("before filter:\t" + onRefCoord + "\t" + offset + "\t" + negStrand + "\t" + pileup.getLocation().getStart() + "\t" + (char)p.getBase());
					System.err.println("refBase: " + refBase);
					
					if(BAC.pairedEndMode)
						System.err.println("isGoodBase: " + ((GATKSAMRecord)p.getRead()).isGoodBase(offset) + "\tsecondOfPair: " + "\tchanged: " + samRecord.getSecondOfPairFlag());
		                     
				}
            }
            
           // double[] cytosineMethyStatus = new double[2];
           // cytosineMethyStatus[0]=0.5;
          //  cytosineMethyStatus[1]=0.5;
            BisulfiteDiploidSNPGenotypeLikelihoods GL = checkCytosineStatus(pileup, cts, BAC.cTypeThreshold, tracker, ref, (BisulfiteDiploidSNPGenotypePriors)priors);
            
            if(GL == null)
            	return refAllele;
            
         //   BisulfiteDiploidSNPGenotypeLikelihoods GL = new BisulfiteDiploidSNPGenotypeLikelihoods(tracker, ref, (BisulfiteDiploidSNPGenotypePriors)priors, BAC, cytosineMethyStatus.clone());
           // if((pileup.getLocation().getStart()) == testLoc)
            //	GL.VERBOSE=true;
            if((pileup.getLocation().getStart()) == testLoc){
     			// System.err.println("CYTOSINE_STATUS[0]: " + CYTOSINE_STATUS[0] + "\tCYTOSINE_STATUS[1]: " + CYTOSINE_STATUS[1] + "\tCYTOSINE_STATUS[2]: " + CYTOSINE_STATUS[2] + "\tCYTOSINE_STATUS[3]: " + CYTOSINE_STATUS[3]);
     			 for(String cytosineType : cts.cytosineListMap.keySet()){
     					String[] tmpKey = cytosineType.split("-");
     					Double[] value = cts.cytosineListMap.get(cytosineType);
     					System.err.println("tmpKey[0]" + tmpKey[0] + "\tvalue[0]" + value[0] + "\tvalue[1]" + value[1] + "\tvalue[2]" + value[2] + "\tvalue[3]" + value[3]);
     			 }
     			
     		 }
            
           // GL.setPriors(tracker, ref, BAC.PCR_error, BAC.novelDbsnpHet, BAC.validateDbsnpHet);
            
           // int nGoodBases = GL.add(pileup, true, true);
           // if ( nGoodBases == 0 )
             //   continue;

            double[] likelihoods_befor = GL.getLikelihoods();
            double[] posterior_befor = GL.getPosteriors();
            double[] prio = GL.getPriors();
           // double[] likelihoods = normalization(likelihoods_befor,likelihoods_befor);
          //  double[] posterior = normalization(posterior_befor,likelihoods_befor);
            double[] likelihoods = GL.getLikelihoods();
            double[] posterior = GL.getPosteriors();
            
            
            initializeBestAndAlternateAlleleFromPosterior(posterior, pileup.getLocation().getStart());
            
            if ( (alternateAllele == null && bestAllele == refBase) || (bestAllele == null) ) {
               
                if ( BAC.OutputMode == BisulfiteGenotyperEngine.OUTPUT_MODE.EMIT_VARIANTS_ONLY )
                    return refAllele;

            }
            
            Allele AlleleA, AlleleB;
            
            
            if(alternateAllele == null || BaseUtils.basesAreEqual(alternateAllele,refBase) || alternateAllele == bestAllele){
            	AlleleA = Allele.create(refBase, true);
            	AlleleB = Allele.create(bestAllele, false);
            	
            		alternateAllele = bestAllele;
            	bestAllele = refBase;
            	
            	
            }
            else if(BaseUtils.basesAreEqual(bestAllele,refBase)){
            	AlleleA = Allele.create(bestAllele, true);
            	AlleleB = Allele.create(alternateAllele, false);
            	
            }
            else{
            	AlleleA = Allele.create(bestAllele, false);
            	AlleleB = Allele.create(alternateAllele, false);
            	if(AlleleA.equals(refAllele, true)){
            		AlleleA = Allele.create(bestAllele, true);
            	
            	}
            	
            	if(AlleleB.equals(refAllele, true)){
            		AlleleB = Allele.create(alternateAllele, true);
            		
            	}
            	
            	
            }
            DiploidGenotype AAGenotype = DiploidGenotype.createHomGenotype(bestAllele);
            DiploidGenotype ABGenotype = DiploidGenotype.createDiploidGenotype(bestAllele, alternateAllele);
            DiploidGenotype BBGenotype = DiploidGenotype.createHomGenotype(alternateAllele);
            
            
            if((pileup.getLocation().getStart()) == testLoc){
            	System.err.println("sample: " + sample.getKey());
            	System.err.println("sample location: " + pileup.getPileupString((char)refBase));
            	System.err.println("sample: " + sample.getValue().getLocation().getStart());
            	System.err.println("refBase: " + refBase + " bestAllele: " + bestAllele + " alternateAllele: " + alternateAllele);
            
            	System.err.println("AAGenotype " + likelihoods[AAGenotype.ordinal()] + "\t" + prio[AAGenotype.ordinal()] + "\t" + posterior[AAGenotype.ordinal()]);
            	System.err.println("ABGenotype " + likelihoods[ABGenotype.ordinal()] + "\t" + prio[ABGenotype.ordinal()] + "\t" + posterior[ABGenotype.ordinal()]);
            	System.err.println("BBGenotype " + likelihoods[BBGenotype.ordinal()] + "\t" + prio[BBGenotype.ordinal()] + "\t" + posterior[BBGenotype.ordinal()]);
            	System.err.println("AAGenotype before normalize " + likelihoods_befor[AAGenotype.ordinal()] + "\t" + prio[AAGenotype.ordinal()] + "\t" + posterior_befor[AAGenotype.ordinal()]);
            	System.err.println("ABGenotype before normaliz " + likelihoods_befor[ABGenotype.ordinal()] + "\t" + prio[ABGenotype.ordinal()] + "\t" + posterior_befor[ABGenotype.ordinal()]);
            	System.err.println("BBGenotype before normaliz " + likelihoods_befor[BBGenotype.ordinal()] + "\t" + prio[BBGenotype.ordinal()] + "\t" + posterior_befor[BBGenotype.ordinal()]);
            	System.err.println("Cytosine status: C-neg: " + numCNegStrand + "\tC-pos: " + numCPosStrand + "\tT-neg: " + numTNegStrand + "\tT-pos: " + numTPosStrand);
            	
            }

            	GLs.put(sample.getKey(), new BiallelicGenotypeLikelihoods(sample.getKey(),
            			AlleleA,
            			AlleleB,
            			posterior[AAGenotype.ordinal()],
            			posterior[ABGenotype.ordinal()],
            			posterior[BBGenotype.ordinal()],
                        getFilteredDepth(pileup)));

            
        }

        return refAllele;
	}
	 //original one...
	public double[] normalization(double[] logPosterior, double[] logLikilyhood){
		double sum = 0;
		double[] returnLikilyhood = logPosterior.clone();
		for(int i = 0; i < logLikilyhood.length; i++){
			sum += Math.pow(10,logLikilyhood[i]);
		}
		sum = Math.log10(sum);
		for(int j = 0; j < logLikilyhood.length; j++){
			returnLikilyhood[j] = returnLikilyhood[j] - sum;
		}
		return returnLikilyhood;
	}
	 /*
	public double[] normalization(double[] logPosterior){
		double sum = 0;
		double[] returnLikilyhood = logPosterior.clone();
		for(int i = 0; i < logPosterior.length; i++){
			sum += Math.pow(10,logPosterior[i]);
		}
		sum = Math.log10(sum);
		for(int j = 0; j < logPosterior.length; j++){
			returnLikilyhood[j] = returnLikilyhood[j] - sum;
		}
		return returnLikilyhood;
	}
	*/

	
	protected void initializeBestAndAlternateAlleleFromPosterior(double[] posterior, int location){
		double maxCount = Double.NEGATIVE_INFINITY;
        double secondMaxCount = Double.NEGATIVE_INFINITY;
        DiploidGenotype bestGenotype = DiploidGenotype.createHomGenotype(BaseUtils.A);
        DiploidGenotype secondGenotype = DiploidGenotype.createHomGenotype(BaseUtils.A);
        bestAllele = null;
        alternateAllele = null;
        
        for ( DiploidGenotype g : DiploidGenotype.values() ){
			if(posterior[g.ordinal()] > maxCount){
				secondMaxCount = maxCount;
				maxCount = posterior[g.ordinal()];
				if(bestGenotype.base1 != secondGenotype.base1){
            		secondGenotype = bestGenotype;
            	}
				bestGenotype = g;
			}
			else if (posterior[g.ordinal()] > secondMaxCount && posterior[g.ordinal()] <= maxCount){
	            	secondMaxCount = posterior[g.ordinal()];
	            	secondGenotype = g;
	        }
		}
        if(bestGenotype.isHom()){
        	bestAllele = bestGenotype.base1;
        	if(secondGenotype.isHom()){
        		alternateAllele = secondGenotype.base1;
        	}	
        	else{
        		if(secondGenotype.base1 == bestAllele){
        			alternateAllele = secondGenotype.base2;
        		}
        		else{
        			alternateAllele = secondGenotype.base1;
        		}
        	}
        		
        }
        else{
        	DiploidGenotype temp1 = DiploidGenotype.createHomGenotype(bestGenotype.base1);
        	DiploidGenotype temp2 = DiploidGenotype.createHomGenotype(bestGenotype.base2);
        	if(posterior[temp1.ordinal()] > posterior[temp2.ordinal()]){
        		bestAllele = bestGenotype.base1;
        		alternateAllele = bestGenotype.base2;
        	}
        	else{
        		bestAllele = bestGenotype.base2;
        		alternateAllele = bestGenotype.base1;
        	}
        }
        
		
        if(location == testLoc){
        	for ( DiploidGenotype g : DiploidGenotype.values() ){
        		System.err.println(g.base1 + "-" + g.base2 + ": " + posterior[g.ordinal()]);
        	}
        	System.err.println("bestAllele: " + bestAllele + "\t" + maxCount);
        	if(alternateAllele != null){
        		System.err.println("AlternateAllele: " + "\t" + alternateAllele + "\t" + secondMaxCount);
        	}
        }
	}
	/*
	protected void getBestGenotypeFromPosterior(double[] posterior,HashMap<Integer,methyStatus> cytosineAdjacent, int key, int location){
		double maxCount = Double.NEGATIVE_INFINITY;
        double secondMaxCount = Double.NEGATIVE_INFINITY;
        methyStatus tmpMethyStatus = new methyStatus();
        tmpMethyStatus.genotype = null;
        tmpMethyStatus.ratio = 0.0;
        DiploidGenotype bestGenotype = DiploidGenotype.createHomGenotype(BaseUtils.A);
 
        for ( DiploidGenotype g : DiploidGenotype.values() ){
			if(posterior[g.ordinal()] > maxCount){
				secondMaxCount = maxCount;
				maxCount = posterior[g.ordinal()];
				
				bestGenotype = g;
			}
			else if (posterior[g.ordinal()] > secondMaxCount && posterior[g.ordinal()] <= maxCount){
	            	secondMaxCount = posterior[g.ordinal()];
	        }
		}
        tmpMethyStatus.ratio = 10 * (maxCount - secondMaxCount);
        if(location == BAC.testLocus){
        	//System.err.println("maxCount: " + maxCount + "\tsecondMaxCount: " + secondMaxCount + "\tratio: " + methyStatus.ratio + "\tgenotype: " + bestGenotype);
        	//for(double poster : posterior){
        	//	System.err.println(poster);
        	//}
        }
        if(autoEstimateC && !secondIteration){
        	if( tmpMethyStatus.ratio < this.BAC.cTypeThreshold + this.BAC.STANDARD_CONFIDENCE_FOR_CALLING ){
        		bestGenotype = null;
        	}
        }
        else{
        	if(tmpMethyStatus.ratio < this.BAC.STANDARD_CONFIDENCE_FOR_CALLING){
        		bestGenotype = null;
        	}
        }
        tmpMethyStatus.genotype = bestGenotype;
        cytosineAdjacent.put(key, tmpMethyStatus);
        
      
	}
	*/
	
	protected void getBestGenotypeFromPosterior(double[] posterior,HashMap<Integer,methyStatus> cytosineAdjacent, int key, int location){
		double maxCount = Double.NEGATIVE_INFINITY;
        double secondMaxCount = Double.NEGATIVE_INFINITY;
        methyStatus tmpMethyStatus = new methyStatus();
        tmpMethyStatus.genotype = null;
        tmpMethyStatus.ratio = 0.0;
        DiploidGenotype bestGenotype = DiploidGenotype.createHomGenotype(BaseUtils.A);
 
        for ( DiploidGenotype g : DiploidGenotype.values() ){
			if(posterior[g.ordinal()] > maxCount){
				secondMaxCount = maxCount;
				maxCount = posterior[g.ordinal()];
				
				bestGenotype = g;
			}
			else if (posterior[g.ordinal()] > secondMaxCount && posterior[g.ordinal()] <= maxCount){
	            	secondMaxCount = posterior[g.ordinal()];
	        }
		}
        tmpMethyStatus.ratio = 10 * (maxCount - secondMaxCount);
        if(location == BAC.testLocus){
        	System.err.println("maxCount: " + maxCount + "\tsecondMaxCount: " + secondMaxCount + "\tratio: " + tmpMethyStatus.ratio + "\tgenotype: " + bestGenotype);
        	for(double poster : posterior){
        		System.err.println(poster);
        	}
        }
       // if(autoEstimateC && !secondIteration){
       // 	if( tmpMethyStatus.ratio < this.BAC.cTypeThreshold + this.BAC.STANDARD_CONFIDENCE_FOR_CALLING ){
       // 		bestGenotype = null;
       // 	}
       // }
      //  else{
      //  	if(tmpMethyStatus.ratio < this.BAC.STANDARD_CONFIDENCE_FOR_CALLING){
      //  		bestGenotype = null;
      //  	}
      //  }
        tmpMethyStatus.genotype = bestGenotype;
        cytosineAdjacent.put(key, tmpMethyStatus);
        
      
	}

	protected Integer[] getCytosineStatus(){
		Integer[] value = new Integer[4];
		value[0] = numCNegStrand;
		value[1] = numCPosStrand;
		value[2] = numTNegStrand;
		value[3] = numTPosStrand;
		return value;
	}
	
	protected CytosineTypeStatus getCytosineTypeStatus(){	
		return this.cts;
	}
	
	public BisulfiteDiploidSNPGenotypeLikelihoods checkCytosineStatus(ReadBackedPileup pileup, CytosineTypeStatus cts, double threshold, RefMetaDataTracker tracker,ReferenceContext ref, BisulfiteDiploidSNPGenotypePriors priors){
		double[] cytosineMethyStatus = new double[2]; // 0: methy status in positive strand; 1: methy status in negative strand;
		BisulfiteDiploidSNPGenotypeLikelihoods maxGL = null;
		GenomeLoc location = pileup.getLocation();
		String contig = location.getContig();
		int position = location.getStart();
		double maxRatio = Double.NEGATIVE_INFINITY;
		double[] tmpMethy = new double[2];
		tmpMethy[0] = FALT_METHY_STATUS;
		tmpMethy[1] = FALT_METHY_STATUS;
		HashMap<Integer,methyStatus> cytosineAndAdjacent = new HashMap<Integer,methyStatus>();
		//check adjacent position likelihood
		for(int i = 0 - cts.maxCytosineLength; i <= cts.maxCytosineLength; i++){
			GenomeLoc loc = ref.getGenomeLocParser().createGenomeLoc(contig, position + i );
			if(i == 0)
				continue;
			List<SAMRecord> reads =  new ArrayList<SAMRecord>();;
			List<Integer> elementOffsets = new ArrayList<Integer>();

			for ( PileupElement p : pileup ) {
					int elementOffset = i + p.getOffset();
					if(elementOffset < 0 || elementOffset > p.getRead().getReadLength()-1)
						continue;
					elementOffsets.add(elementOffset);
					reads.add(p.getRead());
			}
			ReadBackedPileup tmpPileup = new ReadBackedPileupImpl(loc,reads,elementOffsets);
			
			ReferenceContext tmpRef = new ReferenceContext(ref.getGenomeLocParser(),loc, ref.getWindow(),ref.getBases());
			 
			
			BisulfiteDiploidSNPGenotypeLikelihoods tmpGL = new BisulfiteDiploidSNPGenotypeLikelihoods(tracker, tmpRef, (BisulfiteDiploidSNPGenotypePriors)priors, BAC, tmpMethy.clone());
			
				
			tmpGL.setPriors(tracker, tmpRef, BAC.PCR_error, BAC.novelDbsnpHet, BAC.validateDbsnpHet, loc);
			if(position == BAC.testLocus){
            	System.err.println("i: " + i + "\ttmpRef: " + tmpRef.getBase());
            	tmpGL.VERBOSE = true;
            }
			int nGoodBases = tmpGL.add(tmpPileup, true, true);
            if ( nGoodBases == 0 )
                continue;
            double[] posteriorNormalized = normalization(tmpGL.getPosteriors(),tmpGL.getLikelihoods());
            
            getBestGenotypeFromPosterior(posteriorNormalized, cytosineAndAdjacent, i ,position);
            
		}
		BisulfiteDiploidSNPGenotypeLikelihoods tmpGL = new BisulfiteDiploidSNPGenotypeLikelihoods(tracker, ref, (BisulfiteDiploidSNPGenotypePriors)priors, BAC, tmpMethy.clone());
		tmpGL.setPriors(tracker, ref, BAC.PCR_error, BAC.novelDbsnpHet, BAC.validateDbsnpHet, location);
		boolean firstSeen = true;
		for(String cytosineType : cts.cytosineListMap.keySet()){
			String[] tmpKey = cytosineType.split("-");
			Double[] value = cts.cytosineListMap.get(cytosineType);
			tmpMethy[0] = value[2];
			tmpMethy[1] = value[2];
			int cytosinePos = Integer.parseInt(tmpKey[1]);
			
			
            double adjacentCytosineSeqLikelihood = 0;
			double adjacentCytosineSeqLikelihoodReverseStrand = 0;
			int i = 1;
			int countMatchedOnFwd = 0;
			int countMatchedOnRvd = 0;
            for(byte base : tmpKey[0].getBytes()){
            	int pos = i - cytosinePos;
            	i++;
            	if(pos == 0)
            		continue;
            	methyStatus tmpMethyStatus = cytosineAndAdjacent.get(pos);
            	if(tmpMethyStatus == null){
            		break;
            	}
            	else if(tmpMethyStatus.genotype == null){
	            	break;
	            }
	            else {
	            	 if(autoEstimateC && !secondIteration){
	                 	if( tmpMethyStatus.ratio < this.BAC.cTypeThreshold + this.BAC.STANDARD_CONFIDENCE_FOR_CALLING ){
	                 		break;
	                 	}
	                 }
	                 else{
	                 	if(tmpMethyStatus.ratio < this.BAC.STANDARD_CONFIDENCE_FOR_CALLING){
	                 		break;
	                 	}
	                 }
	            	 
	            	 if(tmpMethyStatus.genotype.isHet()){
	 	            	break;
	 	             }	
	 	             else{
	 	            	 
	 	            	if(BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(base, tmpMethyStatus.genotype.base1)){
	 	            		countMatchedOnFwd++;
	 	            		adjacentCytosineSeqLikelihood += tmpMethyStatus.ratio;
	 	            	}
	 	            	
	 	             }
	            }
            	
	            if(position == BAC.testLocus){
	            	System.err.println("base: " + (char)base + "\tgenotype: " + (char)tmpMethyStatus.genotype.base1 + "\tcytosinePos: " + cytosinePos + "\tratio: " + tmpMethyStatus.ratio + "\tadjacentCytosineSeqLikelihood: " + adjacentCytosineSeqLikelihood);
	            }
	            
            }
            i = 1;
            for(byte base : tmpKey[0].getBytes()){
            	int pos = cytosinePos - i;
            	i++;
            	if(pos == 0)
            		continue;
            	methyStatus tmpMethyStatus = cytosineAndAdjacent.get(pos);
            	//DiploidGenotype g = tmpMethyStatus.genotype;
	           // double ratio = tmpMethyStatus.ratio;
            	if(tmpMethyStatus == null){
            		break;
            	}
            	else if(tmpMethyStatus.genotype == null){
	            	break;
	            }
            	else {
	            	 if(autoEstimateC && !secondIteration){
	                 	if( tmpMethyStatus.ratio < this.BAC.cTypeThreshold + this.BAC.STANDARD_CONFIDENCE_FOR_CALLING ){
	                 		break;
	                 	}
	                 }
	                 else{
	                 	if(tmpMethyStatus.ratio < this.BAC.STANDARD_CONFIDENCE_FOR_CALLING){
	                 		break;
	                 	}
	                 }
	            	 
	            	 if(tmpMethyStatus.genotype.isHet()){
	 	            	break;
	 	             }	
	 	             else{
	 	            	 
	 	            	if(BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(base, BaseUtilsMore.iupacCodeComplement(tmpMethyStatus.genotype.base1))){
		            		countMatchedOnRvd++;
		            		adjacentCytosineSeqLikelihoodReverseStrand += tmpMethyStatus.ratio;
		            	}
	 	            	
	 	             }
	            }
	            
	            if(position == BAC.testLocus){
	            	System.err.println("base: " + (char)base + "\tgenotype: " + (char)BaseUtilsMore.iupacCodeComplement(tmpMethyStatus.genotype.base1) + "\treveser: " + (char)base + "\tcytosinePos: " + cytosinePos + "\tratio: " + tmpMethyStatus.ratio + "\tadjacentCytosineSeqLikelihoodReverseStrand: " + adjacentCytosineSeqLikelihoodReverseStrand);
	            }
	            
            }
            if((countMatchedOnFwd < tmpKey[0].length() - 1) && (countMatchedOnRvd < tmpKey[0].length() - 1))
            	continue;
            
            if(autoEstimateC && !secondIteration && !firstSeen){
            	
            }
            else{
            	firstSeen = false;
            	tmpGL.clearLikelihoods(tmpMethy.clone());
            	if(position == BAC.testLocus){
            		tmpGL.VERBOSE = true;
            		System.err.println("cytosineType: " + cytosineType );
            	}
            	
     			int nGoodBases = tmpGL.add(pileup, true, true);
                 if ( nGoodBases == 0 )
                     break;
                 double[] posteriorNormalized = normalization(tmpGL.getPosteriors(),tmpGL.getLikelihoods());
                 
                 getBestGenotypeFromPosterior(posteriorNormalized, cytosineAndAdjacent, 0 ,position);
            }
            
           
            
            methyStatus tmpMethyStatus = cytosineAndAdjacent.get(0);
        	if(tmpMethyStatus == null){
        		continue;
        	}
        	else if(tmpMethyStatus.genotype == null){
        		continue;
            }
            else if(tmpMethyStatus.genotype.isHet()){
            	if(tmpKey[0].length() == 1){
            		try {
						maxGL = (BisulfiteDiploidSNPGenotypeLikelihoods) tmpGL.clone();
					} catch (CloneNotSupportedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
            	}
            	else{
            		continue;
            	}
            		
            }	
            else{
            	
            	if(BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(BaseUtils.C, tmpMethyStatus.genotype.base1)){
            		if(autoEstimateC && !secondIteration){
                     	if( tmpMethyStatus.ratio < this.BAC.cTypeThreshold + this.BAC.STANDARD_CONFIDENCE_FOR_CALLING ){
                     		continue;
                     	}
                     }
                     else{
                     	if(tmpMethyStatus.ratio < this.BAC.STANDARD_CONFIDENCE_FOR_CALLING){
                     		continue;
                     	}
                     }
            		countMatchedOnFwd++;
            		adjacentCytosineSeqLikelihood += tmpMethyStatus.ratio;
            	}
            	else if(BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(BaseUtils.C, BaseUtilsMore.iupacCodeComplement(tmpMethyStatus.genotype.base1))){
            		if(autoEstimateC && !secondIteration){
                     	if( tmpMethyStatus.ratio < this.BAC.cTypeThreshold + this.BAC.STANDARD_CONFIDENCE_FOR_CALLING ){
                     		continue;
                     	}
                     }
                     else{
                     	if(tmpMethyStatus.ratio < this.BAC.STANDARD_CONFIDENCE_FOR_CALLING){
                     		continue;
                     	}
                     }
            		countMatchedOnRvd++;
            		adjacentCytosineSeqLikelihoodReverseStrand += tmpMethyStatus.ratio;
            	}
            	else{
            		if(tmpKey[0].length() == 1){
                		try {
    						maxGL = (BisulfiteDiploidSNPGenotypeLikelihoods) tmpGL.clone();
    					} catch (CloneNotSupportedException e) {
    						// TODO Auto-generated catch block
    						e.printStackTrace();
    					}
                	}
                	
            	}
            	
            }
            
            if(countMatchedOnFwd >= tmpKey[0].length()){
				value[3] = 1.0;
				value[0] = adjacentCytosineSeqLikelihood;
				if(adjacentCytosineSeqLikelihood > maxRatio){
					maxRatio = adjacentCytosineSeqLikelihood;
					cytosineMethyStatus[0] = value[2];
					try {
						maxGL = (BisulfiteDiploidSNPGenotypeLikelihoods) tmpGL.clone();
					} catch (CloneNotSupportedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					
				}

			}
			else if(countMatchedOnRvd >= tmpKey[0].length()){
				value[3] = 1.0;
				value[1] = adjacentCytosineSeqLikelihoodReverseStrand;
				if(adjacentCytosineSeqLikelihoodReverseStrand > maxRatio){
					maxRatio = adjacentCytosineSeqLikelihoodReverseStrand;
					cytosineMethyStatus[1] = value[2];
					try {
						maxGL = (BisulfiteDiploidSNPGenotypeLikelihoods) tmpGL.clone();
					} catch (CloneNotSupportedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				
					
			}
			if(position == BAC.testLocus){
            	System.err.println("countMatchedOnFwd: " + countMatchedOnFwd + "\tcountMatchedOnRvd: " + countMatchedOnRvd);
            } 
		}
		
		
		return maxGL;	
	}
	
/*
	public double[] checkCytosineStatus(ReadBackedPileup pileup, CytosineTypeStatus cts, double threshold, RefMetaDataTracker tracker,ReferenceContext ref, BisulfiteDiploidSNPGenotypePriors priors){
		double[] cytosineMethyStatus = new double[2]; // 0: methy status in positive strand; 1: methy status in negative strand;
		GenomeLoc location = pileup.getLocation();
		String contig = location.getContig();
		int position = location.getStart();
		double maxRatio = Double.NEGATIVE_INFINITY;
		for(String cytosineType : cts.cytosineListMap.keySet()){
			String[] tmpKey = cytosineType.split("-");
			Double[] value = cts.cytosineListMap.get(cytosineType);
			double[] tmpMethy = new double[2];
			tmpMethy[0] = value[2];
			tmpMethy[1] = value[2];
			double adjacentCytosineSeqLikelihood = 0;
			double adjacentCytosineSeqLikelihoodReverseStrand = 0;
			int i = 1;
			int countMatchedOnFwd = 0;
			int countMatchedOnRvd = 0;
			for(byte base : tmpKey[0].getBytes()){
				int cytosinePos = Integer.parseInt(tmpKey[1]);
				GenomeLoc loc = ref.getGenomeLocParser().createGenomeLoc(contig, position + i - cytosinePos);	
				if( i != cytosinePos){
					tmpMethy[0] = FALT_METHY_STATUS;
					tmpMethy[1] = FALT_METHY_STATUS;
				}
				List<SAMRecord> reads =  new ArrayList<SAMRecord>();;
				List<Integer> elementOffsets = new ArrayList<Integer>();

				for ( PileupElement p : pileup ) {
						int elementOffset = i - cytosinePos + p.getOffset();
						if(elementOffset < 0 || elementOffset > p.getRead().getReadLength()-1)
							continue;
						elementOffsets.add(elementOffset);
						reads.add(p.getRead());
				}
				
				ReadBackedPileup tmpPileup = new ReadBackedPileupImpl(loc,reads,elementOffsets);
				BisulfiteDiploidSNPGenotypeLikelihoods tmpGL = new BisulfiteDiploidSNPGenotypeLikelihoods(tracker, ref, (BisulfiteDiploidSNPGenotypePriors)priors, BAC, tmpMethy.clone());
				
					
				tmpGL.setPriors(tracker, ref, BAC.PCR_error, BAC.novelDbsnpHet, BAC.validateDbsnpHet);
				int nGoodBases = tmpGL.add(tmpPileup, true, true);
	            if ( nGoodBases == 0 )
	                continue;
	            double[] posteriorNormalized = normalization(tmpGL.getPosteriors().clone(),tmpGL.getLikelihoods().clone());
	            
	            getBestGenotypeFromPosterior(posteriorNormalized,position);
	            DiploidGenotype g = methyStatus.genotype;
	            double ratio = methyStatus.ratio;
	            if(g == null){
	            	break;
	            }
	            else if(g.isHet()){
	            	break;
	            }	
	            else{
	            	if(BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(base, g.base1)){
	            		countMatchedOnFwd++;
	            		adjacentCytosineSeqLikelihood += ratio;
	            	}
	            	
	            }
	            if(position == BAC.testLocus){
	            	System.err.println("base: " + (char)base + "\tgenotype: " + (char)g.base1 + "\tcytosinePos: " + cytosinePos + "\tratio: " + ratio + "\tadjacentCytosineSeqLikelihood: " + adjacentCytosineSeqLikelihood);
	            }
	            i++;
			}
			i = 1;
			for(byte base : tmpKey[0].getBytes()){
				int cytosinePos = Integer.parseInt(tmpKey[1]);
				GenomeLoc loc = ref.getGenomeLocParser().createGenomeLoc(contig, position - (i - cytosinePos));			
				List<SAMRecord> reads =  new ArrayList<SAMRecord>();;
				List<Integer> elementOffsets = new ArrayList<Integer>();

				for ( PileupElement p : pileup ) {
						int elementOffset = p.getOffset() - (i - cytosinePos);
						if(elementOffset < 0 || elementOffset > p.getRead().getReadLength()-1)
							continue;
						elementOffsets.add(elementOffset);
						reads.add(p.getRead());
				}
				
				ReadBackedPileup tmpPileup = new ReadBackedPileupImpl(loc,reads,elementOffsets);
				BisulfiteDiploidSNPGenotypeLikelihoods tmpGL = new BisulfiteDiploidSNPGenotypeLikelihoods(tracker, ref, (BisulfiteDiploidSNPGenotypePriors)priors, BAC, tmpMethy.clone());
			
				tmpGL.setPriors(tracker, ref, BAC.PCR_error, BAC.novelDbsnpHet, BAC.validateDbsnpHet);
				int nGoodBases = tmpGL.add(tmpPileup, true, true);
	            if ( nGoodBases == 0 )
	                continue;
	            double[] posteriorNormalized = normalization(tmpGL.getPosteriors().clone(),tmpGL.getLikelihoods().clone());

	            getBestGenotypeFromPosterior(posteriorNormalized,position);
	            DiploidGenotype g = methyStatus.genotype;
	            double ratio = methyStatus.ratio;
	            if(g == null){
	            	break;
	            }
	            else if(g.isHet()){
	            	break;
	            }	
	            else{

	            	if(BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(base, BaseUtilsMore.iupacCodeComplement(g.base1))){
	            		countMatchedOnRvd++;
	            		adjacentCytosineSeqLikelihoodReverseStrand += ratio;
	            	}
	            }
	            if(position == BAC.testLocus){
	            	System.err.println("base: " + (char)base + "\tgenotype: " + (char)BaseUtilsMore.iupacCodeComplement(g.base1) + "\treveser: " + (char)base + "\tcytosinePos: " + cytosinePos + "\tratio: " + ratio + "\tadjacentCytosineSeqLikelihoodReverseStrand: " + adjacentCytosineSeqLikelihoodReverseStrand);
	            }
	            i++;
			}
			if(countMatchedOnFwd >= tmpKey[0].length()){
				value[3] = 1.0;
				value[0] = adjacentCytosineSeqLikelihood;
				if(adjacentCytosineSeqLikelihood > maxRatio){
					maxRatio = adjacentCytosineSeqLikelihood;
					cytosineMethyStatus[0] = value[2];
					
				}

			}
			else if(countMatchedOnRvd >= tmpKey[0].length()){
				value[3] = 1.0;
				value[1] = adjacentCytosineSeqLikelihoodReverseStrand;
				if(adjacentCytosineSeqLikelihoodReverseStrand > maxRatio){
					maxRatio = adjacentCytosineSeqLikelihoodReverseStrand;
					cytosineMethyStatus[1] = value[2];
				}
					
			}
			if(position == BAC.testLocus){
            	System.err.println("countMatchedOnFwd: " + countMatchedOnFwd + "\tcountMatchedOnRvd: " + countMatchedOnRvd);
            } 
			

		}
		return cytosineMethyStatus;	
	}

	

	

	static class methyStatus{
		static DiploidGenotype genotype;
		static double ratio;
	}
	
	protected void getBestGenotypeFromPosterior(double[] posterior, int location){
		double maxCount = Double.NEGATIVE_INFINITY;
        double secondMaxCount = Double.NEGATIVE_INFINITY;
        methyStatus.genotype = null;
        methyStatus.ratio = 0.0;
        DiploidGenotype bestGenotype = DiploidGenotype.createHomGenotype(BaseUtils.A);
    //    DiploidGenotype secondGenotype = DiploidGenotype.createHomGenotype(BaseUtils.A);

        
        for ( DiploidGenotype g : DiploidGenotype.values() ){
			if(posterior[g.ordinal()] > maxCount){
				secondMaxCount = maxCount;
				maxCount = posterior[g.ordinal()];
				//if(bestGenotype.base1 != secondGenotype.base1){
      //      		secondGenotype = bestGenotype;
            	//}
				bestGenotype = g;
			}
			else if (posterior[g.ordinal()] > secondMaxCount && posterior[g.ordinal()] <= maxCount){
	            	secondMaxCount = posterior[g.ordinal()];
	     //       	secondGenotype = g;
	        }
		}
        methyStatus.ratio = 10 * (maxCount - secondMaxCount);
        if(location == BAC.testLocus){
        	//System.err.println("maxCount: " + maxCount + "\tsecondMaxCount: " + secondMaxCount + "\tratio: " + methyStatus.ratio + "\tgenotype: " + bestGenotype);
        	//for(double poster : posterior){
        	//	System.err.println(poster);
        	//}
        }
        if(autoEstimateC && !secondIteration){
        	if( methyStatus.ratio < this.BAC.cTypeThreshold){
        		bestGenotype = null;
        	}
        }
        else{
        	if(methyStatus.ratio < this.BAC.STANDARD_CONFIDENCE_FOR_CALLING){
        		bestGenotype = null;
        	}
        }
        methyStatus.genotype = bestGenotype;
        
        //return methyStatus;
	}
	*/
	
	public static class methyStatus{
		DiploidGenotype genotype;
		double ratio;
		methyStatus(){
			
		}
	}

}
