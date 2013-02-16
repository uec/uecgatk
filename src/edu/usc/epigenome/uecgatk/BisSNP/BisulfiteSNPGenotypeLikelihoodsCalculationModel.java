package edu.usc.epigenome.uecgatk.BisSNP;


import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import net.sf.samtools.SAMRecord;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.variantcontext.Allele;

import edu.usc.epigenome.uecgatk.BisSNP.BisulfiteDiploidSNPGenotypeLikelihoods;
import edu.usc.epigenome.uecgatk.BisSNP.BisulfiteDiploidSNPGenotypePriors;
import edu.usc.epigenome.uecgatk.BisSNP.BisulfiteEnums.INVERT_DUPS;
import edu.usc.epigenome.uecgatk.BisSNP.BisulfiteEnums.MethylSNPModel;
import edu.usc.epigenome.uecgatk.bisulfiteIndels.BisBAQ;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.BadMateFilter;

import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;

import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;

import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.MathUtils;

import org.broadinstitute.sting.gatk.walkers.genotyper.DiploidGenotype;
import org.broadinstitute.sting.gatk.walkers.genotyper.SNPGenotypeLikelihoodsCalculationModel.BAQedPileupElement;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

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

public class BisulfiteSNPGenotypeLikelihoodsCalculationModel{

	private BisulfiteArgumentCollection BAC;
	
	protected Byte bestAllele = null;
	protected Byte alternateAllele = null;
	protected long testLoc;
	protected int numCNegStrand = 0;
	protected int numTNegStrand = 0;
	protected int numOtherNegStrand = 0;
	protected int numCPosStrand = 0;
	protected int numTPosStrand = 0;
	protected int numOtherPosStrand = 0;
	private HashMap<String, BisulfiteContextsGenotypeLikelihoods> BCGLs = null;
	
	private boolean autoEstimateC = false;
    private boolean secondIteration = false;
    private double FLAT_METHY_STATUS = 0.5;
    private String CYTOSINE_CONV_START = "Xn";
    private boolean useBAQ = false;
    
    private Allele refAllele = null;
	

	
	
	public BisulfiteSNPGenotypeLikelihoodsCalculationModel(BisulfiteArgumentCollection BAC, boolean useBAQ) {
	
		// TODO Auto-generated constructor stub
		this.BAC = BAC;

		//this.cts = cts;
		this.testLoc = BAC.testLocus;
		this.autoEstimateC = autoEstimateC;
		this.secondIteration = secondIteration;
		this.useBAQ = useBAQ;
		if(BAC.sequencingMode == MethylSNPModel.NM){
			FLAT_METHY_STATUS = 0.0;
		}
	//	BCGLs = new HashMap<String, BisulfiteContextsGenotypeLikelihoods>();
	}
	

	public HashMap<String, BisulfiteContextsGenotypeLikelihoods> getBsContextGenotypeLikelihoods(){
		
		return BCGLs;
	}
	
	public Allele getRefAllele(){
		return refAllele;
	}
	

	public void setBsLikelihoods(RefMetaDataTracker tracker,
			ReferenceContext ref,
			Map<String, AlignmentContext> contexts,
			AlignmentContextUtils.ReadOrientation contextType,
			HashMap<String, BisulfiteContextsGenotypeLikelihoods> BCGLs){

		byte refBase = ref.getBase();
		
		if ( !Allele.acceptableAlleleBases(new byte[]{refBase}) )  //illegal reference sequence
			refBase = 'N';
		
		this.refAllele = Allele.create(refBase, true);

			
    //    numCNegStrand = 0;
    //    numTNegStrand = 0;
    //    numCPosStrand = 0;
    //    numTPosStrand = 0;
        

        for ( Map.Entry<String, AlignmentContext> sample : contexts.entrySet() ) {
        	// when pairs of reads overlapped, get rid of both of them when disagree, get rid of bad base qual reads when they agree
            ReadBackedPileup pileup =AlignmentContextUtils.stratify(sample.getValue(),contextType).getBasePileup().getOverlappingFragmentFilteredPileup();
            if(useBAQ){
            	pileup = createBAQedPileup( pileup);
            }
            
            int numCNegStrand = 0;
        	int numTNegStrand = 0;
        	int numANegStrand = 0;
        	int numGNegStrand = 0;
        	int numCPosStrand = 0;
        	int numTPosStrand = 0;
        	int numAPosStrand = 0;
        	int numGPosStrand = 0;
        	
        	int numGGalleleStrand = 0;
        	int numAGalleleStrand = 0;
        	int numOtherGalleleStrand = 0;
        	int numCCalleleStrand = 0;
        	int numTCalleleStrand = 0;
        	int numOtherCalleleStrand = 0;
        	
           // System.err.println(sample + "\t" + cts.toString());
            for ( PileupElement p : pileup ) {
            	SAMRecord samRecord = p.getRead();
            	boolean negStrand = samRecord.getReadNegativeStrandFlag();
            	if(samRecord.getDuplicateReadFlag()){ //get rid of duplicate reads
                	continue;
                }
            	int offset = p.getOffset();
            	if(offset < 0)//is deletion
            		continue;
            	
            	boolean paired = samRecord.getReadPairedFlag();
            	if(paired){
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
    					continue;
    				}
    				
    				// Inverted dups, count only one end
    				if (samRecord.getAlignmentStart() == samRecord.getMateAlignmentStart() && samRecord.getReadNegativeStrandFlag() == samRecord.getMateNegativeStrandFlag())
    				{
    					if(BAC.invDups == INVERT_DUPS.USE_ONLY_1ST_END){
    						if (samRecord.getSecondOfPairFlag()) continue;
    					}
    					else if(BAC.invDups == INVERT_DUPS.NOT_TO_USE){
    						continue;
    					}
    					else{
    						
    					}
    					
       				}
    	        	if (Paired  && !BAC.USE_BADLY_MATED_READS && !samRecord.getProperPairFlag())
    				{
    					continue;
    				}
    	        	
    	        	if((pileup.getLocation().getStart()) == testLoc){
            		//	System.err.println("NegativeStrandFlag: " + samRecord.getReadNegativeStrandFlag() + "\t" + "MateNegativeStrandFlag: " + samRecord.getMateNegativeStrandFlag() + "\tbase: " + samRecord.getReadBases()[offset] + "\t" + "baseQ: " + samRecord.getBaseQualities()[offset]);
            		//	System.err.println("getReadString: " + samRecord.getReadString() + "\tsecond: " + secondOfPair);
            		}
    	        	if(secondOfPair){	        		
    		        	//samRecord.setReadNegativeStrandFlag(!samRecord.getReadNegativeStrandFlag());    
    	        		//observedBase = BaseUtils.simpleComplement(observedBase);;
    	        		negStrand = !negStrand;
    	        	}
    	        	if((pileup.getLocation().getStart()) == testLoc){
            		//	System.err.println("proper paired: " + samRecord.getProperPairFlag() + "\t" + "getMateAlignmentStart: " + samRecord.getMateAlignmentStart() + "\t" + "MateNegativeStrandFlag: " + samRecord.getMateNegativeStrandFlag());
    	        	//	System.err.println("getReadString: " + samRecord.getReadString() + "\t" + "getAlignmentStart: " + samRecord.getAlignmentStart() + "\t" + "getUnclippedEnd: " + samRecord.getUnclippedEnd() + "\t" +"NegativeStrandFlag: " + samRecord.getReadNegativeStrandFlag() + "\tcytosineOffset: " + offset + "\tbase: " + samRecord.getReadBases()[offset] + "\t" + "baseQ: " + samRecord.getBaseQualities()[offset]);
    	        	//	System.err.println("getBase: " + p.getBase() + "\tp.getRead(): " + p.getRead().getReadBases()[offset]);
            		}
            	}
				
	        	
				int alignmentS = samRecord.getAlignmentStart();
				int	onRefCoord = (negStrand) ? samRecord.getUnclippedEnd() : alignmentS;
				
				//summary number of C,T in the positive and negative strand
				//BadBaseFilterBisulfite badReadPileupFilter = new BadBaseFilterBisulfite(ref, BAC);
				GATKSAMRecordFilterStorage GATKrecordFilterStor = new GATKSAMRecordFilterStorage((GATKSAMRecord)p.getRead(), BAC,ref, p.getOffset());
				setGoodConvertedBase(p.getRead(), ref, negStrand);
               //GATKrecordFilterStor.setGoodBases(badReadPileupFilter, true);
				// && isGoodConvertedBase(p.getRead(), offset)
				if(GATKrecordFilterStor.isGoodBase()){
					if(negStrand){
						if(p.getBase()==BaseUtils.G){
							numCNegStrand++;
						}
						else if(p.getBase()==BaseUtils.A){
							numTNegStrand++;
						}
						else if(p.getBase()==BaseUtils.C){
							numGNegStrand++;
						}
						else if(p.getBase()==BaseUtils.T){
							numANegStrand++;
						}
						
					}
					else{
						if(p.getBase()==BaseUtils.C){
							numCPosStrand++;
						}
						else if(p.getBase()==BaseUtils.T){
							numTPosStrand++;
						}
						else if(p.getBase()==BaseUtils.G){
							numGPosStrand++;
						}
						else if(p.getBase()==BaseUtils.A){
							numAPosStrand++;
						}
					}
				}
				
				
				if((pileup.getLocation().getStart()) == testLoc){
					System.err.println("before filter:\t" + onRefCoord + "\t" + offset + "\t" + negStrand + "\t" + pileup.getLocation().getStart() + "\t" + (char)p.getBase() + "\t" + p.getQual());
					System.err.println("refBase: " + refBase + "\tGoodBase: " + GATKrecordFilterStor.isGoodBase());
					//System.err.println(((GATKSAMRecord)p.getRead()).getMappingQuality() + "\t" + ((GATKSAMRecord)p.getRead()).getBaseQualities()[offset] + "\t" + BAC.USE_BADLY_MATED_READS + "\t" + (!BadMateFilter.hasBadMate(((GATKSAMRecord)p.getRead()))) + "\t" + !((GATKSAMRecord)p.getRead()).getNotPrimaryAlignmentFlag() );
					
					if(paired)
						System.err.println("isGoodBase: " + GATKrecordFilterStor.isGoodBase() + "\tsecondOfPair: " + "\t" + samRecord.getSecondOfPairFlag());
					    
				}
            }

			BisulfiteDiploidSNPGenotypePriors priors = new BisulfiteDiploidSNPGenotypePriors();

			HashMap<Integer,double[]> GPsBeforeCytosineTenGenotypes = new HashMap<Integer,double[]>();
			HashMap<Integer,double[]> GPsAfterCytosineTenGenotypes = new HashMap<Integer,double[]>();
			HashMap<String, Double[]> GPsAtCytosineTenGenotypes = new HashMap<String, Double[]>();
			HashMap<String,CytosineParameters> cytosineParametersStatus =  new HashMap<String,CytosineParameters>();
			
            String bestMatchedCytosinePattern = checkCytosineStatus(pileup, BAC, tracker, ref, (BisulfiteDiploidSNPGenotypePriors)priors, GPsBeforeCytosineTenGenotypes, GPsAfterCytosineTenGenotypes, GPsAtCytosineTenGenotypes, cytosineParametersStatus);
            //System.err.println(bestMatchedCytosinePattern);
            char cytosineStrand;
            BisulfiteDiploidSNPGenotypeLikelihoods GL;
           // if(bestMatchedCytosinePattern == null){
           // 	GL = new BisulfiteDiploidSNPGenotypeLikelihoods(tracker, ref, (BisulfiteDiploidSNPGenotypePriors)priors, BAC, 0);
            	
           // }
          //  else{
          //  	GL = new BisulfiteDiploidSNPGenotypeLikelihoods(tracker, ref, (BisulfiteDiploidSNPGenotypePriors)priors, BAC, BAC.cytosineDefined.getContextDefined().get(bestMatchedCytosinePattern).cytosineMethylation);
            	
         //   }
            GL = new BisulfiteDiploidSNPGenotypeLikelihoods(tracker, ref, (BisulfiteDiploidSNPGenotypePriors)priors, BAC, getMethyLevelFromPileup(pileup, ref));
            if(cytosineParametersStatus.containsKey(bestMatchedCytosinePattern)){
            	cytosineStrand = cytosineParametersStatus.get(bestMatchedCytosinePattern).cytosineStrand;
            }
            else{
            	cytosineStrand = '+';
            	if((pileup.getLocation().getStart()) == testLoc && bestMatchedCytosinePattern != null)
            		System.err.println(BAC.cytosineDefined.getContextDefined().get(bestMatchedCytosinePattern).cytosineMethylation + "\t" + bestMatchedCytosinePattern);
            	bestMatchedCytosinePattern = null;
            }
            
            if((pileup.getLocation().getStart()) == testLoc)
            	GL.VERBOSE=true;
            

            GL.setPriors(tracker, ref, BAC.heterozygosity, BAC.novelDbsnpHet, BAC.validateDbsnpHet, ref.getLocus());
           
            
            int nGoodBases = GL.add(pileup, true, true);
            
            if ( nGoodBases == 0 )
                continue;
            if( cytosineStrand == '+'){
            	numGGalleleStrand = numGNegStrand;
            	numAGalleleStrand = numANegStrand;
            	numOtherGalleleStrand = numCNegStrand + numTNegStrand;
            	numCCalleleStrand = numCPosStrand;
            	numTCalleleStrand = numTPosStrand;
            	numOtherCalleleStrand = numGPosStrand + numAPosStrand;
            }
            else{
            	numGGalleleStrand = numGPosStrand;
            	numAGalleleStrand = numAPosStrand;
            	numOtherGalleleStrand = numCPosStrand + numTPosStrand;
            	numCCalleleStrand = numCNegStrand;
            	numTCalleleStrand = numTNegStrand;
            	numOtherCalleleStrand = numGNegStrand + numANegStrand;
            }
            
            double[] prio = GL.getPriors();
            double[] likelihoods = GL.getLikelihoods();
            double[] posterior = MathUtils.normalizeFromLog10(GL.getPosteriors(), true, false);

            initializeBestAndAlternateAlleleFromPosterior(posterior, pileup.getLocation().getStart());
            
            if ( (alternateAllele == null && bestAllele == refBase) || (bestAllele == null) ) {
               
                if ( BAC.OutputMode == BisulfiteEnums.OUTPUT_MODE.EMIT_VARIANTS_ONLY )
                    //return refAllele;
                	return;
            }
            
            Allele AlleleA, AlleleB;

            if(alternateAllele == null || BaseUtils.basesAreEqual(alternateAllele,refBase) || alternateAllele == bestAllele){
            	AlleleA = Allele.create(refBase, true);

            	if(BaseUtils.basesAreEqual(bestAllele,refBase)){
            		for(byte base : BaseUtils.BASES){
            			if(base != refBase){
            				bestAllele = base;
            				break;
            			}
            		}
            	}
            	
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
            
           // if(AlleleA.equals(AlleleB, true) && ((AlleleA.isReference() && AlleleB.isNonReference())||(AlleleB.isReference() && AlleleA.isNonReference()))){
           // 	if(AlleleA.isReference() && AlleleB.isNonReference()){
           // 		AlleleB = AlleleA;
          //  	}
          //  	else{
          //  		AlleleA = AlleleB;
          //  	}
          //  }
            	
            DiploidGenotype AAGenotype = DiploidGenotype.createHomGenotype(bestAllele);
            DiploidGenotype ABGenotype = DiploidGenotype.createDiploidGenotype(bestAllele, alternateAllele);
            DiploidGenotype BBGenotype = DiploidGenotype.createHomGenotype(alternateAllele);
            
            
            if((pileup.getLocation().getStart()) == testLoc){
            	System.err.println("sample: " + sample.getKey());
            	System.err.println("sample location: " + pileup.getPileupString((char)refBase));
            	System.err.println("sample: " + sample.getValue().getLocation().getStart());
            	System.err.println("refBase: " + refBase + " bestAllele: " + bestAllele + " alternateAllele: " + alternateAllele); 
            	System.err.println("AAGenotype: " + AAGenotype.toString() + " ABGenotype: " + ABGenotype.toString() + " BBGenotype: " + BBGenotype.toString()); 
            	System.err.println("AlleleA: " + AlleleA.toString() + " AlleleB: " + AlleleB.toString()); 
            	System.err.println("AAGenotype " + likelihoods[AAGenotype.ordinal()] + "\t" + prio[AAGenotype.ordinal()] + "\t" + posterior[AAGenotype.ordinal()]);
            	System.err.println("ABGenotype " + likelihoods[ABGenotype.ordinal()] + "\t" + prio[ABGenotype.ordinal()] + "\t" + posterior[ABGenotype.ordinal()]);
            	System.err.println("BBGenotype " + likelihoods[BBGenotype.ordinal()] + "\t" + prio[BBGenotype.ordinal()] + "\t" + posterior[BBGenotype.ordinal()]);
            	System.err.println("Cytosine status: C-neg: " + numCNegStrand + "\tC-pos: " + numCPosStrand + "\tT-neg: " + numTNegStrand + "\tT-pos: " + numTPosStrand);
            	
            }
            Set<String> cytosineContexts = BAC.cytosineDefined.getContextDefined().keySet();
           // System.err.println(posterior[AAGenotype.ordinal()] + "\t" + posterior[ABGenotype.ordinal()] + "\t" + posterior[BBGenotype.ordinal()]);
            
           
            	
            BCGLs.put(sample.getKey(), new BisulfiteContextsGenotypeLikelihoods(sample.getKey(),
            			AlleleA,
            			AlleleB,
            			posterior[AAGenotype.ordinal()],
            			posterior[ABGenotype.ordinal()],
            			posterior[BBGenotype.ordinal()],
                        cytosineContexts,numCCalleleStrand,numTCalleleStrand,numOtherCalleleStrand,
                        numGGalleleStrand, numAGalleleStrand, numOtherGalleleStrand, getFilteredDepth(pileup),
                        cytosineParametersStatus, bestMatchedCytosinePattern,GPsBeforeCytosineTenGenotypes,GPsAfterCytosineTenGenotypes,
                        GPsAtCytosineTenGenotypes, pileup, BAC 
                        ));
        }
       // return refAllele;
        this.BCGLs = BCGLs;
        
	}
	
	public String checkCytosineStatus(ReadBackedPileup pileup, BisulfiteArgumentCollection BAC, RefMetaDataTracker tracker,ReferenceContext ref, BisulfiteDiploidSNPGenotypePriors priors, 
			HashMap<Integer,double[]> GPsBeforeCytosineTenGenotypes, HashMap<Integer,double[]> GPsAfterCytosineTenGenotypes, HashMap<String, Double[]> GPsAtCytosineTenGenotypes, 
			HashMap<String,CytosineParameters> cytosineParametersStatus){
		String bestCytosinePattern = null;
		//String bestCytosinePatternRelative = null; //when there is no best cytosine pattern, just return the maximum Liklihood's cytosine pattern which maybe heterozygous
		double maxRatioInCytosinePos = Double.NEGATIVE_INFINITY;
		GenomeLoc location = pileup.getLocation();
		String contig = location.getContig();
		int position = location.getStart();
		double maxRatio = Double.NEGATIVE_INFINITY;
		double tmpMethy = FLAT_METHY_STATUS;
		//tmpMethy[0] = FLAT_METHY_STATUS;// 0: methy status in positive strand; 1: methy status in negative strand;
		//tmpMethy[1] = FLAT_METHY_STATUS;
		HashMap<Integer,methyStatus> cytosineAndAdjacent = new HashMap<Integer,methyStatus>();
		 
		//check adjacent position likelihood
		int maxCytosineLength=BAC.cytosineDefined.getMaxCytosinePatternLen();
		byte[] refBytes = new byte[maxCytosineLength*2+1];
		//byte[] refBytesFwd = new byte[maxCytosineLength];
		//byte[] refBytesRvd = new byte[maxCytosineLength];
		for(int i = 0 - maxCytosineLength, index=0; i <= maxCytosineLength; i++,index++){
			GenomeLoc loc = ref.getGenomeLocParser().createGenomeLoc(contig, position + i );
			if( !ref.getWindow().containsP(loc) )
				continue;
			
			ReferenceContext tmpRef = new ReferenceContext(ref.getGenomeLocParser(),loc, ref.getWindow(),ref.getBases());
			refBytes[index] = tmpRef.getBase();
			
			if(i == 0)
				continue;
			List<GATKSAMRecord> reads =  new ArrayList<GATKSAMRecord>();;
			List<Integer> elementOffsets = new ArrayList<Integer>();

			for ( PileupElement p : pileup ) {
					int elementOffset = i + p.getOffset();
					if(elementOffset < 0 || elementOffset > p.getRead().getReadLength()-1)
						continue;
					elementOffsets.add(elementOffset);
					reads.add(p.getRead());
			}
			ReadBackedPileup tmpPileup = new ReadBackedPileupImpl(loc,reads,elementOffsets);
			
			
			tmpMethy = getMethyLevelFromPileup(tmpPileup, ref);
			BisulfiteDiploidSNPGenotypeLikelihoods tmpGL = new BisulfiteDiploidSNPGenotypeLikelihoods(tracker, tmpRef, (BisulfiteDiploidSNPGenotypePriors)priors, BAC, tmpMethy);
	
			tmpGL.setPriors(tracker, tmpRef, BAC.heterozygosity, BAC.novelDbsnpHet, BAC.validateDbsnpHet, loc);
			if(position == BAC.testLocus){
            	System.err.println("i: " + i + "\ttmpRef: " + tmpRef.getBase());
            	tmpGL.VERBOSE = true;
            }
			int nGoodBases = tmpGL.add(tmpPileup, true, true);
            if ( nGoodBases == 0 )
                continue;
            double[] posteriorNormalized = MathUtils.normalizeFromLog10(tmpGL.getPosteriors(), true, false);
            Integer distanceToCytosine = Math.abs(i);
            if(i<0){
            	GPsBeforeCytosineTenGenotypes.put(distanceToCytosine, posteriorNormalized.clone());
            }
            else{
            	GPsAfterCytosineTenGenotypes.put(distanceToCytosine, posteriorNormalized.clone());
            }
            getBestGenotypeFromPosterior(posteriorNormalized, cytosineAndAdjacent, i ,position);
            
		}
		
		tmpMethy = getMethyLevelFromPileup(pileup, ref);
		BisulfiteDiploidSNPGenotypeLikelihoods tmpGL = new BisulfiteDiploidSNPGenotypeLikelihoods(tracker, ref, (BisulfiteDiploidSNPGenotypePriors)priors, BAC, tmpMethy);
		tmpGL.setPriors(tracker, ref, BAC.heterozygosity, BAC.novelDbsnpHet, BAC.validateDbsnpHet, location);
		boolean firstSeen = true;
		for(String cytosineType : BAC.cytosineDefined.getContextDefined().keySet()){
			boolean heterozygousPattern = false; //heterozygous at cytosine psoition
			boolean heterozygousAtContext = false; //heterozygous at context psoition
			//tmpMethy = BAC.cytosineDefined.getContextDefined().get(cytosineType).cytosineMethylation;
			int cytosinePos = BAC.cytosineDefined.getContextDefined().get(cytosineType).cytosinePosition;
			//if(secondIteration)
			//	System.err.println(cytosineType + "\t" + tmpMethy);
            double adjacentCytosineSeqLikelihood = 0;
			double adjacentCytosineSeqLikelihoodReverseStrand = 0;
			int i = 1;
			int countMatchedOnFwd = 0;
			int countMatchedOnRvd = 0;
			//forward strand
			byte[] basesAlelleAFwd = cytosineType.getBytes();
			byte[] basesAlelleBFwd = cytosineType.getBytes();
			byte[] refWindFwd = new byte[cytosineType.length()];
            for(byte base : cytosineType.getBytes()){
            	int pos = i - cytosinePos;
            	int index = i-1;
            	refWindFwd[index] = refBytes[pos+maxCytosineLength];
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
	            		 if(BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(base, tmpMethyStatus.genotype.base1) && BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(base, tmpMethyStatus.genotype.base2)){//it is correct now, for CpH, if it is A/T heterozygouse SNP is still homozygous H
		 	            		countMatchedOnFwd++;
		 	            		adjacentCytosineSeqLikelihood += tmpMethyStatus.ratio;
		 	            		basesAlelleAFwd[index] = tmpMethyStatus.genotype.base1;
		 	            		basesAlelleBFwd[index] = tmpMethyStatus.genotype.base2;
		 	            }
	            		 else if(BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(base, tmpMethyStatus.genotype.base1) || BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(base, tmpMethyStatus.genotype.base2)){// it is not include into output now.. since CpG and CpA, CpA is not context anymore.. unless we implement numCAlleleA in CytosineParameters..
	            			//isHeterozygousInContextPosition
	            			 heterozygousAtContext = true;
	            			 basesAlelleAFwd[index] = tmpMethyStatus.genotype.base1;
		            		 basesAlelleBFwd[index] = tmpMethyStatus.genotype.base2;
		            		 countMatchedOnFwd++;
		            		 adjacentCytosineSeqLikelihood += tmpMethyStatus.ratio;
	            			 // heterozygousPattern = true;
	            			// basesAlelleAFwd[index] = tmpMethyStatus.genotype.base1;
	            			// basesAlelleBFwd[index] = tmpMethyStatus.genotype.base2;
	            			// countMatchedOnFwd++;
		 	            	// adjacentCytosineSeqLikelihood += tmpMethyStatus.ratio;
	            		 }
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
            //reverse strand
            byte[] basesAlelleARev = cytosineType.getBytes();
			byte[] basesAlelleBRev = cytosineType.getBytes();
			byte[] refWindRvd = new byte[cytosineType.length()];
            for(byte base : cytosineType.getBytes()){
            	int pos = cytosinePos - i;
            	int index = i-1;
            	refWindRvd[index] = refBytes[pos+maxCytosineLength];
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
	            		 if(BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(base, BaseUtilsMore.iupacCodeComplement(tmpMethyStatus.genotype.base1)) && BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(base, BaseUtilsMore.iupacCodeComplement(tmpMethyStatus.genotype.base2))){
	            			 countMatchedOnRvd++;
	            			 adjacentCytosineSeqLikelihoodReverseStrand += tmpMethyStatus.ratio;
	            			 basesAlelleARev[index] = tmpMethyStatus.genotype.base1;
	            			 basesAlelleBRev[index] = tmpMethyStatus.genotype.base2;
		 	            }
	            		 else if(BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(base, BaseUtilsMore.iupacCodeComplement(tmpMethyStatus.genotype.base1)) || BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(base, BaseUtilsMore.iupacCodeComplement(tmpMethyStatus.genotype.base2))){
	            			 heterozygousAtContext = true;
	            			 basesAlelleARev[index] = tmpMethyStatus.genotype.base1;
	            			 basesAlelleBRev[index] = tmpMethyStatus.genotype.base2;
	            			 countMatchedOnRvd++;
	            			 adjacentCytosineSeqLikelihoodReverseStrand += tmpMethyStatus.ratio;
	            			 
	            			 // heterozygousPattern = true;
	            			// basesAlelleARev[index] = tmpMethyStatus.genotype.base1;
	            			// basesAlelleBRev[index] = tmpMethyStatus.genotype.base2;
	            			// countMatchedOnRvd++;
	            			// adjacentCytosineSeqLikelihoodReverseStrand += tmpMethyStatus.ratio;
	            		 }
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
            if((countMatchedOnFwd < cytosineType.length() - 1) && (countMatchedOnRvd < cytosineType.length() - 1))
            	continue;
            
            
            //check at cytosine position now
            if(autoEstimateC && !secondIteration && !firstSeen){
            	
            }
            else{
            	firstSeen = false;
            	
            	tmpGL.clearLikelihoods(tmpMethy);
            	if(position == BAC.testLocus){
            		tmpGL.VERBOSE = true;
            		System.err.println("cytosineType: " + cytosineType );
            	}
            	
     			int nGoodBases = tmpGL.add(pileup, true, true);
                 if ( nGoodBases == 0 )
                     break;
                 double[] posteriorNormalized = MathUtils.normalizeFromLog10(tmpGL.getPosteriors(), true, false);
                 
                 getBestGenotypeFromPosterior(posteriorNormalized, cytosineAndAdjacent, 0 ,position);
            }
            
           
            
            methyStatus tmpMethyStatus = cytosineAndAdjacent.get(0);
        	if(tmpMethyStatus == null){
        		continue;
        	}
        	else if(tmpMethyStatus.genotype == null){
        		continue;
            }
            else if(tmpMethyStatus.genotype.isHet()){ //for CpG, if C is heterozygous, then marked it here.
            	if(BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(BaseUtils.C, tmpMethyStatus.genotype.base1) || BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(BaseUtils.C, tmpMethyStatus.genotype.base2)){
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
            		heterozygousPattern = true;
            		basesAlelleAFwd[cytosinePos-1] = tmpMethyStatus.genotype.base1;
            		basesAlelleBFwd[cytosinePos-1] = tmpMethyStatus.genotype.base2;
            		countMatchedOnFwd++;
            		adjacentCytosineSeqLikelihood += tmpMethyStatus.ratio;
            	}
            	else if(BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(BaseUtils.C, BaseUtilsMore.iupacCodeComplement(tmpMethyStatus.genotype.base1)) || BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(BaseUtils.C, BaseUtilsMore.iupacCodeComplement(tmpMethyStatus.genotype.base2))){
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
            		heterozygousPattern = true;
            		basesAlelleARev[cytosinePos-1] = BaseUtilsMore.iupacCodeComplement(tmpMethyStatus.genotype.base1);
            		basesAlelleBRev[cytosinePos-1] = BaseUtilsMore.iupacCodeComplement(tmpMethyStatus.genotype.base2);
            		countMatchedOnRvd++;
            		adjacentCytosineSeqLikelihoodReverseStrand += tmpMethyStatus.ratio;
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
            		
                	
            	}
            	
            }
        	
        	if(tmpMethyStatus.ratio > maxRatioInCytosinePos){
        		maxRatioInCytosinePos = tmpMethyStatus.ratio;
        		//bestCytosinePatternRelative = cytosineType;
        	}
            
            if(countMatchedOnFwd >= cytosineType.length()){
            	CytosineParameters cps = new CytosineParameters();
            	cps.isCytosinePattern = true;
            	cps.cytosineMethylation = tmpMethy;
				cps.cytosineStrand = '+';
				if(heterozygousPattern || heterozygousAtContext){
					cps.isHeterozygousCytosinePattern = heterozygousPattern;
					cps.isHeterozygousInContextPosition = heterozygousAtContext;
					cps.patternOfAlleleA = new String(basesAlelleAFwd);
					cps.patternOfAlleleB = new String(basesAlelleBFwd);
					
					if((basesAlelleAFwd[cytosinePos-1] == BaseUtils.C && basesAlelleBFwd[cytosinePos-1] == BaseUtils.T) ||(basesAlelleAFwd[cytosinePos-1] == BaseUtils.T && basesAlelleBFwd[cytosinePos-1] == BaseUtils.C)
							|| (basesAlelleAFwd[cytosinePos-1] == BaseUtils.G && basesAlelleBFwd[cytosinePos-1] == BaseUtils.A)|| (basesAlelleAFwd[cytosinePos-1] == BaseUtils.A && basesAlelleBFwd[cytosinePos-1] == BaseUtils.G)){
						cps.isCTHeterozygousLoci = true;
					}
					//System.err.println("fwd:" + basesAlelleAFwd[cytosinePos-1] + "\t" + basesAlelleBFwd[cytosinePos-1] + "\t" + cps.isCTHeterozygousLoci);
				}
				cps.isReferenceCytosinePattern = BaseUtilsMore.searchIupacPatternFromBases(cytosineType.getBytes(), refWindFwd, false);
				//value[0] = adjacentCytosineSeqLikelihood;
				if(adjacentCytosineSeqLikelihood > maxRatio){ 
					maxRatio = adjacentCytosineSeqLikelihood;
					
					bestCytosinePattern = cytosineType;
				}
				cytosineParametersStatus.put(cytosineType, cps);
				
			}
			else if(countMatchedOnRvd >= cytosineType.length()){
				CytosineParameters cps = new CytosineParameters();
            	cps.isCytosinePattern = true;
            	cps.cytosineMethylation = tmpMethy;
				cps.cytosineStrand = '-';
				if(heterozygousPattern || heterozygousAtContext){
					cps.isHeterozygousCytosinePattern = heterozygousPattern;
					cps.isHeterozygousInContextPosition = heterozygousAtContext;
					cps.patternOfAlleleA = new String(basesAlelleARev);
					cps.patternOfAlleleB = new String(basesAlelleBRev);
					
					if((basesAlelleARev[cytosinePos-1] == BaseUtils.C && basesAlelleBRev[cytosinePos-1] == BaseUtils.T) || (basesAlelleARev[cytosinePos-1] == BaseUtils.T && basesAlelleBRev[cytosinePos-1] == BaseUtils.C)
							|| (basesAlelleARev[cytosinePos-1] == BaseUtils.G && basesAlelleBRev[cytosinePos-1] == BaseUtils.A) || (basesAlelleARev[cytosinePos-1] == BaseUtils.A && basesAlelleBRev[cytosinePos-1] == BaseUtils.G)){
						cps.isCTHeterozygousLoci = true;
					}
					//System.err.println("fwd:" + basesAlelleARev[cytosinePos-1] + "\t" + basesAlelleBRev[cytosinePos-1] + "\t" + cps.isCTHeterozygousLoci);
				}
				cps.isReferenceCytosinePattern = BaseUtilsMore.searchIupacPatternFromBases(cytosineType.getBytes(), refWindRvd, true);
				if(adjacentCytosineSeqLikelihoodReverseStrand > maxRatio){
					maxRatio = adjacentCytosineSeqLikelihoodReverseStrand;
					
					bestCytosinePattern = cytosineType;
				}
				cytosineParametersStatus.put(cytosineType, cps);
					
			}
			else{
				CytosineParameters cps = new CytosineParameters();
				cps.isReferenceCytosinePattern = BaseUtilsMore.searchIupacPatternFromBases(cytosineType.getBytes(), refWindRvd, true);
				cps.isCytosinePattern = false;
				cytosineParametersStatus.put(cytosineType, cps);
			}
			if(position == BAC.testLocus){
            	System.err.println("countMatchedOnFwd: " + countMatchedOnFwd + "\tcountMatchedOnRvd: " + countMatchedOnRvd);
            } 
		}
		//if(bestCytosinePattern != null)
			return bestCytosinePattern;
		//else
		//	return bestCytosinePatternRelative;
		//return maxGL;	
	}
	
	//normalize of posterior
	/*
	public double[] normalization(double[] logPosterior){
		double sum = 0;
		double[] returnLikelyhood = logPosterior.clone();
		for(int i = 0; i < logPosterior.length; i++){
			sum += Math.pow(10,logPosterior[i]);
		}
		sum = Math.log10(sum);
		for(int j = 0; j < logPosterior.length; j++){
			returnLikelyhood[j] = returnLikelyhood[j] - sum;
		}
		return returnLikelyhood;
	}
*/
	private void initializeBestAndAlternateAlleleFromPosterior(double[] posterior, int location){
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

	private void getBestGenotypeFromPosterior(double[] posterior,HashMap<Integer,methyStatus> cytosineAdjacent, int key, int location){
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
        tmpMethyStatus.genotype = bestGenotype;
        cytosineAdjacent.put(key, tmpMethyStatus);
        
      
	}


	

	
	//inner class to record genotype and posterior ratio of best and second best genotype
	private class methyStatus{
		DiploidGenotype genotype;
		double ratio;
		methyStatus(){
			
		}
	}
	
	 private int getFilteredDepth(ReadBackedPileup pileup) {
	        int count = 0;
	        for ( PileupElement p : pileup ) {
	            if ( BaseUtils.isRegularBase( p.getBase() ) )
	                count += p.getRepresentativeCount();
	        }

	        return count;
	    }
	 
	 private double getMethyLevelFromPileup(ReadBackedPileup pileup, ReferenceContext ref){
		 int numCNegStrand = 0;
		 int numTNegStrand = 0;
		 int numCPosStrand = 0;
		 int numTPosStrand = 0;
		 for ( PileupElement p : pileup ) {
         	SAMRecord samRecord = p.getRead();
         	boolean negStrand = samRecord.getReadNegativeStrandFlag();
         	//byte observedBase = p.getBase();
         	if(samRecord.getDuplicateReadFlag()){ //get rid of duplicate reads
             	continue;
             }
         	int offset = p.getOffset();
         	if(offset < 0)//is deletion
         		continue;
         	boolean paired = samRecord.getReadPairedFlag();
         	if(paired){
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
 					continue;
 				}
 				
 				// Inverted dups, count only one end
 				if (samRecord.getAlignmentStart() == samRecord.getMateAlignmentStart() && samRecord.getReadNegativeStrandFlag() == samRecord.getMateNegativeStrandFlag())
 				{
					if(BAC.invDups == INVERT_DUPS.USE_ONLY_1ST_END){
						if (samRecord.getSecondOfPairFlag()) continue;
					}
					else if(BAC.invDups == INVERT_DUPS.NOT_TO_USE){
						continue;
					}
					else{
						
					}
    			}
 	        	if (Paired  && !BAC.USE_BADLY_MATED_READS && !samRecord.getProperPairFlag())
 				{
 					continue;
 				}
 	        	
 	        
 	        	if(secondOfPair){	        		
 		        	//samRecord.setReadNegativeStrandFlag(!samRecord.getReadNegativeStrandFlag());  
 	        		//observedBase = BaseUtils.simpleComplement(observedBase);
 	        		negStrand = !negStrand;
 	        	}
 	        	
         	}
				
	        	
				
				GATKSAMRecordFilterStorage GATKrecordFilterStor = new GATKSAMRecordFilterStorage((GATKSAMRecord)p.getRead(), BAC, ref, p.getOffset());
             
				if(GATKrecordFilterStor.isGoodBase() && isGoodConvertedBase(p.getRead(), offset)){
				//if(GATKrecordFilterStor.isGoodBase()){
					if(negStrand){
						if(p.getBase()==BaseUtils.G){
							numCNegStrand++;
						}
						else if(p.getBase()==BaseUtils.A){
							numTNegStrand++;
						}
						
						
					}
					else{
						if(p.getBase()==BaseUtils.C){
							numCPosStrand++;
						}
						else if(p.getBase()==BaseUtils.T){
							numTPosStrand++;
						}
						
					}
				}

         }
		 double methyLevel = (double)(numCNegStrand + numCPosStrand)/(double)(numCNegStrand + numTNegStrand + numCPosStrand + numTPosStrand);
		 if((numCNegStrand + numTNegStrand + numCPosStrand + numTPosStrand) == 0)
			 return 0;
		 else
			 return methyLevel;
		 
	 }
	
	 private boolean isGoodConvertedBase(GATKSAMRecord GATKrecord, int offset){
		 int start = (Integer)GATKrecord.getTemporaryAttribute(CYTOSINE_CONV_START);
		 if( start != -1 && offset >=  (Integer)GATKrecord.getTemporaryAttribute(CYTOSINE_CONV_START))
				return true;
			else
				return false;
	 }
	 
	 
	private void setGoodConvertedBase(GATKSAMRecord GATKrecord, ReferenceContext refContext, boolean negStrand){
			if(GATKrecord.containsTemporaryAttribute(CYTOSINE_CONV_START)){
				
				//System.err.println(GATKrecord.getTemporaryAttribute(CYTOSINE_CONV_START));
				return;
			}
		int cytosineConvertStart = -1;
			byte[] bases = GATKrecord.getReadBases();
			if(refContext == null)
				return;
			//System.err.println(new String(refBases));
			//byte[] refBases = refContext.getBases();
			int convertedCount=0;
			
			if(negStrand){
				for(int i=0; i<bases.length;i++){
					GenomeLoc loc = refContext.getGenomeLocParser().createGenomeLoc(refContext.getLocus().getContig(), GATKrecord.getAlignmentStart()+i);
					if( !refContext.getWindow().containsP(loc) )
						break;
					
					ReferenceContext tmpRef = new ReferenceContext(refContext.getGenomeLocParser(),loc, refContext.getWindow(),refContext.getBases());
					byte refBase = tmpRef.getBase();
					//System.err.print((char)refBase);
					if(BaseUtils.basesAreEqual(refBase, BaseUtils.G)){
						if(BaseUtils.basesAreEqual(bases[i], BaseUtils.A)){
							convertedCount++;
							if(convertedCount >= BAC.minConv){
								cytosineConvertStart=i;
								break;
							}
								
						}
					}
				}
			}
			else{
				for(int i=0; i<bases.length;i++){
					//System.err.println(GATKrecord.getAlignmentStart() + "\t" + i + "\t" + refContext.getWindow().getStart() + "\t" + refContext.getWindow().getStop());
					GenomeLoc loc = refContext.getGenomeLocParser().createGenomeLoc(refContext.getLocus().getContig(), GATKrecord.getAlignmentStart()+i);
					if( !refContext.getWindow().containsP(loc) )
						break;
					
					ReferenceContext tmpRef = new ReferenceContext(refContext.getGenomeLocParser(),loc, refContext.getWindow(),refContext.getBases());
					byte refBase = tmpRef.getBase();
					//System.err.print((char)refBase);
					if(BaseUtils.basesAreEqual(refBase, BaseUtils.C)){
						if(BaseUtils.basesAreEqual(bases[i], BaseUtils.T)){
							convertedCount++;
							if(convertedCount >= BAC.minConv){  // this will cause problem in non-ref C's methylation estimation, since convert on C here is depended on ref C.
								cytosineConvertStart=i;
								break;
							}
								
						}
					}
				}
			}
			GATKrecord.setTemporaryAttribute(CYTOSINE_CONV_START, cytosineConvertStart);
			
		}
	
	
	 public ReadBackedPileup createBAQedPileup( final ReadBackedPileup pileup ) {
	        final List<PileupElement> BAQedElements = new ArrayList<PileupElement>();
	        for( final PileupElement PE : pileup ) {
	            final PileupElement newPE = new BAQedPileupElement( PE );
	            BAQedElements.add( newPE );
	        }
	        return new ReadBackedPileupImpl( pileup.getLocation(), BAQedElements );
	    }

	    public class BAQedPileupElement extends PileupElement {
	        public BAQedPileupElement( final PileupElement PE ) {
	            super(PE.getRead(), PE.getOffset(), PE.isDeletion(), PE.isBeforeDeletion(), PE.isBeforeInsertion(), PE.isNextToSoftClip());
	        }

	        @Override
	        public byte getQual( final int offset ) { return BisBAQ.calcBAQFromTag(getRead(), offset, true); }
	    }

}
