package edu.usc.epigenome.uecgatk.bisulfitesnpmodel;

import java.util.Map;

import net.sf.samtools.SAMRecord;

import org.apache.log4j.Logger;
import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel.BisulfiteDiploidSNPGenotypeLikelihoods;
import org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel.BisulfiteDiploidSNPGenotypePriors;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext.StratifiedContextType;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.genotyper.BiallelicGenotypeLikelihoods;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypePriors;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.genotype.DiploidGenotype;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

public class BisulfiteSNPGenotypeLikelihoodsCalculationModel extends
		NonRefDependSNPGenotypeLikelihoodsCalculationModel {

	protected Byte bestAllele = null;
	protected Byte alternateAllele = null;
	protected final boolean useAlleleFromVCF;
	protected long testLoc;
	protected static double BISULFITE_CONVERSION_RATE;
	protected static double CPG_METHYLATION_RATE = 0.75;
	
	public BisulfiteSNPGenotypeLikelihoodsCalculationModel(
			UnifiedArgumentCollection UAC, Logger logger) {
		super(UAC, logger);
		useAlleleFromVCF = UAC.GenotypingMode == GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES;
		this.testLoc = UAC.testLocus;
		BISULFITE_CONVERSION_RATE = UAC.bsRate;
		CPG_METHYLATION_RATE = UAC.CpgMethy;
		// TODO Auto-generated constructor stub
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


        byte[] refWindow = ref.getBasesAtLocus(2);
        //System.err.println(refWindow.length);
        byte refBase = refWindow[0];
        byte refNextBase = refWindow[1];
       // byte refNextBase = refWindow[0];
        Allele refAllele = Allele.create(refBase, true);
        //this.testLoc = UAC.testLocus;
        //System.out.println(this.testLoc);

        // find the best allele and alternative allele with the largest sum of quality scores
        if ( alternateAlleleToUse != null ) {
            bestAllele = alternateAlleleToUse.getBases()[0];
            alternateAllele = alternateAlleleToUse.getBases()[0];
        } else if ( useAlleleFromVCF ) {
            final VariantContext vcInput = tracker.getVariantContext(ref, "alleles", null, ref.getLocus(), true);
            if ( vcInput == null )
                return null;
            if ( !vcInput.isSNP() ) {
                logger.info("Record at position " + ref.getLocus() + " is not a SNP; skipping...");
                return null;
            }
            if ( !vcInput.isBiallelic() ) {
                logger.info("Record at position " + ref.getLocus() + " is not bi-allelic; choosing the first allele...");
                //return null;
            }
            bestAllele = vcInput.getAlternateAllele(0).getBases()[0];
            alternateAllele = vcInput.getAlternateAllele(0).getBases()[0];
        } else {
            //initializeBestAndAlternateAllele(refWindow, contexts);
        	

        }
        
      /*
       // System.err.println(refBase+"\t" + bestAllele+"\t" + alternateAllele + "\t" + refWindow[0] + "\t" + refWindow[1]);
     // if there are no non-ref bases...
        if ( (alternateAllele == null && bestAllele == refBase) || (bestAllele == null) ) {
            // if we only want variants, then we don't need to calculate genotype likelihoods
            if ( UAC.OutputMode == UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_VARIANTS_ONLY )
                return refAllele;

            // otherwise, choose any alternate allele (it doesn't really matter)
            bestAllele = (byte)(refBase != BaseUtils.A ? BaseUtils.A : BaseUtils.C);
            alternateAllele = (byte)(refBase != BaseUtils.A ? BaseUtils.G : BaseUtils.T);
        }
        
        Allele AlleleA, AlleleB;
        
        if(alternateAllele == null || alternateAllele == refBase || alternateAllele == bestAllele){
        	AlleleA = Allele.create(refBase, true);
        	AlleleB = Allele.create(bestAllele, false);
        	//if(alternateAllele == null)
        		alternateAllele = bestAllele;
        	bestAllele = refBase;
        }
        else if(bestAllele == refBase){
        	AlleleA = Allele.create(bestAllele, true);
        	AlleleB = Allele.create(alternateAllele, false);
        }
        else{
        	AlleleA = Allele.create(bestAllele, false);
        	AlleleB = Allele.create(alternateAllele, false);
        }
*/
       //System.err.println(refBase+"\t"+bestAllele+"\t"+alternateAllele);
        

        for ( Map.Entry<String, StratifiedAlignmentContext> sample : contexts.entrySet() ) {
            ReadBackedPileup pileup = sample.getValue().getContext(contextType).getBasePileup();

            for ( PileupElement p : pileup ) {
            	SAMRecord samRecord = p.getRead();
            	
            	boolean negStrand = samRecord.getReadNegativeStrandFlag();
				int alignmentS = samRecord.getAlignmentStart();
				int	onRefCoord = (negStrand) ? samRecord.getUnclippedEnd() : alignmentS;

				if((pileup.getLocation().getStart()) == testLoc){
					System.out.println("before filter:\t" + onRefCoord + "\t" + p.getOffset() + "\t" + negStrand + "\t" + pileup.getLocation().getStart() + "\t" + (char)p.getBase());
					//System.out.println("deletion: " + p.isDeletion());
					//System.out.println("GATKSAMRecord: " + (p.getRead() instanceof GATKSAMRecord));
					System.out.println("isGoodBase: " + ((GATKSAMRecord)p.getRead()).isGoodBase(p.getOffset()));
		                     
				}
            }
            
            
            // do not use this prior, this prior is flat prior intiated in genotypeEngine, so we actually do not transfer this priors...
            BisulfiteDiploidSNPGenotypeLikelihoods GL = new BisulfiteDiploidSNPGenotypeLikelihoods(tracker, ref, (BisulfiteDiploidSNPGenotypePriors)priors, UAC.PCR_error, UAC.bsRate, UAC.CpgMethy);
            if((pileup.getLocation().getStart()) == testLoc)
            	GL.VERBOSE=true;
            int nGoodBases = GL.add(pileup, true, true, refNextBase);
            if ( nGoodBases == 0 )
                continue;

            double[] likelihoods = GL.getLikelihoods();
            double[] posterior = GL.getPosteriors();
            double[] prio = GL.getPriors();
            
            initializeBestAndAlternateAlleleFromPosterior(posterior, pileup.getLocation().getStart());
            
            if ( (alternateAllele == null && bestAllele == refBase) || (bestAllele == null) ) {
                // if we only want variants, then we don't need to calculate genotype likelihoods
                if ( UAC.OutputMode == UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_VARIANTS_ONLY )
                    return refAllele;

                // otherwise, choose any alternate allele (it doesn't really matter)
                bestAllele = (byte)(refBase != BaseUtils.A ? BaseUtils.A : BaseUtils.C);
                alternateAllele = (byte)(refBase != BaseUtils.A ? BaseUtils.G : BaseUtils.T);
            }
            
            Allele AlleleA, AlleleB;
            
            if(alternateAllele == null || alternateAllele == refBase || alternateAllele == bestAllele){
            	AlleleA = Allele.create(refBase, true);
            	AlleleB = Allele.create(bestAllele, false);
            	//if(alternateAllele == null)
            		alternateAllele = bestAllele;
            	bestAllele = refBase;
            }
            else if(bestAllele == refBase){
            	AlleleA = Allele.create(bestAllele, true);
            	AlleleB = Allele.create(alternateAllele, false);
            }
            else{
            	AlleleA = Allele.create(bestAllele, false);
            	AlleleB = Allele.create(alternateAllele, false);
            }

            DiploidGenotype AAGenotype = DiploidGenotype.createHomGenotype(bestAllele);
            DiploidGenotype ABGenotype = DiploidGenotype.createDiploidGenotype(bestAllele, alternateAllele);
            DiploidGenotype BBGenotype = DiploidGenotype.createHomGenotype(alternateAllele);
            if((pileup.getLocation().getStart()) == testLoc){
            	System.out.println("sample: " + sample.getKey());
            	System.out.println("sample location: " + pileup.getPileupString((char)refBase));
            	System.out.println("sample: " + sample.getValue().getLocation().getStart());
            	System.out.println("refBase: " + refBase + "refNextBase: " + refNextBase + " bestAllele: " + bestAllele + "alternateAllele: " + alternateAllele);
            	System.out.println("nGoodBases " + nGoodBases);
            	System.out.println("AAGenotype " + likelihoods[AAGenotype.ordinal()] + "\t" + prio[AAGenotype.ordinal()] + "\t" + posterior[AAGenotype.ordinal()]);
            	System.out.println("ABGenotype " + likelihoods[ABGenotype.ordinal()] + "\t" + prio[ABGenotype.ordinal()] + "\t" + posterior[ABGenotype.ordinal()]);
            	System.out.println("BBGenotype " + likelihoods[BBGenotype.ordinal()] + "\t" + prio[BBGenotype.ordinal()] + "\t" + posterior[BBGenotype.ordinal()]);
            }
            
            	GLs.put(sample.getKey(), new BiallelicGenotypeLikelihoods(sample.getKey(),
            			AlleleA,
            			AlleleB,
                        likelihoods[AAGenotype.ordinal()],
                        likelihoods[ABGenotype.ordinal()],
                        likelihoods[BBGenotype.ordinal()],
                        getFilteredDepth(pileup)));

            
        }

        return refAllele;
	}
	
	protected void initializeBestAndAlternateAlleleFromPosterior(double[] posterior, int location){
		int[] qualCounts = new int[4];
		
        
        for ( byte altAllele : BaseUtils.BASES ) {
        	for ( DiploidGenotype g : DiploidGenotype.values() ){
    			if(g.base1 == altAllele){
    				int index = BaseUtils.simpleBaseToBaseIndex(altAllele);
	                if ( index >= 0 )
	                    qualCounts[index] += posterior[g.ordinal()];
    			}
    			if(g.base2 == altAllele){
    				int index = BaseUtils.simpleBaseToBaseIndex(altAllele);
	                if ( index >= 0 )
	                    qualCounts[index] += posterior[g.ordinal()];
    			}
    		}
        }
        
        int maxCount = Integer.MIN_VALUE;
        int secondMaxCount = Integer.MIN_VALUE;
        bestAllele = null;
        alternateAllele = null;
        for ( byte altAllele : BaseUtils.BASES ) {
            int index = BaseUtils.simpleBaseToBaseIndex(altAllele);
            if ( qualCounts[index] > maxCount ) {
            	secondMaxCount = maxCount;
            	maxCount = qualCounts[index];
            	if(bestAllele != null){
            		alternateAllele = bestAllele;
            	}
                bestAllele = altAllele;
            }
            else if (qualCounts[index] > secondMaxCount && qualCounts[index] <= maxCount){
            	secondMaxCount = qualCounts[index];
            	alternateAllele = altAllele;
            }
            //System.err.println();
        }
		
        if(location == testLoc){
        	for ( DiploidGenotype g : DiploidGenotype.values() ){
        		System.err.println(g.base1 + "-" + g.base2 + ": " + posterior[g.ordinal()]);
        	}
        	for ( byte altAllele : BaseUtils.BASES ) {
        		int index = BaseUtils.simpleBaseToBaseIndex(altAllele);
        		System.err.println(altAllele + ": " + qualCounts[index]);
        	}
        	System.err.println("bestAllele: " + bestAllele + "\t" + maxCount);
        	if(alternateAllele != null){
        		System.err.println("AlternateAllele: " + "\t" + alternateAllele + "\t" + secondMaxCount);
        	}
        }
	}
	
	protected void initializeBestAndAlternateAllele(byte[] refWindow, Map<String, StratifiedAlignmentContext> contexts) {
        int[] qualCounts = new int[4];
        byte ref = refWindow[0];
        byte refNextBase = refWindow[1];
        int location = 0;
        for ( Map.Entry<String, StratifiedAlignmentContext> sample : contexts.entrySet() ) {
            // calculate the sum of quality scores for each base
            ReadBackedPileup pileup = sample.getValue().getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getBasePileup();
            for ( PileupElement p : pileup ) {
                // ignore deletions and filtered bases
            	SAMRecord samRecord = p.getRead();
            	
            	boolean negStrand = samRecord.getReadNegativeStrandFlag();
				int alignmentS = samRecord.getAlignmentStart();
				int	onRefCoord = (negStrand) ? samRecord.getUnclippedEnd() : alignmentS;

				/*if((pileup.getLocation().getStart()) == testLoc){
					System.out.println("before filter:\t" + onRefCoord + "\t" + p.getOffset() + "\t" + negStrand + "\t" + pileup.getLocation().getStart() + "\t" + (char)p.getBase());
					System.out.println("deletion: " + p.isDeletion());
					System.out.println("GATKSAMRecord: " + (p.getRead() instanceof GATKSAMRecord));
					System.out.println("isGoodBase: " + ((GATKSAMRecord)p.getRead()).isGoodBase(p.getOffset()));
		                     
				}*/
                if ( p.isDeletion() ||
                     (p.getRead() instanceof GATKSAMRecord && !((GATKSAMRecord)p.getRead()).isGoodBase(p.getOffset())) )
                    continue;

                	//if((p.getRead().getReadNegativeStrandFlag() && p.getBase() == BaseUtils.A) || (!p.getRead().getReadNegativeStrandFlag() && p.getBase() == BaseUtils.T))
                	//	continue;

                	
                	
					if((pileup.getLocation().getStart()) == testLoc){
						System.out.println(onRefCoord + "\t" + p.getOffset() + "\t" + negStrand + "\t" + pileup.getLocation().getStart() + "\t" + (char)p.getBase());
						location = pileup.getLocation().getStart();
					}
					
				if((p.getRead().getReadNegativeStrandFlag() && refNextBase == BaseUtils.C) || (!p.getRead().getReadNegativeStrandFlag() && refNextBase == BaseUtils.G)){
					 int index = BaseUtils.simpleBaseToBaseIndex(p.getBase());
		                if ( index >= 0 )
		                    qualCounts[index] += p.getQual();
				}
				else{
					//need to recheck! sincer 95% conversion of C do not means 95% back conversion of T
					if((p.getRead().getReadNegativeStrandFlag() && ref == BaseUtils.A) || (!p.getRead().getReadNegativeStrandFlag() && ref == BaseUtils.T)){
						if(p.getRead().getReadNegativeStrandFlag()){
							int indexA = BaseUtils.simpleBaseToBaseIndex(BaseUtils.A);
							int indexG = BaseUtils.simpleBaseToBaseIndex(BaseUtils.G);
			                qualCounts[indexG] += p.getQual() * BISULFITE_CONVERSION_RATE;
			                qualCounts[indexA] += p.getQual() * (1-BISULFITE_CONVERSION_RATE);
						}
						else{
							int indexC = BaseUtils.simpleBaseToBaseIndex(BaseUtils.C);
							int indexT = BaseUtils.simpleBaseToBaseIndex(BaseUtils.T);
			                qualCounts[indexC] += p.getQual() * BISULFITE_CONVERSION_RATE;
			                qualCounts[indexT] += p.getQual() * (1-BISULFITE_CONVERSION_RATE);
						}
					}
					else{
						int index = BaseUtils.simpleBaseToBaseIndex(p.getBase());
		                if ( index >= 0 )
		                    qualCounts[index] += p.getQual();
						
					}
					
				}
               
            }
        }

        // set the best base with maximum quality score sum
        int maxCount = 0;
        int secondMaxCount = 0;
        bestAllele = null;
        alternateAllele = null;
        for ( byte altAllele : BaseUtils.BASES ) {
            //if ( altAllele == ref )
              //  continue;
            //if( altAllele == 'T')
            	//continue;
            int index = BaseUtils.simpleBaseToBaseIndex(altAllele);
            if ( qualCounts[index] > maxCount ) {
            	secondMaxCount = maxCount;
            	maxCount = qualCounts[index];
            	if(bestAllele != null){
            		alternateAllele = bestAllele;
            	}
                bestAllele = altAllele;
            }
            else if (qualCounts[index] > secondMaxCount && qualCounts[index] < maxCount){
            	secondMaxCount = qualCounts[index];
            	alternateAllele = altAllele;
            }
           // System.err.println();
        }
        
        if(location == testLoc){
        	System.out.println("refAllele: " + ref);
        	System.out.println("bestAllele: " + bestAllele + "\t" + maxCount);
        	if(alternateAllele != null){
        		System.out.println("AlternateAllele: " + "\t" + alternateAllele + "\t" + secondMaxCount);
        	}
        }

    }
	
}
