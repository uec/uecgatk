package org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel;

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
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.genotype.DiploidGenotype;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
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

public class BisulfiteDiploidSNPGenotypeLikelihoods implements Cloneable  {
	protected BisulfiteDiploidSNPGenotypePriors priors = null;

	protected BisulfiteArgumentCollection BAC;
    protected double BISULFITE_CONVERSION_RATE;
    protected double OVER_CONVERSION_RATE;
    protected String DETERMINED_CYTOSINE_TYPE_POS = ""; //record cytosine pattern in positive strand which owns the maximum posterior probability
    protected double DETERMINED_CYTOSINE_TYPE_C_LIKELIHOOD_POS = Double.NEGATIVE_INFINITY; //record cytosine pattern likelihood in positive strand
    protected double DETERMINED_CYTOSINE_TYPE_C_METHY_POS = 0; //record cytosine pattern methylation level in positive strand which owns the maximum posterior probability
    protected String DETERMINED_CYTOSINE_TYPE_NEG = ""; //record cytosine pattern in negative strand which owns the maximum posterior probability
    protected double DETERMINED_CYTOSINE_TYPE_C_LIKELIHOOD_NEG = Double.NEGATIVE_INFINITY; //record cytosine pattern likelihood in negative strand
    protected double DETERMINED_CYTOSINE_TYPE_C_METHY_NEG = 0; //record cytosine pattern methylation level in negative strand which owns the maximum posterior probability

    public double PROB_OF_REFERENCE_ERROR = 1e-6; //record the reference genome error
    
    protected double[] log10Likelihoods = null;
    protected double[] log10Posteriors = null;
    
    public boolean VERBOSE = false;
   
	public BisulfiteDiploidSNPGenotypeLikelihoods(
			RefMetaDataTracker tracker, ReferenceContext ref, BisulfiteDiploidSNPGenotypePriors priors, BisulfiteArgumentCollection BAC, double[] cytosineMethyStatus) {
		this.priors = priors;
		this.BAC = BAC;
		this.BISULFITE_CONVERSION_RATE = BAC.bsRate;
		this.OVER_CONVERSION_RATE = BAC.overRate;
		this.PROB_OF_REFERENCE_ERROR = BAC.referenceGenomeErr;
        this.DETERMINED_CYTOSINE_TYPE_C_METHY_POS = cytosineMethyStatus[0];
        this.DETERMINED_CYTOSINE_TYPE_C_METHY_NEG = cytosineMethyStatus[1];
		
	}
	
	/**
     * set prior calculated
     */
	public void setPriors(RefMetaDataTracker tracker, ReferenceContext ref, double heterozygousity, double novelDbsnpHet, double validateDbsnpHet, GenomeLoc loc){
		this.priors.setPriors(tracker, ref, heterozygousity, PROB_OF_REFERENCE_ERROR, novelDbsnpHet, validateDbsnpHet, loc, BAC.tiVsTv);
        setToZeroBs();
	}
	
	/**
     * clear the likelihood to primary status. this could save time on prior calculation.
     */
	public void clearLikelihoods(double[] cytosineMethyStatus){
		
		this.DETERMINED_CYTOSINE_TYPE_C_METHY_POS = cytosineMethyStatus[0];
        this.DETERMINED_CYTOSINE_TYPE_C_METHY_NEG = cytosineMethyStatus[1];
        setToZeroBs();
	}
	
	/**
     * add pileup into likelihood calculation.
     */
	public int add(ReadBackedPileup pileup, boolean ignoreBadBases, boolean capBaseQualsAtMappingQual) {
        int n = 0;
        
        ReadBackedPileup pileupPositiveStrand = pileup.getPositiveStrandPileup();
        for(PileupElement p : pileupPositiveStrand )
        	n += add(p, ignoreBadBases, capBaseQualsAtMappingQual, false);
        
        ReadBackedPileup pileupNegativeStrand = pileup.getNegativeStrandPileup();
        for(PileupElement p : pileupNegativeStrand )
        	n += add(p, ignoreBadBases, capBaseQualsAtMappingQual, true);
        
        if ( VERBOSE ) {
        	System.out.println("summary:");
        	for ( DiploidGenotype g : DiploidGenotype.values() ) { System.out.printf("%s\t", g); }
            System.out.println();
            for ( DiploidGenotype g : DiploidGenotype.values() ) { System.out.printf("%.2f\t", log10Likelihoods[g.ordinal()]); }
            System.out.println();
            for ( DiploidGenotype g : DiploidGenotype.values() ) { System.out.printf("%.2f\t", getPriors()[g.ordinal()]); }
            System.out.println();
        }
 
        return n;
    }
	

	/**
     * add pileupelement into likelihood calculation. and call likelihood calculation function
     */
	
	public int add(PileupElement p, boolean ignoreBadBases, boolean capBaseQualsAtMappingQual, boolean negStrand) {
        
        byte observedBase = 0, qualityScore = 0;
        if ( !usableBase(p, ignoreBadBases) ){
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
        	if (Paired  && !BAC.USE_BADLY_MATED_READS && !samRecord.getProperPairFlag())
    		{
        		return 0;
    		}
        	if(secondOfPair){
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

	/**
     * the main likelihood calculation function for any input base with quality information and strand information
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
					if(BaseUtils.basesAreEqual(observedBase, BaseUtils.C)){
						likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) * ( (1.0-DETERMINED_CYTOSINE_TYPE_C_METHY_POS) * (1.0-BISULFITE_CONVERSION_RATE) + DETERMINED_CYTOSINE_TYPE_C_METHY_POS * (1.0-OVER_CONVERSION_RATE)) : pOfBase1 * (error/3.0);
						likelihood += observedBase == g.base2 ? pOfBase2 * (1.0-error) *  ( (1.0-DETERMINED_CYTOSINE_TYPE_C_METHY_POS) * (1.0-BISULFITE_CONVERSION_RATE) + DETERMINED_CYTOSINE_TYPE_C_METHY_POS * (1.0-OVER_CONVERSION_RATE) ) : pOfBase2 * (error/3.0);	
						if ( VERBOSE ) {
						//	System.out.println("flag1: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
						}
					}
					else if(BaseUtils.basesAreEqual(observedBase, BaseUtils.T)){
						likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) : (g.base1 == BaseUtils.C ? pOfBase1 * (error/3.0 + (1.0-DETERMINED_CYTOSINE_TYPE_C_METHY_POS) * BISULFITE_CONVERSION_RATE + DETERMINED_CYTOSINE_TYPE_C_METHY_POS * OVER_CONVERSION_RATE) : pOfBase1 * (error/3.0));
						likelihood += observedBase == g.base2 ? pOfBase1 * (1.0-error) : (g.base2 == BaseUtils.C ? pOfBase2 * (error/3.0 + (1.0-DETERMINED_CYTOSINE_TYPE_C_METHY_POS) * BISULFITE_CONVERSION_RATE + DETERMINED_CYTOSINE_TYPE_C_METHY_POS * OVER_CONVERSION_RATE) : pOfBase2 * (error/3.0));

						if ( VERBOSE ) {
						//	System.out.println("flag3: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
						}
					}
					else{
						likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) : pOfBase1 * (error/3.0);
						likelihood += observedBase == g.base2 ? pOfBase2 * (1.0-error) : pOfBase2 * (error/3.0);
						if ( VERBOSE ) {
							//System.out.println("flag5: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
						}
					}
				}
				else{
						if(BaseUtils.basesAreEqual(observedBase, BaseUtils.G)){
							likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) * ( (1.0-DETERMINED_CYTOSINE_TYPE_C_METHY_NEG) * (1.0-BISULFITE_CONVERSION_RATE) + DETERMINED_CYTOSINE_TYPE_C_METHY_NEG * (1.0-OVER_CONVERSION_RATE) ) : pOfBase1 * (error/3.0);
							likelihood += observedBase == g.base2 ? pOfBase2 * (1.0-error) *  ( (1.0-DETERMINED_CYTOSINE_TYPE_C_METHY_NEG) * (1.0-BISULFITE_CONVERSION_RATE) + DETERMINED_CYTOSINE_TYPE_C_METHY_NEG * (1.0-OVER_CONVERSION_RATE) ) : pOfBase2 * (error/3.0);
							if ( VERBOSE ) {
								//System.out.println("flag6: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
							}
							
						}
						else if(BaseUtils.basesAreEqual(observedBase, BaseUtils.A)){
							likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) : (g.base1 == BaseUtils.G ? pOfBase1 * (error/3.0 + (1.0-DETERMINED_CYTOSINE_TYPE_C_METHY_NEG) * BISULFITE_CONVERSION_RATE + DETERMINED_CYTOSINE_TYPE_C_METHY_NEG * OVER_CONVERSION_RATE) : pOfBase1 * (error/3.0));
							likelihood += observedBase == g.base2 ? pOfBase1 * (1.0-error) : (g.base2 == BaseUtils.G ? pOfBase2 * (error/3.0 + (1.0-DETERMINED_CYTOSINE_TYPE_C_METHY_NEG) * BISULFITE_CONVERSION_RATE + DETERMINED_CYTOSINE_TYPE_C_METHY_NEG * OVER_CONVERSION_RATE) : pOfBase2 * (error/3.0));
							if ( VERBOSE ) {
								//System.out.println("flag9: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
							}
						}
						else{
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
	
	/**
     * identify the base is able to use, copy from GATK since private method there
     */
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
     *   Criterion : observed base isn't a A,C,T,G or lower case equivalent
     *
     * @param observedBase observed base
     * @return true if the base is a bad base
     */
    protected static boolean badBase(byte observedBase) {
        return BaseUtils.simpleBaseToBaseIndex(observedBase) == -1;
    }
	
    public double[] getPriors() {
        return this.priors.getPriors();
    }

    protected void setToZeroBs() {
    	log10Likelihoods = new double[DiploidGenotype.values().length];
    	for ( DiploidGenotype g : DiploidGenotype.values() ) {
            int i = g.ordinal();
            log10Likelihoods[i] = log10(1-PROB_OF_REFERENCE_ERROR);
        }
        log10Posteriors = this.priors.getPriors().clone();     // posteriors are all set to the priors
    }

    protected Object clone() throws CloneNotSupportedException {
    	BisulfiteDiploidSNPGenotypeLikelihoods c = (BisulfiteDiploidSNPGenotypeLikelihoods)super.clone();
        c.priors = priors;
        c.log10Likelihoods = log10Likelihoods.clone();
        c.log10Posteriors = log10Posteriors.clone();
        return c;
    }

}
