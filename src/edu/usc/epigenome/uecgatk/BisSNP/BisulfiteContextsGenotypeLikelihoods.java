/**
 * 
 */
package edu.usc.epigenome.uecgatk.BisSNP;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;

import net.sf.samtools.SAMRecord;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.variantcontext.Allele;

import edu.usc.epigenome.uecgatk.BisSNP.BisulfiteEnums.INVERT_DUPS;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time 2012 Mar 19, 2012 3:38:42 PM
 * 
 */
public class BisulfiteContextsGenotypeLikelihoods {

	private String sample;
	private int numOfCReadsInBisulfiteCStrand = 0;
	private int numOfTReadsInBisulfiteCStrand = 0;
	private int numOfOtherReadsInBisulfiteCStrand = 0;
	private int numOfGReadsInGenotypeGStrand = 0;
	private int numOfAReadsInGenotypeGStrand = 0;
	private int numOfOtherReadsInGenotypeGStrand = 0;
	private int totalDepth = 0;
	private String bestMatchedCytosinePattern;
	private HashMap<String,CytosineParameters> cytosineParameters;
	//private int cytosinePos;
	private Set<String> cytosineContexts;
	
	private HashMap<Integer,double[]> GPsBeforeCytosineTenGenotypes;
	private HashMap<Integer,double[]> GPsAfterCytosineTenGenotypes;
	private HashMap<String, Double[]> GPsAtCytosineTenGenotypes;
	private double[] GPsAtCytosineNormalizedByThreeGenotypes;
    private Allele A, B;
    private ReadBackedPileup pileup;
    private BisulfiteArgumentCollection BAC;
    private int readsFwdRef = 0;
    private int readsFwdAlt = 0;
    private int readsRevRef = 0;
    private int readsRevAlt = 0;
    
    private int readsAlleleA = 0;
    private int readsAlleleB = 0;
    private double totalBaseQualAlleleA = Double.NaN;
    private double totalBaseQualAlleleB = Double.NaN;
    private double rmsBaseQual = Double.NaN;
    private double rmsMapQual = Double.NaN;
    private int mapQual0 = 0;
    
    private ReferenceContext ref = null;
  //Reads supporting ALT. Number of 1) forward ref alleles; 2) reverse ref; 3) forward non-ref; 4) reverse non-ref alleles"
    //Average base quality for reads supporting alleles. For each allele, in the same order as listed
	/**
	 * 
	 */
	public BisulfiteContextsGenotypeLikelihoods(String sample, Allele A,
			Allele B, double log10aaPosteriorLikelihoods, double log10abPosteriorLikelihoods,
			double log10bbPosteriorLikelihoods,Set<String> cytosineContexts, int numOfCReadsInBisulfiteCStrand,int numOfTReadsInBisulfiteCStrand,
			int numOfOtherReadsInBisulfiteCStrand, int numOfGReadsInGenotypeGStrand, int numOfAReadsInGenotypeGStrand, int numOfOtherReadsInGenotypeGStrand,
			int totalDepth, HashMap<String,CytosineParameters> cytosineParameters, String bestMatchedCytosinePattern, HashMap<Integer,double[]> GPsBeforeCytosineTenGenotypes, HashMap<Integer,double[]> GPsAfterCytosineTenGenotypes,
			HashMap<String, Double[]> GPsAtCytosineTenGenotypes, ReadBackedPileup pileup, BisulfiteArgumentCollection BAC ) {
		// TODO Auto-generated constructor stub
		this.sample = sample;
		this.numOfCReadsInBisulfiteCStrand = numOfCReadsInBisulfiteCStrand;
		this.numOfTReadsInBisulfiteCStrand = numOfTReadsInBisulfiteCStrand;
		this.numOfOtherReadsInBisulfiteCStrand = numOfOtherReadsInBisulfiteCStrand;
		this.numOfGReadsInGenotypeGStrand = numOfGReadsInGenotypeGStrand;
		this.numOfAReadsInGenotypeGStrand = numOfAReadsInGenotypeGStrand;
		this.numOfOtherReadsInGenotypeGStrand = numOfOtherReadsInGenotypeGStrand;
		this.totalDepth = totalDepth;
		this.bestMatchedCytosinePattern = bestMatchedCytosinePattern;
		this.cytosineParameters = cytosineParameters;
		//this.cytosinePos;
		this.cytosineContexts = cytosineContexts;
		
		this.GPsBeforeCytosineTenGenotypes = GPsBeforeCytosineTenGenotypes;
		this.GPsAfterCytosineTenGenotypes = GPsAfterCytosineTenGenotypes;
		this.GPsAtCytosineTenGenotypes = GPsAtCytosineTenGenotypes;
		//this.GPsAtCytosineNormalizedByThreeGenotypes = GPsAtCytosineNormalizedByThreeGenotypes;
		this.GPsAtCytosineNormalizedByThreeGenotypes = new double[]{log10aaPosteriorLikelihoods, log10abPosteriorLikelihoods, log10bbPosteriorLikelihoods};
	    this.A = A;
	    this.B = B;
		this.pileup = pileup;
		this.BAC = BAC;
		setupSupportingReadsInfo();
	}

	public String getSample() {
        return sample;
    }

    public double getAALikelihoods() {
        return GPsAtCytosineNormalizedByThreeGenotypes[0];
    }

    public double getABLikelihoods() {
        return GPsAtCytosineNormalizedByThreeGenotypes[1];
    }

    public double getBBLikelihoods() {
        return GPsAtCytosineNormalizedByThreeGenotypes[2];
    }

    public double[] getLikelihoods() {
        return GPsAtCytosineNormalizedByThreeGenotypes;
    }

    public Allele getAlleleA() {
        return A;
    }

    public Allele getAlleleB() {
        return B;
    }

    public int getDepth() {
        return totalDepth;
    }
    
    public double getMethylationLevel(){
    	double methy = (double)numOfCReadsInBisulfiteCStrand/(double)(numOfCReadsInBisulfiteCStrand + numOfTReadsInBisulfiteCStrand);

    	return methy;
    }
    
    public int getNumOfCReadsInBisulfiteCStrand(){
    	return numOfCReadsInBisulfiteCStrand;
    }
	
    public int getNumOfTReadsInBisulfiteCStrand(){
    	return numOfTReadsInBisulfiteCStrand;
    }
    
    public int getNumOfOtherReadsInBisulfiteCStrand(){
    	return numOfOtherReadsInBisulfiteCStrand;
    }
	
    public int getNumOfGReadsInGenotypeGStrand(){
    	return numOfGReadsInGenotypeGStrand;
    }
	
    public int getNumOfAReadsInGenotypeGStrand(){
    	return numOfAReadsInGenotypeGStrand;
    }
    
    public int getNumOfOtherReadsInGenotypeGStrand(){
    	return numOfOtherReadsInGenotypeGStrand;
    }
    
    public String getBaseCountStatusAsString(){
    	
    	String cStatus = numOfCReadsInBisulfiteCStrand + "," + numOfTReadsInBisulfiteCStrand + "," + numOfOtherReadsInBisulfiteCStrand + "," + numOfGReadsInGenotypeGStrand + "," + numOfAReadsInGenotypeGStrand + "," + numOfOtherReadsInGenotypeGStrand;
    	return cStatus;
    }
    
    public String getBestMatchedCytosinePattern(){
    	return bestMatchedCytosinePattern;
    }
    
    public HashMap<String,CytosineParameters> getCytosineParameters(){
    	return cytosineParameters;
    }
    
    public HashMap<Integer,double[]> getGPsBeforeCytosineTenGenotypes(){
    	return GPsBeforeCytosineTenGenotypes;
    }
    
    public void setGPsBeforeCytosineTenGenotypes(HashMap<Integer,double[]> GPsBeforeCytosineTenGenotypes){
    	this.GPsBeforeCytosineTenGenotypes = (HashMap<Integer, double[]>)(GPsBeforeCytosineTenGenotypes.clone());
    }
    
    public HashMap<Integer,double[]> getGPsAfterCytosineTenGenotypes(){
    	return GPsAfterCytosineTenGenotypes;
    }
    
    public void setGPsAfterCytosineTenGenotypes(HashMap<Integer,double[]> GPsAfterCytosineTenGenotypes){
    	this.GPsAfterCytosineTenGenotypes = (HashMap<Integer, double[]>)(GPsAfterCytosineTenGenotypes.clone());
    }
    
    public HashMap<String, Double[]> getGPsAtCytosineTenGenotypes(){
    	return GPsAtCytosineTenGenotypes;
    }
    
    public String getAveBaseQualAsString(){
    	String aveBaseQual = null;
    	if(!Double.isNaN(totalBaseQualAlleleA/readsAlleleA) && !Double.isNaN(totalBaseQualAlleleA/readsAlleleA)){
    		aveBaseQual = String.format("%.1f", totalBaseQualAlleleA/readsAlleleA) + "," + String.format("%.1f", totalBaseQualAlleleB/readsAlleleB);
    	}
    	else if(!Double.isNaN(totalBaseQualAlleleA/readsAlleleA)){
    		aveBaseQual =String.format("%.1f", totalBaseQualAlleleA/readsAlleleA);
    	}
    	else if(!Double.isNaN(totalBaseQualAlleleB/readsAlleleB)){
    		aveBaseQual = String.format("%.1f", totalBaseQualAlleleB/readsAlleleB);
    	}
    	else{
    		aveBaseQual = VCFConstants.MISSING_VALUE_v4;
    	}
    	return aveBaseQual;
    }
    
    public String getDP4AsString(){
    	String DP4String = readsFwdRef + "," + readsRevRef + "," + readsFwdAlt + "," + readsRevAlt;
    	return DP4String;
    } 
    
    public double getBQ(){
    	return rmsBaseQual;
    }
    
    public double getMQ(){
    	return rmsMapQual;
    }
    
    public int getMQ0(){
    	return mapQual0;
    }
	
    private void setupSupportingReadsInfo(){
    	int[] rmsBaseQualTotal = new int[pileup.depthOfCoverage()];
    	int[] rmsMapQualTotal = new int[pileup.depthOfCoverage()];
    	int index = 0;
    	
    	for ( PileupElement p : pileup ) {
        	SAMRecord samRecord = p.getRead();
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
		        	samRecord.setReadNegativeStrandFlag(!samRecord.getReadNegativeStrandFlag());        		
	        	}
	        	
        	}
			
        	boolean negStrand = samRecord.getReadNegativeStrandFlag();
			
			
			 //Reads supporting ALT. Number of 1) forward ref alleles; 2) reverse ref; 3) forward non-ref; 4) reverse non-ref alleles"
            //Average base quality for reads supporting alleles. For each allele, in the same order as listed
        	rmsMapQualTotal[index] = p.getMappingQual();
			rmsBaseQualTotal[index] = p.getQual();
			index++;
			GATKSAMRecordFilterStorage GATKrecordFilterStor = new GATKSAMRecordFilterStorage((GATKSAMRecord)p.getRead(), BAC, ref, p.getOffset());
            //GATKrecordFilterStor.setGoodBases(badReadPileupFilter, true);
			if(p.getMappingQual()==0){
				//if(mapQual0==-1){
				//	mapQual0=1;
				//}
				//else{
					mapQual0++;
				//}
			}
			if(GATKrecordFilterStor.isGoodBase()){
				//should different by GPs
				
				
				if(negStrand){

					if(BaseUtilsMore.iupacCodeEqual(p.getBase(), A.getBases()[0], negStrand)){
						readsAlleleA++;
						
						if(Double.isNaN(totalBaseQualAlleleA)){
							totalBaseQualAlleleA = p.getQual();
						}
						else{
							totalBaseQualAlleleA += p.getQual();
						}
						if(A.isReference()){
							readsRevRef++;
						}
						else{
							readsRevAlt++;
						}
					}
					else if(BaseUtilsMore.iupacCodeEqual(p.getBase(), B.getBases()[0], negStrand)){
						readsAlleleB++;
						
						if(Double.isNaN(totalBaseQualAlleleB)){
							totalBaseQualAlleleB = p.getQual();
						}
						else{
							totalBaseQualAlleleB += p.getQual();
						}
						if(B.isReference()){
							readsRevRef++;
						}
						else{
							readsRevAlt++;
						}
					}

				}
				else{
					if(BaseUtilsMore.iupacCodeEqual(p.getBase(), A.getBases()[0], negStrand)){
						readsAlleleA++;
						
						if(Double.isNaN(totalBaseQualAlleleA)){
							totalBaseQualAlleleA = p.getQual();
						}
						else{
							totalBaseQualAlleleA += p.getQual();
						}
						if(A.isReference()){
							readsFwdRef++;
						}
						else{
							readsFwdAlt++;
						}
					}
					else if(BaseUtilsMore.iupacCodeEqual(p.getBase(), B.getBases()[0], negStrand)){
						readsAlleleB++;
						
						if(Double.isNaN(totalBaseQualAlleleB)){
							totalBaseQualAlleleB = p.getQual();
						}
						else{
							totalBaseQualAlleleB += p.getQual();
						}
						if(B.isReference()){
							readsFwdRef++;
						}
						else{
							readsFwdAlt++;
						}
					}
				}
			}
			
        }
    	
    	rmsBaseQual = MathUtils.rms(rmsBaseQualTotal);
    	rmsMapQual = MathUtils.rms(rmsMapQualTotal);
    }
    
    
}
