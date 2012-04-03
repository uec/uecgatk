/**
 * 
 */
package edu.usc.epigenome.uecgatk.BisSNP;

import java.util.HashMap;
import java.util.Set;

import org.broadinstitute.sting.utils.variantcontext.Allele;

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
	/**
	 * 
	 */
	public BisulfiteContextsGenotypeLikelihoods(String sample, Allele A,
			Allele B, double log10aaPosteriorLikelihoods, double log10abPosteriorLikelihoods,
			double log10bbPosteriorLikelihoods,Set<String> cytosineContexts, int numOfCReadsInBisulfiteCStrand,int numOfTReadsInBisulfiteCStrand,
			int numOfOtherReadsInBisulfiteCStrand, int numOfGReadsInGenotypeGStrand, int numOfAReadsInGenotypeGStrand, int numOfOtherReadsInGenotypeGStrand,
			int totalDepth, HashMap<String,CytosineParameters> cytosineParameters, String bestMatchedCytosinePattern, HashMap<Integer,double[]> GPsBeforeCytosineTenGenotypes, HashMap<Integer,double[]> GPsAfterCytosineTenGenotypes,
			HashMap<String, Double[]> GPsAtCytosineTenGenotypes ) {
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
    	return (double)numOfCReadsInBisulfiteCStrand/(double)(numOfCReadsInBisulfiteCStrand + numOfTReadsInBisulfiteCStrand);
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
	
    
}
