/**
 * 
 */
package edu.usc.epigenome.uecgatk.BisSNP;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time 2012 Mar 20, 2012 11:07:32 AM
 * 
 */
public class CytosineParameters {

	/**
	 * 
	 */
	public int cytosinePosition;
	public double cytosineMethylation;
	public char cytosineStrand;
	public boolean methylationAutoestimated;
	public boolean isCytosinePattern = false;
	
	//implement heterozygous, need to do it in the future, useful for Allele specific analysis.. 
	public boolean isHeterozygousCytosinePattern = false;
	public int numOfCReadsInCytosinePosInBisulfiteCStrandAlleleA = 0;
	public int numOfTReadsInCytosinePosInBisulfiteCStrandAlleleA = 0;
	public int numOfOtherReadsInCytosinePosInBisulfiteCStrandAlleleA = 0;
	public int numOfCReadsInCytosinePosInBisulfiteCStrandAlleleB = 0;
	public int numOfTReadsInCytosinePosInBisulfiteCStrandAlleleB = 0;
	public int numOfOtherReadsInCytosinePosInBisulfiteCStrandAlleleB = 0;
	public String patternOfAlleleA = null;
	public String patternOfAlleleB = null;
	
	public CytosineParameters() {
		// TODO Auto-generated constructor stub
	}
	
		

}
