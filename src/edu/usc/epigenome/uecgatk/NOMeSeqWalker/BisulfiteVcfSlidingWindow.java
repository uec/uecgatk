/**
 * within five CpG dinucleotides (each CpG dinucleotide contains two cytosines,
one on each strand). Only those cytosines covered by at least three cytosine
or thymine reads were counted, and each CpG dinucleotide was assigned a
weighting factor defined as the span (in bp) between the next CpG dinucleotide
upstream and the next CpG dinucleotide downstream. A weighted average
was calculated for each window, and those windows with an average DNA
methylation of less than 5% in both tumor and adjacent normal tissue were
categorized as methylation resistant. Those windows with methylation of less
than 5% in the adjacent normal tissue and greater than 35% in the tumor were
characterized as methylation prone, and those windows with methylation of
greater than 35% in the adjacent normal tissue and less than 5% in the tumor
were characterized as methylation loss. Two or more overlapping regions from
a single methylation class were merged into one. For enrichment of functional
annotations within these regions (Fig. 3b,c), elements of each methylation
class within 500 bp were merged into a single locus. Elements were not merged
for motif analysis (Fig. 3d).
 */
package edu.usc.epigenome.uecgatk.NOMeSeqWalker;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import org.broadinstitute.sting.utils.GenomeLoc;

import edu.usc.epigenome.uecgatk.BisSNP.BisulfiteVCFConstants;
import edu.usc.epigenome.uecgatk.BisSNP.BisulfiteVariantCallContext;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Feb 28, 2013 12:05:58 PM
 * 
 */
public class BisulfiteVcfSlidingWindow {
	
	private HashMap<String, CytosineStatus> cytosineFilters = null; //record the different cytosine status in this window
	private String mainCytosines = "HCG";
	private int windowLenDefined = 200;
	private int windowLen = 0;
	
	private LinkedList<BisulfiteVariantCallContext> window = null;
	
	public BisulfiteVcfSlidingWindow(Integer minCTs, Integer minCs, int windowLenDefined, String mainCytosines){
		
	}
	
	public BisulfiteVcfSlidingWindow(List<String> cytosinePattern, List<Integer> minCTs, List<Integer> minCs, int windowLenDefined, String mainCytosines){ //the main cytosine pattern to construct window
		this.windowLenDefined = windowLenDefined;
		this.mainCytosines = mainCytosines;
		cytosineFilters = new HashMap<String, CytosineStatus>();
		if(cytosinePattern.size() != minCTs.size() ||cytosinePattern.size() != minCs.size() ){
			throw new IllegalArgumentException("Number of cytosines defined:" + cytosinePattern.size() + " is different from the number of minCT criteria: " + minCTs.size() +"or number of minCs criteria: " + minCs.size());
		}
		for(int i = 0; i < cytosinePattern.size() ; i ++){
			CytosineStatus tmp = new CytosineStatus(minCTs.get(i), minCs.get(i));
			cytosineFilters.put(cytosinePattern.get(i),tmp);
		}
		window = new LinkedList<BisulfiteVariantCallContext>();
	}
	
	public BisulfiteVcfSlidingWindow(List<String> cytosinePattern, int minCT, int minNumC, int windowLenDefined, String mainCytosines){ //the main cytosine pattern to construct window
		this.windowLenDefined = windowLenDefined;
		this.mainCytosines = mainCytosines;
		cytosineFilters = new HashMap<String, CytosineStatus>();

		for(int i = 0; i < cytosinePattern.size() ; i++){
			CytosineStatus tmp = new CytosineStatus(minCT, minNumC);
			cytosineFilters.put(cytosinePattern.get(i),tmp);
		}
		window = new LinkedList<BisulfiteVariantCallContext>();
	}
	
	public boolean add(BisulfiteVariantCallContext bvc, boolean slidingWindow){ //return true, if it is really added into the window
		
		if(bvc != null && bvc.getVariantContext() != null && !bvc.getVariantContext().getGenotypes().isEmpty()){
			String cytosineStr = bvc.getVariantContext().getGenotype(0).getAttributeAsString(BisulfiteVCFConstants.BEST_C_PATTERN, ".");
			//System.err.println(cytosineStr + "\t" + cytosineFilters.containsKey(cytosineStr) + "\t" + bvc.getSummaryAcrossRG().numC + "\t" + bvc.getSummaryAcrossRG().numT);
			if(cytosineFilters.containsKey(cytosineStr) && (bvc.getSummaryAcrossRG().numC + bvc.getSummaryAcrossRG().numT) >= cytosineFilters.get(cytosineStr).minCT){

				add(cytosineStr, bvc);

				while(getWindowLen()> this.windowLenDefined && getCytosineNum(mainCytosines) > cytosineFilters.get(mainCytosines).minCs){
					pop();
				}
				return true;
				
			}
			
			
		}
		return false;
	}
	
	
	public void add(String cytosinePattern, BisulfiteVariantCallContext bvc){
		if(!window.isEmpty() && !window.peekLast().ref.getLocus().onSameContig(bvc.ref.getLocus())){
			window.clear();
			windowLen = 0;
			clearCytosineStat();
		}
		
			window.offerLast(bvc);

		windowLen = bvc.ref.getLocus().distance(window.peekFirst().ref.getLocus());
		CytosineStatus cs = cytosineFilters.get(cytosinePattern);
		
		cs.numCytosinePattern++;
		cs.numC += bvc.getSummaryAcrossRG().numC;
		cs.numT += bvc.getSummaryAcrossRG().numT;
		double methy = (double)bvc.getSummaryAcrossRG().numC/(double)(bvc.getSummaryAcrossRG().numC + bvc.getSummaryAcrossRG().numT);
		cs.totalMethyCytosinePattern += methy;
		cs.methyVector.offerLast(methy);
		cytosineFilters.put(cytosinePattern, cs);
	}
	
	public void pop(){
		BisulfiteVariantCallContext bvcPopOut = window.pollFirst();
		String cytosineStr = bvcPopOut.getVariantContext().getGenotype(0).getAttributeAsString(BisulfiteVCFConstants.BEST_C_PATTERN, ".");
		windowLen = window.peekLast().ref.getLocus().distance(window.peekFirst().ref.getLocus());
		CytosineStatus cs = cytosineFilters.get(cytosineStr);
		cs.numCytosinePattern--;
		cs.numC -= bvcPopOut.getSummaryAcrossRG().numC;
		cs.numT -= bvcPopOut.getSummaryAcrossRG().numT;
		double methy = (double)bvcPopOut.getSummaryAcrossRG().numC/(double)(bvcPopOut.getSummaryAcrossRG().numC + bvcPopOut.getSummaryAcrossRG().numT);
		cs.totalMethyCytosinePattern -= methy;
		cs.methyVector.pollFirst();
		cytosineFilters.put(cytosineStr, cs);

	}
	
	
	
	public int getWindowLen(){
		return windowLen;
	}
	
	public LinkedList<BisulfiteVariantCallContext> getWindow(){
		return window;
	}
	
	public double[] getMethyVector(String cytosinePattern){
		if(cytosineFilters.containsKey(cytosinePattern)){
			double[] vector =new double[cytosineFilters.get(cytosinePattern).methyVector.size()];
			int i = 0;
			for(double v : cytosineFilters.get(cytosinePattern).methyVector){
				vector[i] = v;
				i++;
			}
			return vector;
		}
		return null;
	}
	
	public double getCytosineMethyMean(String cytosinePattern){
		
		return getCytosineMethySum(cytosinePattern)/getCytosineNum(cytosinePattern);
	}
	
	public double getCytosineMethyWeightMean(String cytosinePattern){
		int preBvcDist = 0;
		BisulfiteVariantCallContext preBvc = null;
		double methySum = 0.0;
		GenomeLoc startPos = null;
		for(BisulfiteVariantCallContext bvc : window){
			if(bvc.getVariantContext().getGenotype(0).getAttributeAsString(BisulfiteVCFConstants.BEST_C_PATTERN, ".").equalsIgnoreCase(cytosinePattern)){
				if(preBvc == null){
					preBvc = bvc;
					preBvcDist = 0;
					startPos = bvc.ref.getLocus();
				}
				else{
					int dist = preBvc.ref.getLocus().distance(bvc.ref.getLocus());
					methySum += ((double)(preBvcDist + dist)/2.0) * (double)preBvc.getSummaryAcrossRG().numC/(double)(preBvc.getSummaryAcrossRG().numC + preBvc.getSummaryAcrossRG().numT);
					preBvc = bvc;
					preBvcDist = dist;
				}
			}
			
		}
		if(preBvc != null){
			methySum += ((double)preBvcDist/2.0) * (double)preBvc.getSummaryAcrossRG().numC/(double)(preBvc.getSummaryAcrossRG().numC + preBvc.getSummaryAcrossRG().numT);
			methySum /= preBvc.ref.getLocus().distance(startPos);
		}
		else{
			methySum = Double.NaN;
		}
		
		return methySum;
	}
	
	public int getCytosineNum(String cytosinePattern){
		if(cytosineFilters.containsKey(cytosinePattern)){
			return cytosineFilters.get(cytosinePattern).numCytosinePattern;
		}
		return 0;
	}
	
	public int getCytosineCTReads(String cytosinePattern){
		if(cytosineFilters.containsKey(cytosinePattern)){
			return (cytosineFilters.get(cytosinePattern).numC + cytosineFilters.get(cytosinePattern).numT);
		}
		return 0;
	}
	
	public int getCytosineCReads(String cytosinePattern){
		if(cytosineFilters.containsKey(cytosinePattern)){
			return cytosineFilters.get(cytosinePattern).numC;
		}
		return 0;
	}
	
	public double getCytosineMethySum(String cytosinePattern){
		if(cytosineFilters.containsKey(cytosinePattern)){
			return cytosineFilters.get(cytosinePattern).totalMethyCytosinePattern;
		}
		return Double.NaN;
	}
	
	public void clearCytosineStat(){
		for(String pattern : cytosineFilters.keySet()){
			CytosineStatus cs = cytosineFilters.get(pattern);
			cs.methyVector.clear();
			cs.numCytosinePattern = 0;
			cs.numC = 0;
			cs.numT = 0;
			cs.totalMethyCytosinePattern = 0;
			cytosineFilters.put(pattern, cs);
		}
	}
	
	public class CytosineStatus{
		private int minCT = 3;
		private int minCs = 5;  //minimum number of Cs in 
		
		public int numC = 0;
		public int numT = 0;
		
		public int numCytosinePattern = 0;
		public double totalMethyCytosinePattern = 0;
		
		LinkedList<Double> methyVector = null;

		 public CytosineStatus(int minCT, int minCs){
			 methyVector = new LinkedList<Double>();
			 this.minCT = minCT;
			 this.minCs = minCs;
		 }
	}
}
