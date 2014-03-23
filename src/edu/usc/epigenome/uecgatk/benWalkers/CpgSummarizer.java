/**
 * 
 */
package edu.usc.epigenome.uecgatk.benWalkers;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.LinkedList;
import java.util.List;

import javax.swing.text.html.HTML;

import com.googlecode.charts4j.Color;

import edu.usc.epigenome.genomeLibs.MatUtils;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;
import edu.usc.epigenome.uecgatk.benWalkers.cytosineWalkers.Cpg;


/**
 * @author benb
 *
 * The client streams CpGs into the streamCpg function, and the summarizer can
 * give various statistics.
 * 
 */
public abstract class CpgSummarizer implements Cloneable {
	
	protected MethylDbQuerier querier = null;
//	protected String sampleName = "UNKNOWN_SAMPLE";
//	protected String featName = "UNKNOWN_FEAT";
	protected int nCpgsSeen = 0;
	
	// Arbitrary summary statistics for a single value.
	protected double numVals = 0.0;
	protected double valsTotal = 0.0;
	
	protected double valsWeightingTotal = 0.0;
	protected double weightingTotal = 0.0;
	protected double valsSquareTotal = 0.0;
	protected double valsMin = Double.NaN;
	protected double valsMax = Double.NaN;
	
	// Used for chart drawing
	protected String label = "UNKNOWN_SUMMARIZER";
	protected String desc = "NO_DESC";
	protected double naturalMin = Double.NaN;
	protected double naturalMax = Double.NaN;
	protected Color color = Color.BLACK;
	
//	protected Cpg lastCpgSeen = null;
//	int curWeight = -1;
	
	/**
	 * We need the querier criteria for certain output
	 */
	public CpgSummarizer() {
		super();
		init(new MethylDbQuerier());
	}

	
	
	/**
	 * We need the querier criteria for certain output
	 */
	public CpgSummarizer(MethylDbQuerier inQuerier) {
		super();
		init(inQuerier);
	}

//	/**
//	 * We need the querier criteria for certain output
//	 */
//	public CpgSummarizer(MethylDbQuerier inQuerier, String sampleName, String featName) {
//		super();
//		this.setFeatName(featName);
//		this.setSampleName(sampleName);
//		init(inQuerier);
//	}

	// This should be called by all constructors
	public void init(MethylDbQuerier inQuerier)
	{
		this.querier = inQuerier;
		
		this.setDesc(String.format("minReads=%d nonconvFilt=%s, maxOppAfrac=%f",
				querier.getMinCTreads(), ""+querier.getUseNonconversionFilter(), querier.getMaxOppstrandAfrac()));

	}
	
	public void reset()
	{
		nCpgsSeen = 0;
		numVals = 0.0;
		valsTotal = 0.0;
		valsSquareTotal = 0.0;
		valsMin = Double.NaN;
		valsMax = Double.NaN;

//		this.lastCpgSeen = null;
		this.valsWeightingTotal = 0.0;
		this.weightingTotal = 0.0;
		
		
	}
	
	
	

	@Override
	protected Object clone() throws CloneNotSupportedException {
		return super.clone();
	}



	static public CpgSummarizer sumSummarizers(CpgSummarizer a, CpgSummarizer b)
	{
		//System.err.printf("Adding summarizers: %s,%s\n",a,b);
		CpgSummarizer out = null;
		try {
			out = a.getClass().cast(a.clone());
			out.nCpgsSeen += b.nCpgsSeen;
			out.numVals += b.numVals;
			out.valsTotal += b.valsTotal;
			out.valsWeightingTotal += b.valsWeightingTotal;
			out.weightingTotal += b.weightingTotal;
			out.valsSquareTotal += b.valsSquareTotal;
			out.valsMin = MatUtils.nanMin(a.valsMin, b.valsMin);
			out.valsMax = MatUtils.nanMax(a.valsMax, b.valsMax);
			out.naturalMin = MatUtils.nanMin(a.naturalMin, b.naturalMin);
			out.naturalMax = MatUtils.nanMax(a.naturalMax, b.naturalMax);
			} 
		catch (CloneNotSupportedException e) 
		{
			System.err.printf("FAILED adding summarizers: %s,%s\n",a,b);
			e.printStackTrace();
		}

		return out;
	}
	
	/**
	 * Override this!
	 * @param cpg
	 */
	public void streamCpg(Cpg cpg)
	{
		nCpgsSeen++;
		
//		// We'll use this weighting when a value is captured
//		if (this.lastCpgSeen != null)
//		{
//			
//			// If the last CpG is simply the opposite strand one,
//			// the dist will be 1.  In this case, use the weighting to
//			// the previous one.
//			int dist = cpg.chromPos - this.lastCpgSeen.chromPos;
//			if (dist == 1)
//			{
//				// Do nothing
//			}
//			else
//			{
//				this.curWeight = dist;
//
//				if (this.curWeight > 10000)
//				{
//					// If it's >10000, it's a jump between elements generally.  Just count it as one
//					this.curWeight = 1; 
//					//System.err.printf("%d\tDistance between two CpGs>100000\t%d\t%d\n",curWeight, lastCpgSeen.chromPos, cpg.chromPos);
//				}
//				else if (this.curWeight < 0)
//				{
//					// This happens at the beginning of the chromosome when the last cpg is from the last chrom
//					this.curWeight = 1;
//					//System.err.printf("%d\tDistance between two CpGs<0\t%d\t%d\n",curWeight, lastCpgSeen.chromPos, cpg.chromPos);
//				}
//			}
//
//			//System.err.println("Setting cpg weight: " + this.curWeight);
//			cpg.setCpgWeight(this.curWeight);
//
//		}
//		this.lastCpgSeen = cpg;
		
	}
	
	abstract public void removeCpg(Cpg cpg);
		
	
	/**
	 * Use this to update summary stats
	 * @param val
	 */
	public void streamValue(double val, double weight)
	{
		if (!Double.isNaN(val))
		{
			this.numVals++;
			this.valsTotal += val;
			this.valsSquareTotal += Math.pow(val, 2.0);
			if (Double.isNaN(this.valsMin) || (val < this.valsMin)) this.valsMin = val;
			if (Double.isNaN(this.valsMax) || (val > this.valsMax)) this.valsMax = val;
			
			// Keep track of one version weighted by base pairs covered by the CpG
			if (weight > 0)
			{
				//System.err.println("\t\tAdding weight: " + curWeight);
				this.weightingTotal += weight;
				this.valsWeightingTotal += (weight * val);
				if (Double.isNaN(this.weightingTotal) || Double.isNaN(this.valsWeightingTotal))
				{
					System.err.printf("Running totals weighting: %f\tvalsWeightingTotal=%f\n", this.weightingTotal, this.valsWeightingTotal);
				}
				
			}
		}
		
//		System.err.printf("Summarizer streaming, val =%.2f\n",this.getValMean(false));

	}

//	protected void removeValue(double val)
//	{
//		removeValue(val, Double.NaN);
//	}
	
	protected void removeValue(double val, double weight)
	{
		if (!Double.isNaN(val))
		{
			this.numVals--;
			this.valsTotal -= val;
			this.valsSquareTotal -= Math.pow(val, 2.0);
			
			if (!Double.isNaN(weight) && (weight > 0))
			{
				//System.err.println("\t\tRemoving weight: " + weight);
				this.weightingTotal -= weight;
				this.valsWeightingTotal -= (weight * val);
			}
		}
	}

	public String getDesc() {
		return desc;
	}



	public void setDesc(String desc) {
		this.desc = desc;
	}



	public String getSampleName() {
		return this.querier.getMethylTablePrefix();
	}

//	public void setSampleName(String sampleName) {
//		this.sampleName = sampleName;
//	}

	public String getFeatName() {
		return querier.getFeatFilterTypesList();
	}

//	public void setFeatName(String featName) {
//		this.featName = featName;
//	}

	public int getnCpgsSeen() {
		return nCpgsSeen;
	}

	public double getValMean()
	{
		return getValMean(false);
	}
	
	/**
	 * @param useSpaceWeightedMean If this is set, the value for each Cpg is weighted by the amount of 
	 *  sequence covered by the CpG.
	 * @return
	 */
	public double getValMean(boolean useSpaceWeightedMean)
	{
		return (useSpaceWeightedMean) ? (this.valsWeightingTotal / this.weightingTotal) : (this.valsTotal / this.numVals);
	}
	
	public double getValStdev()
	{
		return Math.sqrt((valsSquareTotal/numVals)-Math.pow(this.getValMean(),2));
	}
	
	public double getValsTotal() {
		return valsTotal;
	}

	
	public double getNumVals() {
		return numVals;
	}

	public double getValsSquareTotal() {
		return valsSquareTotal;
	}

	public double getValsMin() {
		return valsMin;
	}

	public double getValsMax() {
		return valsMax;
	}

	public double getNaturalMin() {
		return naturalMin;
	}

	public double getNaturalMax() {
		return naturalMax;
	}

	public String getLabel() {
		return label;
	}

	public void setLabel(String label) {
		this.label = label;
	}

	public Color getColor() {
		return color;
	}

	public void setColor(Color color) {
		this.color = color;
	}

	public void setNaturalMin(double naturalMin) {
		this.naturalMin = naturalMin;
	}

	public void setNaturalMax(double naturalMax) {
		this.naturalMax = naturalMax;
	}
	

	/**** HTML tables , Static ****/
	
	public static String htmlTableStart()
	{
		StringWriter sw = new StringWriter(500);
		PrintWriter pw = new PrintWriter(sw);
		
		pw.print("<TABLE BORDER=1>");
		return sw.toString();
	}
		
	public static String htmlTableFinish()
	{
		StringWriter sw = new StringWriter(500);
		PrintWriter pw = new PrintWriter(sw);
		
		pw.print("</TABLE>");
		return sw.toString();
	}

	public static String htmlTable(List<CpgSummarizer> summarizers)
	{
		StringWriter sw = new StringWriter(100000);
		PrintWriter pw = new PrintWriter(sw);
		
		pw.print(htmlTableStart());
		for (CpgSummarizer summ : summarizers)
		{
			pw.print(summ.htmlTableRow());
		}
		pw.print(htmlTableFinish());
		
		return sw.toString();
	}
	

	/**** HTML tables , non-static ****/
	public String htmlTableRow()
	{
		StringWriter sw = new StringWriter(30000);
		PrintWriter pw = new PrintWriter(sw);

		pw.println("<TR>");
		pw.printf("<TD>%s</TD>\n", this.getLabel());
		pw.printf("<TD>%s</TD>\n", this.getSampleName());
		pw.printf("<TD>%s</TD>\n", this.getFeatName());
		pw.printf("<TD>%d</TD>\n", this.getnCpgsSeen());
		pw.printf("<TD>%.3f</TD>\n", this.getValMean());
		pw.printf("<TD>%.3f</TD>\n", this.getValMean(true));
		pw.printf("<TD>%.3f</TD>\n", this.getValStdev());
		pw.printf("<TD>%s</TD>\n", this.getDesc());
		pw.print("</TR>");
		
		
		return sw.toString();
	}
	
	

}
