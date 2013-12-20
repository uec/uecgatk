/**
 * 
 */
package edu.usc.epigenome.uecgatk.benWalkers;

import java.io.PrintWriter;
import java.util.List;

import org.usckeck.genome.ChromFeatures;

import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerParams;
import edu.usc.epigenome.uecgatk.benWalkers.cytosineWalkers.Cpg;


/**
 * @author benb
 * 
 * Writes an optional bedGraph file, which is better for UCSC.
 *
 */
public class CpgWalkerWindowWigWriter extends CpgWalker {

	protected PrintWriter wigWriter = null; 
	protected PrintWriter bedgraphWriter = null; 
	protected boolean csvMode = false;
	protected int curChrNum = -1;
	
	private static int[][] regionBuffer = new int[10000][2];
	
	// Tracking for bedgraph
	int lastWindMidpoint = -1;
	int lastSpanEnd = -1;
	double lastMeth = -1.0;
	

	public CpgWalkerWindowWigWriter(CpgWalkerParams inWalkParams) {
		super(inWalkParams);
	}



	
	//////////////////////////
	//	
	// Overridden CpgWalker
	//	
	//////////////////////////
	
	/**
	 * @return the writer
	 */
	public PrintWriter getWigWriter() {
		return wigWriter;
	}

	/**
	 * @param writer the writer to set
	 */
	public void setWigWriter(PrintWriter writer) {
		this.wigWriter = writer;
	
	}
	
	
	
	/**
	 * @return the csvMode
	 */
	public boolean getCsvMode() {
		return csvMode;
	}




	/**
	 * @param csvMode the csvMode to set
	 */
	public void setCsvMode(boolean csvMode) {
		this.csvMode = csvMode;
	}




	/**
	 * @return the bedgraphWriter
	 */
	public PrintWriter getBedgraphWriter() {
		return bedgraphWriter;
	}

	/**
	 * @param bedgraphWriter the bedgraphWriter to set
	 */
	public void setBedgraphWriter(PrintWriter bedgraphWriter) {
		this.bedgraphWriter = bedgraphWriter;
	}


	public void checkWriter()
	{
		try
		{
			if (this.wigWriter == null) throw new Exception();
		}
		catch (Exception e)
		{
			System.err.println("Fatal error, CpgWalkerWindowWigWriter trying to write to a file but hasn't been assigned one yet");
			e.printStackTrace();
			System.exit(1);
		}
	}


	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalker#init()
	 */
	@Override
	protected void init() {
		super.init();
		
		try
		{
			System.err.printf("Initializing CpgWalkerWigWriter: %s\n", this.walkParams.toString());
		}
		catch (Exception e)
		{
		}
	}

	
	
	
	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalker#finishChr()
	 */
	@Override
	public void finishChr() {
		super.finishChr();
		checkWriter();
	}

	
	
	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalker#alertNewChrom()
	 */
	@Override
	protected void alertNewChrom() {
		super.alertNewChrom();
		this.curChrNum = (new ChromFeatures()).chrom_from_public_str(this.curChr);
		System.err.printf("New chrom: %s\n",this.curChr);
		checkWriter();
		wigWriter.printf("variableStep chrom=%s\n", this.curChr);
		lastWindMidpoint = -1;
		lastSpanEnd = -1;
		lastMeth = -1.0;
	}



	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalker#processWindow(java.util.List)
	 */
	@Override
	protected void processWindow(List<Cpg> inWindow) 
	{
		checkWriter();
		super.processWindow(inWindow);

		int windMidpoint = (int)Math.round(((double)(CpgWalker.windEnd(inWindow)+CpgWalker.windStart(inWindow)))/2.0);
		double meth = this.methSummarizer.get(0).getValMean(true);
		int winds = CpgWalker.windStart(inWindow);
		int winde = CpgWalker.windEnd(inWindow);
//		System.err.printf("Processing window CpgWalkerWigWriter: %s\tmeth=%.2f\n", CpgWalker.windStr(inWindow), meth);
		
		// First wig writer
		wigWriter.printf("%d\t%d\n",windMidpoint, (int)Math.round(100.0*meth));

		// Then bedgraph
		if (bedgraphWriter != null)
		{
			int midpointMidpoint = (int)Math.floor(((double)lastWindMidpoint+(double)windMidpoint)/2.0);
			int spanStart = lastSpanEnd + 1;
			int spanEnd = (this.lastWindMidpoint>0) ? midpointMidpoint : -1; // This isn't really correct
			
			// Don't do it if we're on the first one of the chromosome.
			if (this.lastSpanEnd > 0)
			{
				if (this.csvMode)
				{
					bedgraphWriter.printf("%d,%d,%d,%d,%d\n",this.curChrNum, winds, winde, (int)Math.round(100.0*lastMeth), inWindow.size());
				}
				else
				{
					bedgraphWriter.printf("%s\t%d\t%d\t%d\n",this.curChr, spanStart, spanEnd, (int)Math.round(100.0*lastMeth));
				}
			}
//			System.err.printf("\t%s\t%d\t%d\t%d\t%d\n",this.curChr, winds, winde, (int)Math.round(100.0*lastMeth), winde-winds);
			
			this.lastSpanEnd = spanEnd;
			this.lastWindMidpoint = windMidpoint;
			this.lastMeth = meth;			
		}
		
		
	}
	
	
//	private void writeBedgraphRegion(int spanStart, int spanEnd, int meth)
//	{
//		// This doesn's actually work.
//		final int step = 20;
//		final int minBreak = 20;
//		
//		int numRegions = 0;
//		int totalLen = spanEnd-spanStart+1;
//		
//		if (totalLen<minBreak)
//		{
//			this.regionBuffer[numRegions][0] = (spanStart<spanEnd) ? (1+spanStart) : spanStart;
//			this.regionBuffer[numRegions][1] = spanEnd;
//			numRegions++;
//		}
//		else
//		{
//			int numBreaks = (int)Math.ceil( (double)totalLen / (double)step );
//			int breakSize = (int)Math.floor( (double)totalLen / (double)numBreaks );
//			for (int i = 0; i < (numBreaks-1); i++)
//			{
//				int starti = 1 + (spanStart + (i * breakSize));
//				int endi = (i==(numBreaks-2)) ? spanEnd : (starti + breakSize);
//				regionBuffer[numRegions][0] = starti;
//				regionBuffer[numRegions][1] = endi;
//				numRegions++;
//			}
//		}
//		
//		for (int i = 0; i < numRegions; i++)
//		{
//			bedgraphWriter.printf("%s\t%d\t%d\t%d\n",this.curChr, regionBuffer[i][0], regionBuffer[i][1], val);
//		}
//	}
}
