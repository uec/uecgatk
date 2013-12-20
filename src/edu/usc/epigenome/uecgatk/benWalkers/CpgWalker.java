package edu.usc.epigenome.uecgatk.benWalkers;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.logging.Logger;

import org.biojava.bio.seq.StrandedFeature;

import edu.usc.epigenome.genomeLibs.TabularOutput;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbUtils;



import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerParams;
import edu.usc.epigenome.uecgatk.benWalkers.cytosineWalkers.Cpg;

/**
 * @author benb
 *
 * Walks over a sliding window of adjacent Cpgs
 */

/**
 * @author benb
 *
 */
public class CpgWalker implements TabularOutput {

/*
 * 
 * 	
 */
	public CpgWalkerParams walkParams = null;
	protected boolean useSummarizers = true;
	public static final int PROCESS_WINDOW_EVENT = 1;
	protected List<Cpg[]> lastProcessedWindow = null;
	//protected List<String> tables = new ArrayList<String>();
	protected int nTables = 0;
	
	
	// List management
	private LinkedList<Cpg[]> window = new LinkedList<Cpg[]>();
	
	// Some useful summarizers for the window
	protected List<CpgSummarizer> methSummarizer = new ArrayList<CpgSummarizer>();//new CpgMethLevelSummarizer();
	protected List<CpgSummarizer> methSummarizerFw = new ArrayList<CpgSummarizer>();//new CpgMethLevelSummarizerStranded(true);
	protected List<CpgSummarizer> methSummarizerRev = new ArrayList<CpgSummarizer>();//new CpgMethLevelSummarizerStranded(false);
	String curChr = null;
	public String lastChrom = "noChrom";

	// Sends events when it gets a good window
	protected List<ActionListener> listeners = new ArrayList<ActionListener>(10);
	
	/**
	 * 
	 */
//	public CpgWalker(CpgWalkerParams inWalkParams, List<String>inTables) {
//		super();
//		this.walkParams = inWalkParams;
//		this.setnTables(inTables.size());
////		System.out.printf("track name=\"%s\" description=\"%s\" useScore=0 itemRgb=On visibility=4\n",
////				"test", "test");
//
//	}
	
	public CpgWalker(CpgWalkerParams inWalkParams, int inNumTables) {
		super();
		this.walkParams = inWalkParams;
		this.setnTables(inNumTables);

		init();
	}

	/**
	 * Constructor for single table analysis
	 * @param inWalkParams
	 * @param inTables
	 */
	public CpgWalker(CpgWalkerParams inWalkParams) {
		super();
		this.walkParams = inWalkParams;
		
//		List<String> newTables = new ArrayList<String>();
//		String methylTable = (inWalkParams.methylParams!=null) ? inWalkParams.methylParams.methylTablePrefix : "fakeTable";
//		newTables.add(methylTable);
//		this.setTables(newTables);

		this.setnTables(1);
		
		init();

	}

	protected void init()
	{
//		System.out.printf("track name=\"%s\" description=\"%s\" useScore=0 itemRgb=On visibility=4\n",
//		"test", "test");
	
	}
	
	/**
	 * @return the nTables
	 */
	public int getnTables() {
		return nTables;
	}

	/**
	 * @param nTables the nTables to set
	 */
	public void setnTables(int nTables) {
		this.nTables = nTables;
	}

	public int numTables()
	{
		return this.getnTables();
	}




	public boolean useSummarizers() {
		return this.useSummarizers;
	}

	public void useSummarizers(boolean useSummarizers) {
		this.useSummarizers = useSummarizers;
	}

	public void clearWindowListeners()
	{
		listeners = new ArrayList<ActionListener>(10);
	}
	
	public void addWindowListener(ActionListener l)
	{
		listeners.add(l);
	}
	
	/**
	 * @return the curChr
	 */
	public String getCurChr() {
		return curChr;
	}

	/**
	 * @param curChr the curChr to set
	 */
	public void setCurChr(String curChr) {
		this.curChr = curChr;
		resetWindow();
		resetSummarizers();
		this.alertNewChrom();
	}

	protected boolean onNewChrom()
	{
		return !this.lastChrom.equalsIgnoreCase(this.getCurChr());
	}


	/**
	 * Outputs the last domain
	 */
	public void finishChr() {
	}

	
	protected void alertNewChrom()
	{
	}
	
	public void resetWindow()
	{
		window = new LinkedList<Cpg[]>();
	}
	
	public void reset()
	{
	}
	
	
	public void resetSummarizers()
	{
		//System.err.printf("resetting summarizers for %d tables (useSummarizers=%s)\n",this.numTables(),useSummarizers);
		if (useSummarizers)
		{
			methSummarizer = new ArrayList<CpgSummarizer>();
			methSummarizerFw = new ArrayList<CpgSummarizer>();
			methSummarizerRev = new ArrayList<CpgSummarizer>();
			for (int i = 0; i < this.numTables(); i++)
			{
				methSummarizer.add(new CpgMethLevelSummarizer());
				methSummarizerFw.add(new CpgMethLevelSummarizerStranded(true));
				methSummarizerRev.add(new CpgMethLevelSummarizerStranded(false));
			}
		}
	}
	
//	public double methMean(boolean useSpaceWeightedMean)
//	{
//		return this.methMean(0, useSpaceWeightedMean);
//	}
//	
//	public double methMean(int tableIndex, boolean useSpaceWeightedMean)
//	{
//		return this.getSummarizer(tableIndex).getValMean(useSpaceWeightedMean);
//	}
//	
//	public CpgSummarizer getSummarizer()
//	{
//		return getSummarizer(0);
//	}
//	
//	public CpgSummarizer getSummarizer(int tableIndex)
//	{
//		return this.methSummarizer.get(tableIndex);
//	}
	
	public List<Cpg[]> getLastProcessedWindow() {
		return lastProcessedWindow;
	}

	/**
	 * This uses a fixed window size and is much faster
	 * @param cpg Must be streamed in serially.
	 */
	protected void streamCpgFixedWind(Cpg[] cpg)
	{
		int newPos = cpg[0].chromPos;
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).fine(
				String.format("Found Cpg: %d\n", newPos));

		// Add this Cpg to the head of the queue
		window.add(cpg);
		if (useSummarizers)
		{
			for (int t = 0; t < this.numTables(); t++)
			{
				methSummarizer.get(t).streamCpg(cpg[t]);
				methSummarizerFw.get(t).streamCpg(cpg[t]);
				methSummarizerRev.get(t).streamCpg(cpg[t]);
			}
		}


		// Remove cpgs from the tail
		boolean done = false;
		Cpg[] endCpg, secondCpg;
//		while (!done && ((endCpg = window.peek()) != null))
//		{
//			if ((newPos - endCpg[0].chromPos) < this.walkParams.maxScanningWindSize)
//			{
//				done = true;
//			}
		while (!done && (window.size()>1))
		{
			endCpg = window.get(0);
			secondCpg = window.get(1);
			if ((newPos - secondCpg[0].chromPos) < this.walkParams.minScanningWindSize)
			{
				done = true;
			}
			else
			{
				window.remove();
				if (useSummarizers)
				{
					for (int t = 0; t < this.numTables(); t++)
					{
						methSummarizer.get(t).removeCpg(endCpg[t]);
						methSummarizerFw.get(t).removeCpg(endCpg[t]);
						methSummarizerRev.get(t).removeCpg(endCpg[t]);
					}
				}
			}
		}
		
//		System.err.printf("\tFixed wind checking %s (id=%d)\n",this.windStr(),this.hashCode());
		
		// And process the window
		if ((window.size()>=walkParams.minScanningWindCpgs) &&
				((CpgWalker.windEnd(window,false)-CpgWalker.windStart(window,false)) >= walkParams.minScanningWindSize))
		{
//			System.err.println("\t\t Sufficient Cpgs: " + window.size());
//			double mean = this.methSummarizer.getValMean(true);
//			System.err.println("\t\t Mean meth=" + mean);
//			if (mean < 0.7)
//			System.out.println(MethylDbUtils.bedLine("chr11", windStart(), windEnd(), ".", mean));
			this.processWindow(this.window, true);
		}
		else
		{
//			System.err.println("\t\t Not enough Cpgs");
		}
	}
		
	/**
	 * This uses a fixed step, and can have variable window sizing to accommodate
	 * widely varying CG density.
	 * @param cpg Must be streamed in serially.
	 */
	protected void streamCpgVariableWind(Cpg[] cpg)
	{
		int newPos = cpg[0].chromPos;
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).fine(
				String.format("Variable wind found Cpg: %d\n", newPos));

		// Add this Cpg to the head of the queue
		window.add(cpg);


		// Remove cpgs from the tail
		boolean done = false;
		Cpg[] endCpg;
		while (!done && ((endCpg = window.peek()) != null))
		{
			if ((newPos - endCpg[0].chromPos) < this.walkParams.maxScanningWindSize)
			{
				done = true;
			}
			else
			{
				window.remove();
			}
		}

		//System.err.println("\tChecking " + this.windStr());
		
		// And process the window
		//System.err.println(this.windStr());
		if (window.size()>=walkParams.minScanningWindCpgs)
		{
			
			// First we set our summarizers.
			// Take the last minCpgs as the sub-window.  Do it as an iterator since
			// a linked list might be quicker iterating than using "get(i)".
			Iterator<Cpg[]> backIt = window.descendingIterator();
			int i = 0;
			int lastPos = 0;
			if (useSummarizers) this.resetSummarizers();
			boolean minSizeReached = false;
			while ((i<window.size()) && !(minSizeReached && (i>walkParams.minScanningWindCpgs)))
			{
				// ***** REMOVE THIS
				if (!backIt.hasNext()) System.err.println("Why did we run out of window elements?!"); // REMOVE
				// ***** REMOVE THIS

				Cpg[] backCpg = backIt.next();
				//System.err.println("\t(i=" + i + ") cpg=" + backCpg.chromPos);
				lastPos = backCpg[0].chromPos;
				
				if (useSummarizers)
				{
					for (int t = 0; t < this.numTables(); t++)
					{
						methSummarizer.get(t).streamCpg(backCpg[t]);
						methSummarizerFw.get(t).streamCpg(backCpg[t]);
						methSummarizerRev.get(t).streamCpg(backCpg[t]);
					}
				}
				
				int windLen = newPos-lastPos;
				minSizeReached = (windLen>=this.walkParams.minScanningWindSize);
				//System.err.printf("Checking size %d (minSizeReached=%s)\n",windLen,minSizeReached);

				i++;
			}


			// Process this window
	//		System.err.printf("Checking size %d (minSizeReached=%s)\n",cpg.chromPos-lastPos+1,minSizeReached);
			if (minSizeReached && (i>walkParams.minScanningWindCpgs))
			{
				this.processWindow(this.window.subList(this.window.size()-(i-1), this.window.size()), true);
			}
		}
	}

	/**
	 * @param cpg
	 */
	public void streamCpg(Cpg[] cpg)
	{
		//System.err.println("Streaming CpG: " + cpg.toStringExpanded());
		//System.err.println("\tuseFixedStep=" + !this.walkParams.useVariableWindow);
		
		// Sanity checks
		if (cpg.length != this.numTables())
		{
			System.err.printf("Trying to stream CpGIterator with %d samples to CpGWalker with %d samples, exiting.\n",cpg.length,this.numTables());
			System.exit(1);
		}
		
		if (this.curChr == null)
		{
			System.err.printf("Trying to stream CpGIterator with without first calling walker.setCurChrom. Quitting\n");
			System.exit(1);
		}
		
		// Automatically set weights if possible
		if (Double.isNaN(cpg[0].getCpgWeight()))
		{
			int lastPos = windEnd();
			if (lastPos>=1)
			{
				int dist = cpg[0].chromPos - lastPos;
				int windSize = this.window.size();
				
				
				if (dist == 1)
				{
					StrandedFeature.Strand thisStrand = cpg[0].getStrand();
					if ( (thisStrand == StrandedFeature.NEGATIVE) &&  (windSize >= 1) ) 
					{
						Cpg lastCpg = window.get(windSize-1)[0];
						if (lastCpg.getStrand() == StrandedFeature.POSITIVE)
						{
							// This means we are the opposite strand of the previous one.
							if (windSize == 1)
							{
								dist = (int)lastCpg.getCpgWeight(); // cpg[0].chromPos; // This doesn't work because it can be the beginning of the chrom . 
							}
							else
							{
								Cpg secondLast = window.get(windSize-2)[0];
								dist = lastCpg.chromPos - secondLast.chromPos;
							}
						}
					}
				}
				
				//System.err.printf("Setting weight for cpg (cur=%d, last=%d) = %d\n",cpg[0].chromPos,lastPos,dist);
				for (int i = 0; i < cpg.length; i++)
				{
					cpg[i].setCpgWeight((double)dist);
				}
			}
		}
		
		
		
		if (this.walkParams.useVariableWindow)
		{
			this.streamCpgVariableWind(cpg);
		}
		else
		{
			this.streamCpgFixedWind(cpg);
		}
	}

	/**
	 * Include for backwards compatibility with single-table analyses
	 * @param cpg
	 */
	public void streamCpg(Cpg cpg)
	{
		Cpg[] multi = new Cpg[1];
		multi[0]  = cpg;
		streamCpg(multi);
	}
	
	public String lastProcessedWindStr(boolean longVers)
	{
		return windStr(this.getLastProcessedWindow(),longVers);
	}

	public String windStr()
	{
		return windStrMulti(this.window);
	}
	
	public static String windStrMulti(List<Cpg[]> inWind)
	{
		return windStr(inWind, false);
	}

	public static String windStr(List<Cpg> inWind)
	{
		return windStr(singleWindowToMulti(inWind), false);
	}

	public String windStr(boolean longVers)
	{
		return windStr(this.window, longVers);
	}

	
	public static String windStr(List<Cpg[]> inWind, boolean longVers)
	{
		StringBuffer sb = new StringBuffer(10000);
		
		int e = windEnd(inWind,true);
		int s = windStart(inWind,true);
		sb.append(String.format("Wind %d-%d (%d bp) has %d elements: ",
				s, e, e-s+1, inWind.size()));
		if (longVers)
		{
			for (int i = 0; i < inWind.size(); i++)
			{
				if (i>0) sb.append(", ");
				Cpg[] cpg = inWind.get(i);
				sb.append(cpg[0].chromPos);
			}
		}
		
		
		return sb.toString();
	}

	public int windLen()
	{
		return windLen(this.window,true);
	}
	
	public int windStart()
	{
		return windStart(window,true);
	}

	public int windEnd()
	{
		return windEnd(window,true);
	}
	
	public static int windLen(List<Cpg[]> inWind, boolean fake)
	{
		return windEnd(inWind,fake) - windStart(inWind,fake) + 1;
	}
	
	public static int windStart(List<Cpg[]> inWind, boolean fake)
	{
		if (inWind.size()==0) return -1;
		return inWind.get(0)[0].chromPos;
	}

	public static int windEnd(List<Cpg[]> inWind, boolean fake)
	{
		if (inWind.size()==0) return -1;
		return inWind.get(inWind.size()-1)[0].chromPos;
	}
	
	// FOR BACKWARDS COMPATIBILITY	
	public static int windLen(List<Cpg> inWind)
	{
		return windEnd(inWind) - windStart(inWind) + 1;
	}
	
	public static int windStart(List<Cpg> inWind)
	{
		if (inWind.size()==0) return -1;
		return inWind.get(0).chromPos;
	}

	public static int windEnd(List<Cpg> inWind)
	{
		if (inWind.size()==0) return -1;
		return inWind.get(inWind.size()-1).chromPos;
	}

	
	/***
	 ******** OVERRIDE ONE OR BOTH OF THESE *********
	 * @param inWindow
	 */
	protected static List<Cpg> multiWindowToSingle(List<Cpg[]> inMulti)
	{
		// Call into the single table version to provide backwards compatibility
		List<Cpg> singleSampleWindow = new ArrayList<Cpg>(inMulti.size());
		for (int i = 0; i < inMulti.size(); i++)
		{
			singleSampleWindow.add(inMulti.get(i)[0]);
		}
		return singleSampleWindow;
	}
	
	protected static List<Cpg[]> singleWindowToMulti(List<Cpg> inSingle)
	{
		// Call into the single table version to provide backwards compatibility
		List<Cpg[]> multi = new ArrayList<Cpg[]>(inSingle.size());
		for (int i = 0; i < inSingle.size(); i++)
		{
			Cpg[] multiCpg = new Cpg[1];
			multiCpg[0] = inSingle.get(i);
			multi.add(multiCpg);
		}
		return multi;
	}


	protected void processWindow(List<Cpg[]> inWindow, boolean fake) // The fake is just to give it a different signature from the List<Cpg> version
	{
		// Call into the single table version to provide backwards compatibility
		List<Cpg> singleSampleWindow = multiWindowToSingle(inWindow);
		//System.err.println("ProcessWindow, making single sample version: " + this.windStr());
		processWindow(singleSampleWindow);
		singleSampleWindow.clear();
		
		
		// Now bookkeeping
		lastProcessedWindow = inWindow;
		
		// Basic version does nothing except alert listeners
		ActionEvent e = new ActionEvent(this, PROCESS_WINDOW_EVENT, null);
		for (ActionListener l : listeners)
		{
			l.actionPerformed(e);
		}

		this.lastChrom = this.curChr;
	}

	/**
	 * This was the version prior to making the class multi-sample.  It is called for backwards compatibility with old sub-classes
	 * @param inWindow
	 */
	protected void processWindow(List<Cpg> inWindow)
	{
	}
	
	
	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.TabularOutput#headerStr()
	 */
	public String headerStr() throws Exception {
		return null;
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.TabularOutput#toCsvStr()
	 */
	public String toCsvStr() throws Exception {
		return null;
	}

	
}
