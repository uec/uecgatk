package edu.usc.epigenome.uecgatk.benWalkers.cytosineWalkers;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.biojava.bio.seq.StrandedFeature;

import edu.usc.epigenome.uecgatk.benWalkers.CytosineContextCounter;





public class Cpg implements Comparable, Cloneable {

	private static final double DEFAULT_FAIL_MIN_FRAC_A = 0.2;
	private static final int DEFAULT_FAIL_MIN_ABS_A = 2;
	
	// Object vars
	public int chromPos = 0;
	public boolean negStrand = false;
	public short totalReads = 0;
	public short cReads = 0;
	public short cReadsNonconversionFilt = 0;
	public short tReads = 0;
	public short agReads = 0;
	public short totalReadsOpposite = 0;
	public short aReadsOpposite = 0;
	
	protected CytosineContextCounter contextCounter = new CytosineContextCounter();
	
//	private short nextBaseTotalReads = 0;
//	private short nextBaseGreads = 0;
//	
//	private short prevBaseTotalReads = 0;
//	private short prevBaseGreads = 0;
//	
	protected String nextBasesRefUpperCase = "0";
	protected String prevBasesRefUpperCase = "0";
	
	// This is for weighted averages
	protected double cpgWeight = Double.NaN;
	
	// Read info
	protected Map<Integer,CpgRead> readMap = new TreeMap<Integer,CpgRead>();



	/**
	 * 
	 */
	public Cpg() {
		super();
	}

	/**
	 * @param chromPos
	 * @param negStrand
	 * @param totalReads
	 * @param cReads
	 * @param cReadsNonconversionFilt
	 * @param tReads
	 * @param agReads
	 * @param totalReadsOpposite
	 * @param aReadsOpposite
	 */
	public Cpg(int chromPos, boolean negStrand) {
		super();
		this.chromPos = chromPos;
		this.negStrand = negStrand;
	}
	
	
	public Cpg(int chromPos, boolean negStrand, short totalReads, short cReads,
			short cReadsNonconversionFilt, short tReads, short agReads,
			short totalReadsOpposite, short aReadsOpposite, int cpgWeight,
			CytosineContextCounter inCounter, char inPrevBaseRef, char inNextBaseRef)
	 {
		super();
		this.init(chromPos, negStrand, totalReads, cReads, cReadsNonconversionFilt, tReads,
				agReads, totalReadsOpposite, aReadsOpposite, cpgWeight, inCounter,
				Character.toString(inNextBaseRef), Character.toString(inPrevBaseRef));
	 }
	
	public Cpg(int chromPos, boolean negStrand, short totalReads, short cReads,
			short cReadsNonconversionFilt, short tReads, short agReads,
			short totalReadsOpposite, short aReadsOpposite, int cpgWeight,
			short inNextBaseGreads, short inNextBaseTotalReads, char inNextBaseRef)
	 {
		super();
		this.initBackwardsCompatible(chromPos, negStrand, totalReads, cReads, cReadsNonconversionFilt, tReads,
				agReads, totalReadsOpposite, aReadsOpposite, cpgWeight, inNextBaseGreads,
				inNextBaseTotalReads, (short)0, (short)0, Character.toString(inNextBaseRef), "0");
	 }
	
	public Cpg(int chromPos, boolean negStrand, short totalReads, short cReads,
			short cReadsNonconversionFilt, short tReads, short agReads,
			short totalReadsOpposite, short aReadsOpposite, int cpgWeight,
			short inNextBaseGreads, short inNextBaseTotalReads, short inPrevBaseGreads,
			short inPrevBaseTotalReads,String inPrevBasesRef,String inNextBasesRef)
	 {
		super();
		this.initBackwardsCompatible(chromPos, negStrand, totalReads, cReads, cReadsNonconversionFilt, tReads,
				agReads, totalReadsOpposite, aReadsOpposite, cpgWeight, inNextBaseGreads,
				inNextBaseTotalReads, inPrevBaseGreads, inPrevBaseTotalReads, inNextBasesRef, inPrevBasesRef);
	}
	

	
	protected void init(int chromPos, boolean negStrand, short totalReads, short cReads,
			short cReadsNonconversionFilt, short tReads, short agReads,
			short totalReadsOpposite, short aReadsOpposite, int cpgWeight,
			CytosineContextCounter inCounter, String inNextBasesRef, String inPrevBasesRef)
	{
		this.initCommon(chromPos, negStrand, totalReads, cReads, cReadsNonconversionFilt, tReads, agReads, totalReadsOpposite, aReadsOpposite, cpgWeight);
		
		this.setContextCounter(inCounter);
		this.nextBasesRefUpperCase = inNextBasesRef.toUpperCase();
		this.prevBasesRefUpperCase = inPrevBasesRef.toUpperCase();
	}

	protected void initBackwardsCompatible(int chromPos, boolean negStrand, short totalReads, short cReads,
			short cReadsNonconversionFilt, short tReads, short agReads,
			short totalReadsOpposite, short aReadsOpposite, int cpgWeight,
			short inNextBaseGreads, short inNextBaseTotalReads, short inPrevBaseGreads, 
			short inPrevBaseTotalReads, String inNextBasesRef, String inPrevBasesRef)
	{
		this.initCommon(chromPos, negStrand, totalReads, cReads, cReadsNonconversionFilt, tReads, agReads, totalReadsOpposite, aReadsOpposite, cpgWeight);
		
		this.contextCounter.setNextCountsCompatibilty(inNextBaseTotalReads, inNextBaseGreads);
		this.nextBasesRefUpperCase = inNextBasesRef.toUpperCase();
		
		this.contextCounter.setPrevCountsCompatibilty(inPrevBaseTotalReads, inPrevBaseGreads);
		this.prevBasesRefUpperCase = inPrevBasesRef.toUpperCase();


	}

	private void initCommon(int chromPos, boolean negStrand, short totalReads, short cReads,
			short cReadsNonconversionFilt, short tReads, short agReads,
			short totalReadsOpposite, short aReadsOpposite, int cpgWeight)
	{
		this.chromPos = chromPos;
		this.negStrand = negStrand;
		this.totalReads = totalReads;
		this.cReads = cReads;
		this.cReadsNonconversionFilt = cReadsNonconversionFilt;
		this.tReads = tReads;
		this.agReads = agReads;
		this.totalReadsOpposite = totalReadsOpposite;
		this.aReadsOpposite = aReadsOpposite;		
		
		// CpG weight has some issues.  For instance at large gaps you get huge values.  Remove these
		if (cpgWeight < 0)
		{
			System.err.printf("Got negative cpgWeight: pos=%d, weight=%d\n",chromPos,cpgWeight);
			cpgWeight = 100;
		}
		else if (cpgWeight > 5000)
		{
			System.err.printf("Got extra-large cpgWeight: pos=%d, weight=%d\n",chromPos,cpgWeight);
			cpgWeight = 100;
		}
		
		// It's actually divided by two since we're looking at both strands.
		this.cpgWeight = (double)cpgWeight/2.0;	
	}
	
	
	
	
	/**
	 * @return the contextCounter
	 */
	public CytosineContextCounter getContextCounter() {
		return contextCounter;
	}

	/**
	 * @param contextCounter the contextCounter to set
	 */
	public void setContextCounter(CytosineContextCounter contextCounter) {
		this.contextCounter = contextCounter;
	}

	/*** Read information ***/
	public void addRead(CpgRead read)
	{
		readMap.put(read.readId, read);
		
		this.agReads += read.agRead;
		this.tReads += read.tRead;
		this.cReads += read.cRead;
		this.totalReads += (read.agRead + read.tRead + read.cRead);
		this.cReadsNonconversionFilt += read.cReadNonconversionFilt;

		this.contextCounter.addContext(read.getContext());
		
		// Debugging
//		System.err.printf("Adding read with context \"%s\". New contextCounter:\n%s\n",read.getContext(),this.contextCounter.toString());
	}
	
	public Map<Integer,CpgRead> getReads()
	{
		return readMap;
	}
	
	
	/*** Overridden Comparable methods
	 */

	/* (non-Javadoc)
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	public int compareTo(Object o) {

		Cpg other = (Cpg)o;
		
		return (this.chromPos - other.chromPos);
	}


	/* (non-Javadoc)
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object obj) {
		return 	(this.compareTo(obj) == 0);
	}
	
	public static Cpg compositeCpg(List<Cpg> inCpgs) throws CloneNotSupportedException
	{
		Cpg out = null;
		
		for (Cpg inCpg : inCpgs)
		{
			if (out == null)
			{
				// First one
				out = (Cpg)inCpg.clone();
			}
			else
			{
				out.totalReads += inCpg.totalReads;
				out.cReads += inCpg.cReads;
				out.cReadsNonconversionFilt += inCpg.cReadsNonconversionFilt;
				out.tReads += inCpg.tReads;
				out.agReads += inCpg.agReads;
				out.totalReadsOpposite += inCpg.totalReadsOpposite;
				out.aReadsOpposite += inCpg.aReadsOpposite;
			}
		
		}
		
		return out;
	}
	

	public double fracNextBaseG()
	{
		return (double)this.contextCounter.getGCount(1) / (double)this.contextCounter.getTotalCount(1);
	}

	public double fracPrevBaseG()
	{
		return (double)this.contextCounter.getGCount(-1) / (double)this.contextCounter.getTotalCount(-1);
	}

	public double fracOppositeA()
	{
		return (double)this.aReadsOpposite / ((double)this.aReadsOpposite + (double)this.totalReadsOpposite);
	}
	
	public boolean passesOppositeAFilterDefault()
	{
		return passesOppositeAFilter(DEFAULT_FAIL_MIN_FRAC_A, DEFAULT_FAIL_MIN_ABS_A);
	}
	
	public boolean passesOppositeAFilter(double failMinFracA, int failMinAbsA)
	{
		double fracOppA = this.fracOppositeA();
		boolean failFrac = (!Double.isNaN(fracOppA)) && (fracOppA>=failMinFracA);
		boolean failAbs = (this.aReadsOpposite>=failMinAbsA);
		return !(failFrac && failAbs);
	}

	public double fracReadsCorT()
	{
		return (double)totalReadsCorT(false) / ((double)totalReadsCorT(false)+(double)agReads);
	}
	
	public double fracMeth(boolean useNonconvFilt)
	{
		int cs = (useNonconvFilt) ? this.cReads : (this.cReads + this.cReadsNonconversionFilt); 
		return (double)cs / ((double)cs + (double)this.tReads);
	}
	
	public String variableStepWigLine(boolean useNonconvFilter)
	{
		return String.format("%d\t%f",chromPos, this.fracMeth(useNonconvFilter));
	}

	public static void outputCpgsToFile(PrintWriter pw, Map<Integer,Cpg> cpgMap, String prefix, String sampleName, String chr, int minCphCoverage, double minCphMethFrac)
	throws IOException
	{
		outputCpgsToFile(false, pw, cpgMap, prefix, sampleName, chr, minCphCoverage, minCphMethFrac);
	}
	
	public static void outputCpgsToFile(boolean readFormat, PrintWriter pw, Map<Integer,Cpg> cpgMap, String prefix, String sampleName, String chr, int minCphCoverage, double minCphMethFrac)
	throws IOException
	{
		//System.err.println("About to write " + cpgMap.size() + " Cytosines to file");
		Iterator<Cpg> cpgIt = cpgMap.values().iterator();
		CPG: while (cpgIt.hasNext())
		{
			Cpg cpg = cpgIt.next();
			if ((cpg.getNextBaseRef() != 'G'))// && (cpg.getNextBaseRef() != '0'))
			{
				if ((cpg.totalReads < minCphCoverage) || (cpg.fracMeth(true)<minCphMethFrac))
				{
					continue CPG;
				}
			}
			
			String line = (readFormat) ? cpg.toReadFormatString() : cpg.toString();
			pw.print(line);
			if (!readFormat) pw.println();
		}
	}
	
	// The default version now demands the reads to be matching. Be careful this doesn't screw anything up.
	public String context(int numPrev, int numPost)
	{
		return this.context(numPrev, numPost,0.899,true);
	}
	

	
	//*************************/
	//**** CRITICAL LOGIC *****/
	//*************************/
	public String context(int numPrev, int numPost, double minFracReadsMatching, boolean bisulfiteMatching)
	{
		int numPrevAvail = this.contextCounter.numPrevBases();
		if (numPrev > numPrevAvail)
		{
//			System.err.printf("Cpg::context() Requesting %d prevBases when only %d exist\n", numPrev, numPrevAvail);
			numPrev = numPrevAvail;
		}
		
		int numPostAvail = this.contextCounter.numNextBases();
		if (numPost > numPostAvail)
		{
//			System.err.printf("Cpg::context() Requesting %d nextBases when only %d exist\n", numPost, numPostAvail);
			numPost = numPostAvail;
		}
		
		StringBuffer sb = new StringBuffer(numPrev+1+numPost);
		for (int i = -numPrev; i <= numPost; i++)
		{
			if (i == 0)
			{
				sb.append("C");
			}
			else
			{
				boolean done = false;
				char comparisonBase = 'N';
				if ((i==1) || (i==-1))
				{
					comparisonBase = Character.toUpperCase((i==1) ? this.getNextBaseRef() : this.getPrevBaseRef());
					done = (comparisonBase != 'N');
				}

				if (!done)
				{
					comparisonBase = GatkBaseUtils.CharFromBaseByte(this.contextCounter.majorityBase(i));
					done = true;
				}

				byte comparisonByte = GatkBaseUtils.BaseByteFromChar(comparisonBase);
				double fracTargetRef = this.contextCounter.getMatchingFrac(i, comparisonByte, bisulfiteMatching);
//				System.err.printf("\tcomparisonBase=%c\tcomparisonByte=%d\tfracMatching=%.2f (coord=%d, relPos=%s)\n", 
//						comparisonBase, comparisonByte, fracTargetRef, this.chromPos, i);
				if (Double.isNaN(fracTargetRef) || (fracTargetRef>=minFracReadsMatching))
				{
					sb.append(comparisonBase);
				}
				else
				{
					sb.append("N");
				}
			}
		}
		

//		// Prev
//		String out = "";
//		double fracPrevTarget=Double.NaN, fracNextTarget=Double.NaN;
//		if (numPrev > 0)
//		{
//			fracPrevTarget = this.fracPrevBaseG();
//			if (Double.isNaN(fracPrevTarget)) fracPrevTarget = 1.0;
//			else if (!BaseUtils.basesAreEqual((byte)this.getPrevBaseRef(),BaseUtils.G)) fracPrevTarget = 1 - fracPrevTarget;
//			char base = (fracPrevTarget>=minFracReadsMatching) ? this.getPrevBaseRef() : 'N';
//			out += base;
//		}
//
//		out += "C";
//
//		// Next
//		if (numPost > 0)
//		{
//			 fracNextTarget = this.fracNextBaseG();
//			if (Double.isNaN(fracNextTarget)) fracNextTarget = 1.0;
//			else if (!BaseUtils.basesAreEqual((byte)this.getNextBaseRef(),BaseUtils.G)) fracNextTarget = 1 - fracNextTarget;
//			char base = (fracNextTarget>=minFracReadsMatching) ? this.getNextBaseRef() : 'N';
//			out += base;
//		}
//
////		System.err.printf("%s", this.toReadFormatString());
////		System.err.printf("\tCONTEXT(%d): Context %s (frac matching prevRef %c = %.2f, frac matching nextRef %c = %.2f)\n",
////				this.chromPos, out, this.getPrevBaseRef(), fracPrevTarget, this.getNextBaseRef(), fracNextTarget);
		
		return sb.toString();
	}
	
	public String context()
	{
		return this.context(1,1);
	}

	public String context(double minFracReadsMatching)
	{
		return this.context(1,1,minFracReadsMatching,true);
	}

	public boolean isCph(boolean onlyUseRef)
	{
		return isCph(onlyUseRef, 0.5);
	}
	
	public boolean isCph(boolean onlyUseRef, double maxNextBaseGfrac)
	{
		boolean iscph = false;
		
		iscph = (this.getNextBaseRef() != 'G'); 

		if (!onlyUseRef)
		{
			if (this.contextCounter.getTotalCount(1) > 0)
			{
				iscph = (this.fracNextBaseG() < maxNextBaseGfrac);
			}
		}
		
		return iscph;
	}

	public static PrintWriter outputChromToFile(Map<Integer,Cpg> cpgMap, String prefix, String sampleName, String chr, int minCphCoverage, double minCphMethFrac)
	throws IOException
	{
		return outputChromToFile(false, cpgMap, prefix, sampleName, chr,  minCphCoverage, minCphMethFrac);
	}
		
	public static void outputChromToFile(Map<Integer,Cpg> cpgMap, String prefix, String sampleName, String chr)
	throws IOException
	{
		outputChromToFile(cpgMap, prefix, sampleName, chr, 0, 0.0);
	}
	
	public static PrintWriter outputChromToFile(boolean readFormat, Map<Integer,Cpg> cpgMap, String prefix, String sampleName, String chr, 
			int minCphCoverage, double minCphMethFrac)
	throws IOException
	{
		
		String fn = prefix + sampleName + "_" + chr + ".txt";
		PrintWriter writer = new PrintWriter(new File(fn));
		
		outputCpgsToFile(readFormat, writer, cpgMap, prefix, sampleName, chr, minCphCoverage, minCphMethFrac);
		
		//writer.close();
		return writer;
	}

	@Override
	public String toString() {

		if (true)
		{
			return this.toStringExpanded();
		}
		else
		{
		return String.format("%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%d", 
				chromPos,
				(negStrand) ? '-' : '+',
				totalReads,
				cReads,
				cReadsNonconversionFilt,
				tReads,
				agReads,
				totalReadsOpposite,
				aReadsOpposite
				);
		}
	}
	
	public String toStringExpanded() 
	{

		return String.format("%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%.2f\t%.2f\t%.2f", 
				chromPos,
				(negStrand) ? '-' : '+',
				totalReads,
				cReads,
				cReadsNonconversionFilt,
				tReads,
				agReads,
				totalReadsOpposite,
				aReadsOpposite,
				this.contextCounter.getGCount(1),
				this.contextCounter.getTotalCount(1),
				nextBasesRefUpperCase,
				prevBasesRefUpperCase,
				100*this.fracMeth(true),
				100*this.fracNextBaseG(),
				this.cpgWeight
				);
	}
	
	public String toReadFormatString() {
		
		StringBuffer sb = new StringBuffer(100*this.readMap.size());
		Iterator<CpgRead> readIt = this.readMap.values().iterator();
		while (readIt.hasNext())
		{
			CpgRead read = readIt.next();
			

			sb.append(String.format("%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%c\t%s\t%c\t%.2f\n", 
					chromPos,
					((negStrand) ? '-' : '+'),
					read.readId,
					
					read.cRead,
					read.cReadNonconversionFilt,
					read.tRead,
					read.agRead,
					
					totalReadsOpposite,
					aReadsOpposite,
					
					read.getContext().prevBaseGread(),
					read.getContext().contextPrevBase(),

					read.getContext().nextBaseGread(),
					read.getContext().contextNextBase(),
					
					this.cpgWeight
			));
		}

//		  `chromPos` INT UNSIGNED NOT NULL,
//		  `strand` enum('+','-') NOT NULL,
//		  `readId` INT UNSIGNED	NOT NULL,
//
//		  `cRead` SMALLINT UNSIGNED NOT	NULL,
//		  `cReadNonconversionFilt` SMALLINT UNSIGNED NOT NULL,
//		  `tRead` SMALLINT UNSIGNED NOT NULL,
//		  `agRead` SMALLINT UNSIGNED NOT NULL,
//
//		  `totalReadsOpposite` SMALLINT UNSIGNED NOT NULL,
//		  `aReadsOpposite` SMALLINT UNSIGNED NOT NULL,
//
//		  `nextBaseGread` SMALLINT UNSIGNED NULL default '0',
//		  `cpgWeight` SMALLINT UNSIGNED NOT NULL,

		return sb.toString();
	}
	
	
	@Override
	protected Object clone() throws CloneNotSupportedException {
		Cpg out = (Cpg)super.clone();
		
		//System.err.printf("Cloning cpg, old tReads=%d, new tReads=%d\n",  this.tReads, out.tReads);
		return out;
	}

	public double getCpgWeight() {
		return cpgWeight;
	}

	public void setCpgWeight(double cpgWeight) {
		this.cpgWeight = cpgWeight;
	}

	public String getNextBasesRef() {
		return nextBasesRefUpperCase;
	}

	public void setNextBasesRef(String inNextBasesRef) {
		this.nextBasesRefUpperCase = inNextBasesRef.toUpperCase();
	}

	/**
	 * @return the prevBaseRefUpperCase
	 */
	public String getPrevBasesRef() {
		return prevBasesRefUpperCase;
	}

	/**
	 * @param prevBasesRefUpperCase the prevBaseRefUpperCase to set
	 */
	public void setPrevBasesRef(String inPrevBasesRef) {
		this.prevBasesRefUpperCase = inPrevBasesRef.toUpperCase();
	}

	public char getNextBaseRef() {
		return nextBasesRefUpperCase.charAt(0);
	}

	public void setNextBaseRef(char inNextBaseRef) {
		this.setNextBasesRef(Character.toString(inNextBaseRef));
	}

	/**
	 * @return the prevBaseRefUpperCase
	 */
	public char getPrevBaseRef() {
		return prevBasesRefUpperCase.charAt(0);
	}

	/**
	 * @param prevBasesRefUpperCase the prevBaseRefUpperCase to set
	 */
	public void setPrevBaseRef(char inPrevBaseRef) {
		this.setPrevBasesRef(Character.toString(inPrevBaseRef));
	}
	
	public StrandedFeature.Strand getStrand()
	{
		return (this.negStrand) ? StrandedFeature.NEGATIVE : StrandedFeature.POSITIVE;
	}
	
	public String getStrandStr()
	{
		return (this.negStrand) ? "-" : "+";
	}
	
	public short totalReadsC(boolean useNonconvFilt)
	{
		int cs = (useNonconvFilt) ? this.cReads : (this.cReads + this.cReadsNonconversionFilt); 
		return (short)cs;
	}
	
	public short totalReadsCorT(boolean useNonconvFilt)
	{
		short cs = totalReadsC(useNonconvFilt);
		return (short)(cs + this.tReads);
	}

	public Cpg downsample(double downsamplingFactor) {
		// TODO Auto-generated method stub
		
		
		Cpg out = null;
		try
		{
			out = (Cpg)this.clone();
			
			// For positive strand, we have all the counts so we can sample every flavor separately
			out.cReads = (short)MatUtils.randomBinomialGeneratedCount(this.cReads, downsamplingFactor);
			out.cReadsNonconversionFilt = (short)MatUtils.randomBinomialGeneratedCount(this.cReadsNonconversionFilt, downsamplingFactor);
			out.tReads = (short)MatUtils.randomBinomialGeneratedCount(this.tReads, downsamplingFactor);
			out.agReads = (short)MatUtils.randomBinomialGeneratedCount(this.agReads, downsamplingFactor);
			out.totalReads = (short)(out.cReads + out.cReadsNonconversionFilt + out.tReads + out.agReads);
			
			// For reverse strand and next base, we only have counts for one variant, so we have to resample the total
			// and then sample the single variant based on the observed fraction
			out.totalReadsOpposite = (short)MatUtils.randomBinomialGeneratedCount(this.totalReadsOpposite, downsamplingFactor);
			out.aReadsOpposite = (short)MatUtils.randomBinomialGeneratedCount(out.totalReadsOpposite, this.fracOppositeA());
			out.contextCounter.setNextCountsCompatibilty(MatUtils.randomBinomialGeneratedCount(this.contextCounter.getTotalCount(1), downsamplingFactor),
					MatUtils.randomBinomialGeneratedCount(this.contextCounter.getTotalCount(1), this.fracNextBaseG()));
		}
		catch (Exception e)
		{
			System.err.printf("Problem downsampling: %s\n",e.toString());
			e.printStackTrace();
			System.exit(1);
		}

		// - - - REMOVE --- REMOVE --- REMOVE !!!
//		System.out.printf("%d,%.2f,%.2f,%.2f,%.2f\n", this.totalReads, this.fracMeth(true), out.fracMeth(true), this.fracMeth(false), out.fracMeth(false));

		
//		System.err.printf("Before: %s\n",this.toStringExpanded());
//		System.err.printf("After: %s\n\n",out.toStringExpanded());
			
		
		
//		System.err.println("Cpg::downsample not yet supported");
//		System.exit(1);
		
		return out;
	}
	
	
	
	
}
