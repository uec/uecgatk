package edu.usc.epigenome.uecgatk.benWalkers.cytosineReadWalkers;

import java.util.Collections;
import java.util.Set;
import java.util.TreeMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

import org.apache.commons.math.fraction.Fraction;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.utils.collections.Pair;

import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.uecgatk.FractionNonidentical;
import edu.usc.epigenome.uecgatk.benWalkers.ReadWalkerToBisulfiteCytosineReadWalker;


/**
 * @author benb
 * 
 * Pair<Frac HCG, Frac GCH>
 * Map< Pair<Frac HCG, Frac GCH>, count>
 * 
 */
public class GnomeSeqFracMethFinderCytosineReadWalker extends
		ReadWalkerToBisulfiteCytosineReadWalker<Integer, Long> {

    @Argument(fullName = "hcgMin", shortName = "hcgmin", doc = "Minimum methylation across HCGs in a read (0-1, default 0)", required = false)
    public double hcgMin = 0.0;
    @Argument(fullName = "hcgMax", shortName = "hcgmax", doc = "Maximum methylation across HCGs in a read (0-1, default 1)", required = false)
    public double hcgMax = 1.0;
    @Argument(fullName = "gchMin", shortName = "gchmin", doc = "Minimum methylation across GCHs in a read (0-1, default 0)", required = false)
    public double gchMin = 0.0;
    @Argument(fullName = "gchMax", shortName = "gchmax", doc = "Maximum methylation across GCHs in a read (0-1, default 1)", required = false)
    public double gchMax = 1.0;

    @Argument(fullName = "minNumHcg", shortName = "minhcg", doc = "Minimum number of HCGs in the read to count (default 3)", required = false)
    public int minNumHcg = 3;
    @Argument(fullName = "minNumGch", shortName = "mingch", doc = "Minimum number of GCHs in the read to count (default 3)", required = false)
    public int minNumGch = 3;

    public GnomeSeqFracMethFinderCytosineReadWalker() {
		super();
	}

	@Override
	protected void alertNewContig(String newContig) {
	}

	/***************************************************
	 * cytosine read walker overrides
	 ***************************************************/
	
	@Override
	public Long treeReduce(Long arg0, Long arg1)
	{
		// **** FIX LATER ***
		System.err.println("GnomeSeqAutocorrByReadWalker does not yet implement multi-threaded mode. Use -nt 1");
		System.exit(1);
		return null;
	}

	@Override
	protected Integer processReadCytosines(List<Cpg> cs) {
		
		FracPair outPair = new FracPair(false);
		for (Cpg c : cs)
		{
			int positionInRead = c.chromPos;
			String cContext = c.context();
			boolean isMethylated = (c.cReads>0);
			
			if (cContext.equalsIgnoreCase("GCH"))
			{
				outPair.incrementGCH(isMethylated);
			}
			else if (cContext.equalsIgnoreCase("HCG"))
			{
				outPair.incrementHCG(isMethylated);
			}
			else
			{
				// Ignore other context
				//System.err.printf("Ignoring cytosine: %d, %s, %s\n", positionInRead, cContext, isMethylated);
			}
		}

		// Check if we're in range and output. HCG is first
		if ( (outPair.hcg().getDenominator() >= this.minNumHcg) && (outPair.gch().getDenominator() >= this.minNumGch) )
		{
			if ( (outPair.hcg().doubleValue() >= this.hcgMin) && (outPair.hcg().doubleValue() <= this.hcgMax ) &&
					(outPair.gch().doubleValue() >= this.gchMin) && (outPair.gch().doubleValue() <= this.gchMax ) )
			{
				// Bonanza.  Just use the middle of the cytosines as the point.  They should be ordered.
				int midIndex = (int)Math.round((double)cs.size()/2.0);
				int point = cs.get(midIndex).chromPos;
				String chr = this.prevContig;
				
				// Use the methyldb schema
				out.printf("%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
						chr, point, "+", 1, 1, 0, 0, 0, 0, 0, 1, 1, 0);
			}
		}
		
		return 1;
	}

	@Override
	public Long reduceInit() {
		return 0L;
	}

	@Override
	public Long reduce(Integer inNew, Long inOld) {
		return inNew + inOld;
	}

	@Override
	public void initialize() {
		super.initialize();
	}

	@Override
	public void onTraversalDone(Long result) {
		super.onTraversalDone(result);
		
		//logger.info(String.format("Found %d total cytosines\n",result));
		logger.info(String.format("onTraversalDone NUmber seen=%d\n",result));
	}

	
	// -- // -- // -- // -- // -- 
	
	
	public class FracPair extends Pair<FractionNonidentical,FractionNonidentical> implements Comparable<FracPair>
	{
		public FracPair(boolean mergeIdentical) {
			super(new FractionNonidentical(), new FractionNonidentical());
			this.first.setUseIdentical(mergeIdentical);
			this.second.setUseIdentical(mergeIdentical);
		}

		public FracPair(FractionNonidentical hcg, FractionNonidentical gch, boolean mergeIdentical) {
			super(hcg, gch);
			this.first.setUseIdentical(mergeIdentical);
			this.second.setUseIdentical(mergeIdentical);
		}
		
		public void incrementHCG(boolean methylated)
		{
			increment(false, methylated);
		}
		
		public void incrementGCH(boolean methylated)
		{
			increment(true, methylated);
		}
		
		public FractionNonidentical hcg()
		{
			return this.first;
		}
		
		public FractionNonidentical gch()
		{
			return this.second;
		}

		private void increment(boolean gch, boolean methylated)
		{
			FractionNonidentical frac = (gch) ? this.second : this.first;
			frac.incDenominator();
			if (methylated) frac.incNumerator();
		}

		@Override
		public String toString() 
		{
			String str = String.format("HCG=%s, GCH=%s", this.first.toString(), this.second.toString());
			return str;
		}

		@Override
		public int compareTo(FracPair arg0) {
			int result = this.first.compareTo(arg0.first);
			
			if (result == 0)
			{
				result = this.second.compareTo(arg0.second);
			}
			//System.err.printf("CompareTo(%s,%s)=%d\n",this,arg0,result);
			return result;
		}
		
		
		
	}
	
}
