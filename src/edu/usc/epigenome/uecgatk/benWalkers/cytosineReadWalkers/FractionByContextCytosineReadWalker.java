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

import edu.usc.epigenome.genomeLibs.IupacPatterns;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.uecgatk.FractionNonidentical;
import edu.usc.epigenome.uecgatk.benWalkers.ReadWalkerToBisulfiteCytosineReadWalker;
import edu.usc.epigenome.uecgatk.benWalkers.ReadWithCpgMeths;


/**
 * @author benb
 * 
 * Pair<Frac HCG, Frac GCH>
 * Map< Pair<Frac HCG, Frac GCH>, count>
 * 
 */
public class FractionByContextCytosineReadWalker extends
		ReadWalkerToBisulfiteCytosineReadWalker<FractionByContextCytosineReadWalker.FracPair, Map<FractionByContextCytosineReadWalker.FracPair,Integer>> {

    @Argument(fullName = "mergeEqualVals", shortName = "m", doc = "Merges methylation states with the same vals, i.e. 1/2 and 2/4", required = false)
    public boolean mergeEqualVals = false;

    protected static IupacPatterns patternMap = null;

    public FractionByContextCytosineReadWalker() {
		super();
	}

	@Override
	protected void alertNewContig(String newContig) {
	}

	/***************************************************
	 * cytosine read walker overrides
	 ***************************************************/
	
	@Override
	public Map<FracPair,Integer> treeReduce(Map<FracPair,Integer> arg0, Map<FracPair, Integer> arg1)
	{
		// **** FIX LATER ***
		System.err.println("GnomeSeqAutocorrByReadWalker does not yet implement multi-threaded mode. Use -nt 1");
		System.exit(1);
		return null;
	}

	@Override
	protected FracPair processReadCytosines(ReadWithCpgMeths read) {
		
		List<Cpg> cs = read;
		
		FracPair outPair = new FracPair(this.mergeEqualVals);
		for (Cpg c : cs)
		{
			String contextOrig = c.getPrevBasesRef() + "C" + c.getNextBasesRef();
			String cContext = this.patternMap.firstMatch(contextOrig);
			
			boolean isMethylated = (c.cReads>0);
			
			if (cContext == null)
			{
			}
			else if (cContext.equalsIgnoreCase("GCH"))
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
				System.err.printf("Ignoring cytosine: context=%s, %s\n", cContext, isMethylated);
			}
		}

		return outPair;
	}

	@Override
	public Map<FracPair,Integer> reduceInit() {
		//logger.info(String.format("reduceInit\n"));
		//return null;
		Map<FracPair,Integer> out = new TreeMap<FracPair,Integer>();
		return out;
	}

	@Override
	public Map<FracPair,Integer> 
	reduce(FracPair inNew, 
			Map<FracPair,Integer> inOld) {
		//logger.info(String.format("reduce(%s,%s)\n",arg0,arg1));
		
		Map<FracPair,Integer> out = new TreeMap<FracPair,Integer>();
		out.putAll(inOld);
		
		Integer newVal = new Integer(1);
		if (out.containsKey(inNew))
		{
			newVal = out.get(inNew);
			newVal++;
		}
		out.put(inNew, newVal);
		
		return out;
	}

	@Override
	public void initialize() {
		super.initialize();

		patternMap = new IupacPatterns();
		for (String pat : new String[] {"GCH", "HCG"})
		{
			patternMap.register(pat);
		}

	
	}

	@Override
	public void onTraversalDone(Map<FracPair,Integer> result) {
		super.onTraversalDone(result);
		
		//logger.info(String.format("Found %d total cytosines\n",result));
		logger.info(String.format("onTraversalDone mergeEqualVals=%s\n",this.mergeEqualVals));

		// First get all the fractions. 
		Set<FractionNonidentical> allFracs = new TreeSet<FractionNonidentical>();
		for (FracPair pair : result.keySet())
		{
			FractionNonidentical frac = pair.first;
			allFracs.add(frac);
			frac = pair.second;
			allFracs.add(frac);
//			int thisReads = result.get(pair);
//			totalReads += thisReads;
//			out.printf("%s: count=%d\n", pair, thisReads);
		}
		
		int totalReads = 0;
		for (FractionNonidentical fracY : allFracs)
		{
			out.printf("%.3f,%d,%d,", fracY.doubleValue(), fracY.getNumerator(), fracY.getDenominator());
			int onX = 0;
			for (FractionNonidentical fracX : allFracs)
			{
				FracPair pair = new FracPair(fracX, fracY, this.mergeEqualVals);
				int thisReads = 0;
				if (result.containsKey(pair))
				{
					thisReads = result.get(pair);
				}
				
				if (onX > 0)
				{
					out.print(",");
				}
				out.printf("%d", thisReads);
				onX++;
				totalReads += thisReads;
			}
			out.println();
		}
		logger.info(String.format("Total reads=%d\n", totalReads));
	}

	
	
	
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
