package edu.usc.epigenome.uecgatk.benWalkers.cytosineReadWalkers;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math.fraction.Fraction;
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
public class FractionByContextCytosineReadWalker extends
		ReadWalkerToBisulfiteCytosineReadWalker<FractionByContextCytosineReadWalker.FracPair, Map<FractionByContextCytosineReadWalker.FracPair,Integer>> {

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
	protected FracPair processReadCytosines(List<Cpg> cs) {
		
		FracPair outPair = new FracPair();
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
				out.printf("Ignoring cytosine: %d, %s, %s\n", positionInRead, cContext, isMethylated);
			}
		}

		return outPair;
	}

	@Override
	public Map<FracPair,Integer> reduceInit() {
		//logger.info(String.format("reduceInit\n"));
		//return null;
		Map<FracPair,Integer> out = new HashMap<FracPair,Integer>();
		return out;
	}

	@Override
	public Map<FracPair,Integer> 
	reduce(FracPair inNew, 
			Map<FracPair,Integer> inOld) {
		//logger.info(String.format("reduce(%s,%s)\n",arg0,arg1));
		
		Map<FracPair,Integer> out = new HashMap<FracPair,Integer>();
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
	}

	@Override
	public void onTraversalDone(Map<FracPair,Integer> result) {
		super.onTraversalDone(result);
		
		//logger.info(String.format("Found %d total cytosines\n",result));
		for (FracPair pair : result.keySet())
		{
			out.printf("%s: count=%d\n", pair, result.get(pair));
		}
	}

	
	
	
	public class FracPair extends Pair<FractionNonidentical,FractionNonidentical>
	{
		public FracPair() {
			super(new FractionNonidentical(), new FractionNonidentical());
		}

		public FracPair(FractionNonidentical hcg, FractionNonidentical gch) {
			super(hcg, gch);
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
		
		
		
	}
	
}
