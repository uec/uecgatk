package edu.usc.epigenome.uecgatk.benWalkers.cytosineWalkers;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import org.broadinstitute.sting.utils.collections.Pair;

import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers.CpgMethLevelSummarizer;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers.CpgMethLevelSummarizer;
import edu.usc.epigenome.uecgatk.benWalkers.LocusWalkerToBisulfiteCytosineWalker;

/**
 * @author benb
 * 
 * Maptype = Pair<String,Double>. String is the cytosine context, Double is the percent meth.
 * Reducetype = Map<String,CpgMethLevelSummarizer>. The key is the cytosine context, the value is a meth summarizer
 *
 */
public class MethLevelAveragesWalker 
	extends LocusWalkerToBisulfiteCytosineWalker<Pair<String,Double>, Map<String,CpgMethLevelSummarizer>> {


//	public MethLevelAveragesWalker() {
//		super();
//	}


	/**
	 * locus walker overrides
	 */


	/**
	 * Provides an initial value for the reduce function.  Hello walker counts loci,
	 * so the base case for the inductive step is 0, indicating that the walker has seen 0 loci.
	 * @return 0.
	 */
	@Override
	public Map<String,CpgMethLevelSummarizer> reduceInit()
	{ 
		Map<String,CpgMethLevelSummarizer> out = 
			new HashMap<String,CpgMethLevelSummarizer>();
		return out;
	}



	@Override
	public Map<String, CpgMethLevelSummarizer> treeReduce(
			Map<String, CpgMethLevelSummarizer> lhs,
			Map<String, CpgMethLevelSummarizer> rhs) 
	{
		Map<String,CpgMethLevelSummarizer> out = 
			new HashMap<String,CpgMethLevelSummarizer>();
		
		Set<String> keys = lhs.keySet();
		keys.addAll(rhs.keySet());
		for (String key : keys)
		{
			if (lhs.containsKey(key) && !rhs.containsKey(key))
			{
				out.put(key, lhs.get(key));
			}
			else if (!lhs.containsKey(key) && rhs.containsKey(key))
			{
				out.put(key, lhs.get(key));
			}
			else
			{
				// They must both contain the key.  Merge
				out.put(key, (CpgMethLevelSummarizer)CpgMethLevelSummarizer.sumSummarizers(lhs.get(key), rhs.get(key)));
			}
		}
		
		return out;
	}

	/**
	 * Retrieves the final result of the traversal.
	 * @param result The ultimate value of the traversal, produced when map[n] is combined with reduce[n-1]
	 *               by the reduce function. 
	 */
	@Override
	public void onTraversalDone(Map<String, CpgMethLevelSummarizer> result) 
	{
		out.println("Number of cytosines viewed is: " + result);
		
		for (String key : result.keySet())
		{
			CpgMethLevelSummarizer summarizer = result.get(key);
			out.printf("%s:\t%d\t%.2f%%\n",key,(int)summarizer.getNumVals(), summarizer.getValMean()*100);
		}
	}

	
	/***************************************************
	 * cytosine walker overrides
	 ***************************************************/
	
	@Override
	protected Pair<String,Double> processCytosine(Cpg thisC)
	{
		String context = thisC.context();
		double meth = thisC.fracMeth(false);
		return new Pair<String,Double>(context, meth);
	}


	@Override
	protected Map<String, CpgMethLevelSummarizer> reduceCytosines(
			Pair<String, Double> value, Map<String, CpgMethLevelSummarizer> sum) 
	{
		String context = value.first;
		CpgMethLevelSummarizer summarizer = null;
		if (sum.containsKey(context))
		{
			summarizer = sum.get(context);
		}
		else
		{
			summarizer = new CpgMethLevelSummarizer();
			sum.put(context, summarizer);
		}
		// Stream cytosine
		summarizer.streamValue(value.second,1);
		
		return sum;
	}
	
}
