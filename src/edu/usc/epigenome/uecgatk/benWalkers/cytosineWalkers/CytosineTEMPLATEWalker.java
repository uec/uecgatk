package edu.usc.epigenome.uecgatk.benWalkers.cytosineWalkers;

import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.uecgatk.benWalkers.LocusWalkerToBisulfiteCytosineWalker;


/**
 * @author benb
 * 
 */
public class CytosineTEMPLATEWalker extends LocusWalkerToBisulfiteCytosineWalker<Integer, Long> {



	/**
	 * locus walker overrides
	 */


	/**
	 * Provides an initial value for the reduce function.  Hello walker counts loci,
	 * so the base case for the inductive step is 0, indicating that the walker has seen 0 loci.
	 * @return 0.
	 */
	@Override
	public Long reduceInit()
	{ 
		Long out = 0L;
		return out;
	}



	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.TreeReducible#treeReduce(java.lang.Object, java.lang.Object)
	 * 
	 * Must override this to get multithreaded capability
	 */
	@Override
	public Long treeReduce(Long lhs, Long rhs) 
	{
		return lhs + rhs;
	}

	/**
	 * Retrieves the final result of the traversal.
	 * @param result The ultimate value of the traversal, produced when map[n] is combined with reduce[n-1]
	 *               by the reduce function. 
	 */
	@Override
	public void onTraversalDone(Long result) 
	{
		out.printf("Saw %d cytosines\n", result);
	}

	
	/***************************************************
	 * cytosine walker overrides
	 ***************************************************/
	
	@Override
	protected Integer processCytosine(Cpg thisC)
	{
		String context = thisC.context();
		double meth = thisC.fracMeth(false);
		return 1;
	}


	@Override
	protected Long reduceCytosines(Integer value, Long sum) 
	{
		return value.longValue() + sum;
	}
	
}
