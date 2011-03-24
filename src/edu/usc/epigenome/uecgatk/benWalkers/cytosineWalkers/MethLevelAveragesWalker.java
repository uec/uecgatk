package edu.usc.epigenome.uecgatk.benWalkers.cytosineWalkers;

import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.uecgatk.benWalkers.LocusWalkerToBisulfiteCytosineWalker;

public class MethLevelAveragesWalker extends LocusWalkerToBisulfiteCytosineWalker<Integer, Long> {



	/**
	 * Provides an initial value for the reduce function.  Hello walker counts loci,
	 * so the base case for the inductive step is 0, indicating that the walker has seen 0 loci.
	 * @return 0.
	 */
	@Override
	public Long reduceInit() { return 0L; }





	@Override
	protected Long reduceCytosines(Integer value, Long sum) {
		return value + sum;
	}

	@Override
	public Long treeReduce(Long lhs, Long rhs) {
		// TODO Auto-generated method stub
		return lhs + rhs;
	}

	/**
	 * Retrieves the final result of the traversal.
	 * @param result The ultimate value of the traversal, produced when map[n] is combined with reduce[n-1]
	 *               by the reduce function. 
	 */
	@Override
	public void onTraversalDone(Long result) {
		out.println("Number of cytosines viewed is: " + result);
	}

	@Override
	protected Integer processCytosine(Cpg thisC)
	{
		return 1;
	}
}
