package edu.usc.epigenome.uecgatk.benWalkers.cytosineWalkers;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;

import edu.usc.epigenome.uecgatk.benWalkers.CpgBackedByGatk;
import edu.usc.epigenome.uecgatk.benWalkers.LocusWalkerToBisulfiteCytosineWalker;


/**
 * @author benb
 * 
 * Generates MethylDb table:
 * 
   CREATE TABLE `TEMPLATE_CHR` (
  `chromPos` INT UNSIGNED NOT NULL,
  `strand` enum('+','-') NOT NULL,
  `totalReads` SMALLINT UNSIGNED NOT NULL,
  `cReads` SMALLINT UNSIGNED NOT NULL,
  `cReadsNonconversionFilt` SMALLINT UNSIGNED NOT NULL,
  `tReads` SMALLINT UNSIGNED NOT NULL,
  `agReads` SMALLINT UNSIGNED NOT NULL,

  `totalReadsOpposite` SMALLINT UNSIGNED NOT NULL,
  `aReadsOpposite` SMALLINT UNSIGNED NOT NULL,

  `nextBaseGreads` SMALLINT UNSIGNED NULL default '0',
  `nextBaseTotalReads` SMALLINT UNSIGNED NULL default '0',

  `cpgWeight` SMALLINT UNSIGNED NOT NULL,

  -- Must have a primary key to use updatable rows
  PRIMARY KEY chromPos(chromPos)
 )
	
 * 
 */
public class CytosineToMethylDbWalker extends LocusWalkerToBisulfiteCytosineWalker<Integer, Long> {



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
	protected void alertNewContig(String newContig) 
	{
	}

	
	@Override
	protected Integer processCytosine(CpgBackedByGatk thisC)
	{
		String context = thisC.context();
		AlignmentContext ac = thisC.getAlignmentContext();
		double meth = thisC.fracMeth(false);
		return 1;
	}





	@Override
	protected Long reduceCytosines(Integer value, Long sum) 
	{
		return value.longValue() + sum;
	}
	
}
