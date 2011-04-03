package edu.usc.epigenome.uecgatk.benWalkers.cytosineReadWalkers;

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
		ReadWalkerToBisulfiteCytosineReadWalker<Pair<FractionNonidentical,FractionNonidentical>, Map<Pair<FractionNonidentical,FractionNonidentical>,Integer>> {

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
	public Map<Pair<FractionNonidentical,FractionNonidentical>,Integer> treeReduce(Map<Pair<FractionNonidentical,FractionNonidentical>,Integer> arg0, 
			Map<Pair<FractionNonidentical,FractionNonidentical>,Integer> arg1) {
		// **** FIX LATER ***
		System.err.println("GnomeSeqAutocorrByReadWalker does not yet implement multi-threaded mode. Use -nt 1");
		System.exit(1);
		return null;
	}

	@Override
	protected Pair<FractionNonidentical,FractionNonidentical> processReadCytosines(List<Cpg> cs) {
		
		int count = 0;
		FractionNonidentical hcg = new FractionNonidentical();
		int denom = 0;
		for (Cpg c : cs)
		{
			int positionInRead = c.chromPos;
			String cContext = c.context();
			boolean isMethylated = (c.cReads>0);
			// out.printf("Processing cytosine: %d, %s, %s\n", positionInRead, cContext, isMethylated);


		}
		return new Integer(count);
	}

	@Override
	public Map<Pair<FractionNonidentical,FractionNonidentical>,Integer> reduceInit() {
		//logger.info(String.format("reduceInit\n"));
		//return null;
	}

	@Override
	public Map<Pair<FractionNonidentical,FractionNonidentical>,Integer> 
	reduce(Pair<FractionNonidentical,FractionNonidentical> arg0, 
			Map<Pair<FractionNonidentical,FractionNonidentical>,Integer> arg1) {
		if (arg0 == null) arg0 = 0;
		//logger.info(String.format("reduce(%s,%s)\n",arg0,arg1));
		return (new Long(arg0))+arg1;
	}

	@Override
	public void initialize() {
		super.initialize();
	}

	@Override
	public void onTraversalDone(Map<Pair<FractionNonidentical,FractionNonidentical>,Integer> result) {
		super.onTraversalDone(result);
		
		logger.info(String.format("Found %d total cytosines\n",result));
	}

	
	
}
