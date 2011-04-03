package edu.usc.epigenome.uecgatk.benWalkers.cytosineReadWalkers;

import java.util.List;

import org.broadinstitute.sting.jna.lsf.v7_0_6.LibBat.newDebugLog;
import org.broadinstitute.sting.utils.collections.Pair;

import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.uecgatk.benWalkers.CpgBackedByGatk;
import edu.usc.epigenome.uecgatk.benWalkers.ReadWalkerToBisulfiteCytosineReadWalker;


/**
 * @author benb
 * 
 */
public class CytosineReadTEMPLATEWalker extends
		ReadWalkerToBisulfiteCytosineReadWalker<Integer, Long> {

	public CytosineReadTEMPLATEWalker() {
		super();
	}

	@Override
	protected void alertNewContig(String newContig) {
	}

	/***************************************************
	 * cytosine read walker overrides
	 ***************************************************/
	
	@Override
	public Long treeReduce(Long arg0, Long arg1) {
		// **** FIX LATER ***
		System.err.println("GnomeSeqAutocorrByReadWalker does not yet implement multi-threaded mode. Use -nt 1");
		System.exit(1);
		return null;
	}

	@Override
	protected Integer processReadCytosines(List<Cpg> cs) {
		
		int count = 0;
		for (Cpg c : cs)
		{
			int positionInRead = c.chromPos;
			String cContext = c.context();
			boolean isMethylated = (c.cReads>0);
			out.printf("Processing cytosine: %d, %s, %s\n", positionInRead, cContext, isMethylated);
			count++;
		}
		return new Integer(count);
	}

	@Override
	public Long reduceInit() {
		//logger.info(String.format("reduceInit\n"));
		return 0L;
	}

	@Override
	public Long reduce(Integer arg0, Long arg1) {
		if (arg0 == null) arg0 = 0;
		//logger.info(String.format("reduce(%s,%s)\n",arg0,arg1));
		return (new Long(arg0))+arg1;
	}

	@Override
	public void initialize() {
		super.initialize();
	}

	@Override
	public void onTraversalDone(Long result) {
		super.onTraversalDone(result);
		
		logger.info(String.format("Found %d total cytosines\n",result));
	}

	
	
}
