package edu.usc.epigenome.uecgatk.benWalkers.cytosineReadWalkers;

import edu.usc.epigenome.uecgatk.benWalkers.CpgBackedByGatk;
import edu.usc.epigenome.uecgatk.benWalkers.ReadWalkerToBisulfiteCytosineReadWalker;

public class CytosineReadTEMPLATEWalker extends
		ReadWalkerToBisulfiteCytosineReadWalker<Integer, Integer> {

	public CytosineReadTEMPLATEWalker() {
		super();
		// TODO Auto-generated constructor stub
	}

	@Override
	protected void alertNewContig(String newContig) {
		// TODO Auto-generated method stub
		
	}

	@Override
	protected Integer processReadCytosines(CpgBackedByGatk thisC) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	protected Integer reduceReadCytosines(Integer value, Integer sum) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Integer reduce(Integer value, Integer sum) {
		// TODO Auto-generated method stub
		return super.reduce(value, sum);
	}

	
	
	
}
