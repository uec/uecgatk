package edu.usc.epigenome.uecgatk.benWalkers;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerParams;
import edu.usc.epigenome.uecgatk.benWalkers.cytosineWalkers.Cpg;

import java.util.*;

public abstract class CpgWalkerAllpairs extends CpgWalker {
	
	protected boolean samestrandOnly = false;
	protected boolean oppstrandOnly = false;

	/**
	 * @param inWalkParams
	 */
	public CpgWalkerAllpairs(CpgWalkerParams inWalkParams, boolean inSamestrandOnly) {
		super(inWalkParams);
		this.samestrandOnly = inSamestrandOnly;
	}

	public CpgWalkerAllpairs(CpgWalkerParams inWalkParams, boolean inSamestrandOnly, boolean inOppstrandOnly) {
		super(inWalkParams);
		this.samestrandOnly = inSamestrandOnly;
		this.oppstrandOnly = inOppstrandOnly;
	}

	@Override
	protected void processWindow(List<Cpg> inWindow) {
		
		// Compare head node to each prior node in the window.
		int n = inWindow.size();
		if (n < 2) return;
		
		//System.err.println(this.windStr());
		
		Cpg head = inWindow.get(n-1);
		for (int i = 0; i < (n-1); i++)
		{
			Cpg prior = inWindow.get(i);
			
			boolean process = true;
			if (this.samestrandOnly)
			{
				//System.err.printf("Using same strand only: strand1=%s,\tstrand2=%s\n",prior.getStrandStr(), head.getStrandStr());
				process = (prior.getStrand() == head.getStrand());
			}
			else if (this.oppstrandOnly)
			{
				//System.err.printf("Using opp strand only: strand1=%s,\tstrand2=%s\n",prior.getStrandStr(), head.getStrandStr());
				process = (prior.getStrand() != head.getStrand());
			}
			
			if (process)
			{
				//System.err.printf("\tProcessing pair: %d, %d\n", prior.chromPos, head.chromPos);
				
				if (prior.chromPos > head.chromPos)
				{
					// Only valid reason is if they're on different chromosomes.
					System.err.printf("Skipping out of order pair (new chrom?): %d, %d\n", prior.chromPos, head.chromPos);
					this.reset();
					return;
				}
				else if (prior.chromPos == head.chromPos)
				{
					System.err.printf("Why did we get the same CpG twice (%d)??\n",head.chromPos);  
					(new Exception()).printStackTrace();
					System.exit(1);
				}
				else
				{
					recordPair(prior, head);
				}
			}
		}
		
		// And the superclass
		super.processWindow(inWindow);
	}
	
	protected abstract void recordPair(Cpg first, Cpg second);
}
