package edu.usc.epigenome.uecgatk.benWalkers;

import java.util.ArrayList;
import java.util.List;

import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.collections.Pair;

import edu.usc.epigenome.genomeLibs.GatkBaseUtils;

public class CytosineContextCounter 
{
	private final int STARTINGSIZE = 50; 
	
	protected List<CytosineContext> contexts = new ArrayList<CytosineContext>(STARTINGSIZE);
	
	// For backwards compatibility
	protected Pair<Integer,Integer> prevCountsCompatibility = null;
	protected Pair<Integer,Integer> nextCountsCompatibility = null;

	public CytosineContextCounter() {
		super();
		reset();
	}

	public CytosineContextCounter(List<CytosineContext> inContexts) {
		super();
		contexts = inContexts;
	}

	public void reset()
	{
		contexts = new ArrayList<CytosineContext>(STARTINGSIZE);
	}
	
	public void addContext(CytosineContext inContext)
	{
		if ((this.nextCountsCompatibility != null) || (this.prevCountsCompatibility != null))
		{
			System.err.println("Can not call CytosineContextCounter::addContext() if you've already called setPrevCountsCompatibilty()");
			System.exit(1);
		}
		contexts.add(inContext);
	}
	
	public void incrementPrevCountsCompatibilty(int newTotalCounts, int newGCounts)
	{
		if (this.prevCountsCompatibility == null)
		{
			setPrevCountsCompatibilty(0,0);
		}
		this.prevCountsCompatibility.first += newTotalCounts;
		this.prevCountsCompatibility.second += newGCounts;
	}
	
	public void setPrevCountsCompatibilty(int totalCounts, int gCounts)
	{
		if (contexts.size()>0)
		{
			System.err.println("Can not call CytosineContextCounter::setPrevCountsCompatibilty() if you've already called addContext()");
			System.exit(1);
		}
		prevCountsCompatibility = new Pair<Integer,Integer>(totalCounts, gCounts);
	}
	
	public void setNextCountsCompatibilty(int totalCounts, int gCounts)
	{
		if (contexts.size()>0)
		{
			System.err.println("Can not call CytosineContextCounter::setNextCountsCompatibilty() if you've already called addContext()");
			System.exit(1);
		}
		nextCountsCompatibility = new Pair<Integer,Integer>(totalCounts, gCounts);
	}
	
	
	public int numPrevBases()
	{
		if (this.contexts.size()==0) return 0;
		return this.contexts.get(0).numPrevBases();
	}
	
	
	public int numNextBases()
	{
		if (this.contexts.size()==0) return 0;
		return this.contexts.get(0).numNextBases();
	}

	/**
	 * @param relInd: relativePosition
	 * @param matchingBase Should this be IUPAC codes? Not currently.
	 * Update : Now IUPAC handling is done elsewhere, so we don't really need support for it here.
	 */
	public int getGCount(int relInd)
	{
		return getMatchingCount(relInd, BaseUtils.Base.G.base, true);
	}

	public int getGCount(int relInd, boolean bisulfite)
	{
		return getMatchingCount(relInd, BaseUtils.Base.G.base, bisulfite);
	}
	
	public double getMatchingFrac(int relInd, byte matchingBase, boolean bisulfite)
	{
		return (double)getMatchingCount(relInd, matchingBase,bisulfite) / (double)getTotalCount(relInd);
	}
		
	
	public byte[] getRelativeIndexPileup(int relInd)
	{
		byte[] out = new byte[this.contexts.size()];
		for (int i = 0; i < this.contexts.size(); i++)
		{
			CytosineContext context = this.contexts.get(i);
			byte contextBase = context.contextBaseAtRelativeIndex(relInd);
			out[i] = contextBase;
		}
		return out;
	}
	

	public byte majorityBase(int relInd)
	{
		return GatkBaseUtils.mostFrequentBase(getRelativeIndexPileup(relInd));
	}
	
	public int getMatchingCount(int relInd, byte matchingBase, boolean bisulfite)
	{

		// Why was this here??
		//		if (matchingBase == BaseUtils.Base.G.base)
		//		{

		if ((relInd == 1) && (this.nextCountsCompatibility != null))
		{
			return this.nextCountsCompatibility.second;
		}
		else if ((relInd == -1) && (this.prevCountsCompatibility != null))
		{
			return this.prevCountsCompatibility.second;
		}

		int total = 0;
		for (CytosineContext context : this.contexts)
		{
			byte contextBase = context.contextBaseAtRelativeIndex(relInd);
			// ** // ** // ** // ** DOES NOT WORK WITH IUPAC CODES
			if (contextBase == matchingBase)
			{
				total++;
			}
			else if (bisulfite && (matchingBase == BaseUtils.Base.C.base))
			{
				// If it's a cytosine in bisulfite space, we can also match a T
				if (contextBase == BaseUtils.Base.T.base) total++;
			}
//			System.err.printf("\t\tcontextBase=%d\tmatchingBase=%d\ttotalMatch=%d\n", contextBase,matchingBase,total);
		}
		
		return total;
	}

	public int getTotalCount(int relInd)
	{
		if ((relInd == 1) && (this.nextCountsCompatibility != null))
		{
			return this.nextCountsCompatibility.first;
		}
		else if ((relInd == -1) && (this.prevCountsCompatibility != null))
		{
			return this.prevCountsCompatibility.first;
		}
		
		int total = 0;
		for (CytosineContext context : this.contexts)
		{
			byte contextBase = context.contextBaseAtRelativeIndex(relInd);
			// ** // ** // ** // ** DOES NOT WORK WITH IUPAC CODES
			if(contextBase ==  BaseUtils.Base.A.base || contextBase ==  BaseUtils.Base.C.base || contextBase ==  BaseUtils.Base.T.base || contextBase ==  BaseUtils.Base.G.base)
				total++;
//			switch (contextBase)
//			{
//			case BaseUtils.Base.A.base:
//			case BaseUtils.Base.C.base:
//			case BaseUtils.Base.T.base:
//			case BaseUtils.Base.G.base:
//				total++;
//				break;
//			default:
//				break;	
//			}
		}
		
		return total;
	}

	public int getMatchingHCount(int relInd)
	{
		int total = 0;
		for (CytosineContext context : this.contexts)
		{
			byte contextBase = context.contextBaseAtRelativeIndex(relInd);
			if(contextBase ==  BaseUtils.Base.A.base || contextBase ==  BaseUtils.Base.C.base || contextBase ==  BaseUtils.Base.T.base)
				total++;
			
			// ** // ** // ** // ** DOES NOT WORK WITH IUPAC CODES
//			switch (contextBase)
//			{
//			case BaseUtils.Base.A.base:
//			case BaseUtils.Base.C.base:
//			case BaseUtils.Base.T.base:
//				total++;
//				break;
//			default:
//				break;	
//			}
		}
		
		return total;
	}

	public int getMatchingYCount(int relInd)
	{
		int total = 0;
		for (CytosineContext context : this.contexts)
		{
			byte contextBase = context.contextBaseAtRelativeIndex(relInd);
			if(contextBase ==  BaseUtils.Base.C.base || contextBase ==  BaseUtils.Base.T.base)
				total++;
			// ** // ** // ** // ** DOES NOT WORK WITH IUPAC CODES
//			switch (contextBase)
//			{
//			case BaseUtils.Base.C.base:
//			case BaseUtils.Base.T.base:
//				total++;
//				break;
//			default:
//				break;	
//			}
		}
		
		return total;
	}
}
