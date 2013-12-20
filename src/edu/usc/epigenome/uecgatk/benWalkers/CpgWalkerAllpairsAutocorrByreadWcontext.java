package edu.usc.epigenome.uecgatk.benWalkers;

import edu.usc.epigenome.genomeLibs.IupacPatterns;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerParams;
import edu.usc.epigenome.uecgatk.benWalkers.cytosineWalkers.Cpg;


public class CpgWalkerAllpairsAutocorrByreadWcontext extends
		CpgWalkerAllpairsAutocorrByread {

	protected String fromContext = null;
	protected String toContext = null;
	private int numPre = 0;
	private int numPost = 0;
	IupacPatterns patternMap = new IupacPatterns();
	
	public CpgWalkerAllpairsAutocorrByreadWcontext(
			CpgWalkerParams inWalkParams, boolean inSamestrandOnly,
			boolean inOppstrandOnly, boolean inSameRead, boolean inDifferentRead,
			String inFromContext, String inToContext) {
		super(inWalkParams, inSamestrandOnly, inOppstrandOnly, inSameRead,
				inDifferentRead);
		
		fromContext = inFromContext;
		toContext = inToContext;
		patternMap.register(fromContext);
		patternMap.register(toContext);
		
		if (fromContext.length() == 2)
		{
			numPre = 0;
			numPost = 1;
		}
		else if (fromContext.length() == 3)
		{
			numPre = 1;
			numPost = 1;
		}
	}

	@Override
	protected void recordPair(Cpg a, Cpg b) {
		
		
		
		String ac = patternMap.firstMatch(a.context(numPre, numPost));
		String bc = patternMap.firstMatch(b.context(numPre, numPost));
		
		if ((ac==null) || (bc==null))
		{
//			System.err.printf("%s,%s context MISSING (%s,%s) should be (%s,%s)\n",
//					a.chromPos, b.chromPos, ac, bc, fromContext, toContext);
		}
		else if (ac.equals(this.fromContext) && bc.equals(this.toContext))
		{
//			System.err.printf("%s,%s context match (%s,%s)\n",a.chromPos, b.chromPos, ac, bc);
			super.recordPair(a, b);
		}
		else if (ac.equals(this.toContext) && bc.equals(this.fromContext))
		{
			super.recordPair(b, a);
		}

		
//				else
//				{
//					System.err.printf("%s,%s context MISMATCH (%s,%s) should be (%s,%s)\n",
//							a.chromPos, b.chromPos, ac, bc, fromContext, toContext);
//				}
	}
	
	
	

}
