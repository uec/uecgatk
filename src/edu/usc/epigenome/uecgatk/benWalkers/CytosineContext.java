package edu.usc.epigenome.uecgatk.benWalkers;

import org.broadinstitute.sting.utils.BaseUtils;

import edu.usc.epigenome.genomeLibs.GatkBaseUtils;

public class CytosineContext {

	// Uses GATK encoding of bases
	protected byte[] context = null;
	protected int cytosineIndexWithinContext = -1; // 0-offset


	/**
	 * 
	 */
	public CytosineContext() {
		super();
	}


	public CytosineContext(int inNumPre, int inNumPost) {
		super();
		this.resetContext(inNumPre, inNumPost);
	}

	public static CytosineContext fakeContext()
	{
		CytosineContext c = new CytosineContext(1,1);
		
		for (int i = -1; i<= 1; i++)
		{
			c.setContextBaseAtRelativeIndex(i, 'G');
		}
		
//		System.err.printf("fakeContext(%d,%d) output: %s\n", c.numPrevBases(), c.numNextBases(), c.toString());
		
		return c;
	}
	
	/**
	 * @param context
	 * @param cytosineIndexWithinContext
	 */
	public CytosineContext(byte[] context, int cytosineIndexWithinContext) {
		super();
		this.context = context;
		this.cytosineIndexWithinContext = cytosineIndexWithinContext;
	}


	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		String out = null;
		StringBuffer buf = new StringBuffer();
		try
		{
			int curPos = 0;
			for (int i = 0; i < this.numPrevBases(); i++, curPos++)
			{
				buf.append( GatkBaseUtils.CharFromBaseByte(this.contextBaseAtGlobalIndex(curPos) )); 
			}
			buf.append('.');
			// Central base
			buf.append( GatkBaseUtils.CharFromBaseByte(this.contextBaseAtGlobalIndex(curPos++) )); 
			buf.append('.');
			for (int i = 0; i < this.numNextBases(); i++, curPos++)
			{
				buf.append( GatkBaseUtils.CharFromBaseByte(this.contextBaseAtGlobalIndex(curPos) )); 
			}
			out = buf.toString();
		}
		catch (Exception e)
		{
			out = super.toString();
		}
		return out;
	}


	public int numPrevBases()
	{
		return this.cytosineIndexWithinContext;
	}
	
	public int numNextBases()
	{
		return this.cytosineIndexWithinContext;
	}

	public void resetContext(int inNumPre, int inNumPost)
	{
		this.context = CytosineContext.blankContext(inNumPre+inNumPost+1);
		this.cytosineIndexWithinContext = inNumPre;
	}
	
	static public byte[] blankContext(int inLength)
	{
		byte[] out = new byte[inLength];
		for (int i = 0; i < inLength; i++)
		{
			out[i] = BaseUtils.Base.N.base;
		}
		return out;
	}
	
	public byte contextNextBase()
	{
		return this.contextBaseAtRelativeIndex(1);
	}
	
	public byte contextPrevBase()
	{
		return this.contextBaseAtRelativeIndex(-1);
	}

	/**
	 * @param relInd 0 is the cytosine, -1 is previous, +1 is next base, and so on.
	 * @return
	 */
	public byte contextBaseAtRelativeIndex(int relInd)
	{
		int globalInd = this.cytosineIndexWithinContext+relInd;
		return contextBaseAtGlobalIndex(globalInd);
	}
	
	protected byte contextBaseAtGlobalIndex(int inIndex)
	{
		byte base = BaseUtils.Base.N.base;
		try
		{
			if (this.context == null)
			{
				throw new Exception ("Trying to access empty this.context");
			}
			if (this.context.length <= (inIndex))
			{
				throw new Exception (String.format("Trying to access index %d from a %d length this.context array",inIndex, this.context.length));
			}
		}
		catch (Exception e)
		{
			System.err.println(e.toString());
			e.printStackTrace();
			System.exit(1);
		}
		base = this.context[inIndex];
		return base;
	}

	public void setContextBaseAtRelativeIndex(int inIndex, char inBaseChar)
	{
		this.setContextBaseAtRelativeIndex(inIndex, GatkBaseUtils.BaseByteFromChar(inBaseChar));
	}
	
	public void setContextBaseAtRelativeIndex(int inIndex, byte inBase)
	{
		// If we're off the left end, we have to recreate our context array
		if (inIndex < (-this.numPrevBases()))
		{
			System.err.printf("Growing number of prevBases in context from %d to %d", this.numPrevBases(),-inIndex);
			int newLength = (-inIndex) + this.numNextBases() + 1;
			byte[] newContext = CytosineContext.blankContext(newLength);
			System.arraycopy(this.context, 0, newContext, -inIndex-this.numPrevBases(), this.context.length);
			this.context = newContext;
		}

		// And set
		this.setContextBaseAtGlobalIndex(inIndex + this.cytosineIndexWithinContext, inBase);
	}
	
	protected void setContextBaseAtGlobalIndex(int inIndex, byte inBase)
	{
		if (this.context == null)
		{
			context = CytosineContext.blankContext(inIndex+1);
		}
		
		if (this.context.length <= inIndex)
		{
			// Stretch
			System.err.printf("Growing number of nextBases in context from %d to %d", this.numNextBases(),inIndex+1-this.cytosineIndexWithinContext);
			byte[] newContext = CytosineContext.blankContext(inIndex+1);
			System.arraycopy(context, 0, newContext, 0, this.context.length);
			this.context = newContext;
		}
	
		this.context[inIndex] = inBase;
	}
	

	

	public boolean contextBaseIsG(int relInd)
	{
		byte base = contextBaseAtRelativeIndex(relInd);
		return (base == BaseUtils.Base.G.base);
	}
	
	public boolean nextBaseGread()
	{
		return contextBaseIsG(1);
	}
	
	public boolean prevBaseGread()
	{
		return contextBaseIsG(-1);
	}

	
}
