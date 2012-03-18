package edu.usc.epigenome.uecgatk.nmerContexts;

import org.broadinstitute.sting.utils.BaseUtils;

import edu.usc.epigenome.genomeLibs.GatkBaseUtils;

/**
 * @author benb
 *
 * Used to represent an nmer with a particular position being the "reference".  Allows
 * degeneracy via IUPAC codes
 * Our primary purpose is cytosine contexts (upper case is "reference"):
 * Cg
 * Ch
 * Chh
 * Chg
 * gCh
 * gCg
 * 
 * For making multiple contexts conform to the same number of prefix and suffix positions,
 * notice that you can always make a context longer by using n symbols:
 * Cg.lengthen(Chg) = Cgn 
 */
public class NmerWithRef {
//	BaseUtils.Base
	
	// some presets
	public static NmerWithRef CG = fromRefUpper("Cg");
	
	// Internal representation
	// Uses GATK encoding of bases
	protected byte[] context = null;
	protected int cytosineIndexWithinContext = -1; // 0-offset

	
	
	/***************************
	/*** Constructors 
	/***************************/
	
	public static NmerWithRef fromRefUpper(String in)
	{
		return null; // TO DO
	}
	
	/**
	 * @param in Yaping's notation, "GCG-2"
	 * @return
	 */
	public static NmerWithRef fromRefInt(String in)
	{
		return null; // TO DO
	}

	public static NmerWithRef fakeContext(int numPre, int numPost)
	{
		NmerWithRef c = new NmerWithRef(numPre,numPost);
		
		for (int i = -1; i<= 1; i++)
		{
			c.setContextBaseAtRelativeIndex(i, 'n');
		}
		
//		System.err.printf("fakeContext(%d,%d) output: %s\n", c.numPrevBases(), c.numNextBases(), c.toString());
		
		return c;
	}

	
	/**
	 * @param context
	 * @param cytosineIndexWithinContext
	 */
	public NmerWithRef(byte[] context, int cytosineIndexWithinContext) {
		super();
		this.context = context;
		this.cytosineIndexWithinContext = cytosineIndexWithinContext;
	}

	protected NmerWithRef() {
		super();
	}


	protected NmerWithRef(int inNumPre, int inNumPost) {
		super();
		this.resetContext(inNumPre, inNumPost);
	}

	



	/***************************
	/*** Output 
	/***************************/
	
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

	
	
	/***************************
	/*** Information 
	/***************************/
	

	public int numPrevBases()
	{
		return this.cytosineIndexWithinContext;
	}
	
	public int numNextBases()
	{
		return this.cytosineIndexWithinContext;
	}


	public boolean contextBaseIsG(int relInd)
	{
		byte base = contextBaseAtRelativeIndex(relInd);
		return (base == BaseUtils.G);
	}
	
	public boolean nextBaseG()
	{
		return contextBaseIsG(1);
	}
	
	public boolean prevBaseG()
	{
		return contextBaseIsG(-1);
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
		byte base = BaseUtils.N;
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


	/***************************
	/*** Modify 
	/***************************/
	
	
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
			byte[] newContext = NmerWithRef.blankContext(newLength);
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
			context = NmerWithRef.blankContext(inIndex+1);
		}
		
		if (this.context.length <= inIndex)
		{
			// Stretch
			System.err.printf("Growing number of nextBases in context from %d to %d", this.numNextBases(),inIndex+1-this.cytosineIndexWithinContext);
			byte[] newContext = NmerWithRef.blankContext(inIndex+1);
			System.arraycopy(context, 0, newContext, 0, this.context.length);
			this.context = newContext;
		}
	
		this.context[inIndex] = inBase;
	}
	
	
	public void resetContext(int inNumPre, int inNumPost)
	{
		this.context = NmerWithRef.blankContext(inNumPre+inNumPost+1);
		this.cytosineIndexWithinContext = inNumPre;
	}
	
	
	/***************************
	/*** Utility 
	/***************************/
	
	static protected byte[] blankContext(int inLength)
	{
		byte[] out = new byte[inLength];
		for (int i = 0; i < inLength; i++)
		{
			out[i] = BaseUtils.N;
		}
		return out;
	}
	
	


	
	
}
