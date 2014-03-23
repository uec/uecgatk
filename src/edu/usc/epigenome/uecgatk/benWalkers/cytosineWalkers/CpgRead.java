/**
 * 
 */
package edu.usc.epigenome.uecgatk.benWalkers.cytosineWalkers;

import org.broadinstitute.sting.utils.BaseUtils;

import edu.usc.epigenome.uecgatk.benWalkers.CytosineContext;



/**
 * @author benb
 *
 */
public class CpgRead {


	// Object vars
	public int readId = 0;
	
	public short cRead = 0;
	public short cReadNonconversionFilt = 0;
	public short tRead = 0;
	public short agRead = 0;

	protected CytosineContext context = new CytosineContext();
	
	// Deprecated
//	public short nextBaseGread = 0;
//	protected char nextBaseUpperCase = '0';
//
//	public short prevBaseGread = 0;
//	protected char prevBaseUpperCase = '0';
	
	
	
	/**
	 * @param readId
	 */
	public CpgRead(int readId) {
		super();
		this.readId = readId;
		
	}

	public CpgRead(int readId, short cRead, short cReadNonconversionFilt,
			short tRead, short agRead, CytosineContext inContext)
	{
		super();
		this.readId = readId;
		this.cRead = cRead;
		this.cReadNonconversionFilt = cReadNonconversionFilt;
		this.tRead = tRead;
		this.agRead = agRead;
		
		this.context = inContext;
	}
		
		
		/**
	 * @param readId
	 * @param cRead
	 * @param cReadNonconversionFilt
	 * @param tRead
	 * @param agRead
	 */
	public CpgRead(int readId, short cRead, short cReadNonconversionFilt,
			short tRead, short agRead, char nextBaseUpperCase) {
		super();
		this.readId = readId;
		this.cRead = cRead;
		this.cReadNonconversionFilt = cReadNonconversionFilt;
		this.tRead = tRead;
		this.agRead = agRead;
		
		this.context = new CytosineContext(0,1);
		this.context.setContextBaseAtRelativeIndex(1, nextBaseUpperCase);
	}

	
	public CpgRead(int readId, short cRead, short cReadNonconversionFilt,
			short tRead, short agRead, char nextBaseUpperCase, char prevBaseUpperCase) {
		super();
		this.readId = readId;
		this.cRead = cRead;
		this.cReadNonconversionFilt = cReadNonconversionFilt;
		this.tRead = tRead;
		this.agRead = agRead;

		this.context = new CytosineContext(1,1);
		this.context.setContextBaseAtRelativeIndex(-1, prevBaseUpperCase);
		this.context.setContextBaseAtRelativeIndex(1, nextBaseUpperCase);
	}


	
	/**
	 * @return the context
	 */
	public CytosineContext getContext() {
		return context;
	}

	public boolean nextBaseGread()
	{
		return context.nextBaseGread();
	}
	
	public boolean validCg(boolean useNonconvFilt)
	{
		return 
			(agRead == 0) &&
			(nextBaseGread()) &&
			(!useNonconvFilt || (cReadNonconversionFilt==0));
	}
	
	public boolean meth(boolean useNonconvFilt)
	{
		return
			(cRead > 0) ||
			(!useNonconvFilt && (cReadNonconversionFilt > 0));
	}
	

	public boolean validNextBase()
	{
		boolean out = false;
		if ((this.context != null) && (this.context.numNextBases()>=1))
		{
			byte base = this.context.contextBaseAtRelativeIndex(1);
			out = (base != BaseUtils.Base.N.base);
		}
		return out;
	}
	
	public boolean validPrevBase()
	{
		boolean out = false;
		if ((this.context != null) && (this.context.numPrevBases()>=1))
		{
			byte base = this.context.contextBaseAtRelativeIndex(-1);
			out = (base != BaseUtils.Base.N.base);
		}
		return out;
	}

	@Override
	public String toString() {

		return String.format("%020d\t%d\t%d\t%d\t%d\t%d\t%c\t%d\t%c", 
						readId,
						cRead,
						cReadNonconversionFilt,
						tRead,
						agRead,
						this.context.prevBaseGread(),
						this.context.contextNextBase(),
						this.context.nextBaseGread(),
						this.context.contextNextBase()
		);
	}
}
