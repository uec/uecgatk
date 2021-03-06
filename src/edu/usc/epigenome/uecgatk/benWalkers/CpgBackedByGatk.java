package edu.usc.epigenome.uecgatk.benWalkers;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.GenomeLoc;
//ZR import org.jfree.util.Log;



import edu.usc.epigenome.uecgatk.benWalkers.cytosineWalkers.Cpg;
import edu.usc.epigenome.uecgatk.benWalkers.cytosineWalkers.CpgRead;

public class CpgBackedByGatk extends Cpg {

	protected AlignmentContext alignmentContext = null;
	protected RefMetaDataTracker metaData = null;
	protected ReferenceContext refContext = null;
	
	public CpgBackedByGatk() {
	}

	public CpgBackedByGatk(int chromPos, boolean negStrand) {
		super(chromPos, negStrand);
	}

	public CpgBackedByGatk(int chromPos, boolean negStrand, AlignmentContext ac, RefMetaDataTracker meta, ReferenceContext rc) {
		super(chromPos, negStrand);
		this.setAlignmentContext(ac);
		this.setMetaData(meta);
		this.setRefContext(rc);
	}

	public CpgBackedByGatk(int chromPos, boolean negStrand, short totalReads,
			short cReads, short cReadsNonconversionFilt, short tReads,
			short agReads, short totalReadsOpposite, short aReadsOpposite,
			int cpgWeight, short nextBaseGreads, short nextBaseTotalReads,
			char nextBaseRefUpperCase) {
		super(chromPos, negStrand, totalReads, cReads, cReadsNonconversionFilt,
				tReads, agReads, totalReadsOpposite, aReadsOpposite, cpgWeight,
				nextBaseGreads, nextBaseTotalReads, nextBaseRefUpperCase);
	}

	public CpgBackedByGatk(int chromPos, boolean negStrand, short totalReads,
			short cReads, short cReadsNonconversionFilt, short tReads,
			short agReads, short totalReadsOpposite, short aReadsOpposite,
			int cpgWeight, CytosineContextCounter inCounter, char inPrevBaseRef, 
			char inNextBaseRef) {
		super(chromPos, negStrand, totalReads, cReads, cReadsNonconversionFilt,
				tReads, agReads, totalReadsOpposite, aReadsOpposite, cpgWeight,
				inCounter, inPrevBaseRef, inNextBaseRef);
	}

	/**
	 * @return the alignmentContext
	 */
	public AlignmentContext getAlignmentContext() {
		return alignmentContext;
	}

	/**
	 * @param alignmentContext the alignmentContext to set
	 */
	public void setAlignmentContext(AlignmentContext alignmentContext) {
		this.alignmentContext = alignmentContext;
	}

	/**
	 * @return the metaData
	 */
	public RefMetaDataTracker getMetaData() {
		return metaData;
	}

	/**
	 * @param metaData the metaData to set
	 */
	public void setMetaData(RefMetaDataTracker metaData) {
		this.metaData = metaData;
	}

	/**
	 * @return the refContext
	 */
	public ReferenceContext getRefContext() {
		return refContext;
	}

	/**
	 * @param refContext the refContext to set
	 */
	public void setRefContext(ReferenceContext refContext) {
		this.refContext = refContext;
	}
	
	public String getChrom()
	{
		return this.refContext.getLocus().getContig();
	}
	
	public GenomeLoc getGenomeLoc()
	{
		GenomeLoc out = this.refContext.getLocus();

		// ***** REMOVE ***** DEBUGGING ONLY **** 
		if (out.getStart() != this.chromPos)
		{
			System.err.printf("CpgBackedByGatk::getGenomeLoc() refContextPos(%d) != chromPos(%d)\n", out.getStart(), this.chromPos);
			System.exit(1);
		}
		
		return out;
	}

	public void addRead(CpgRead cRead)
	{
		// TODO Auto-generated method stub
		
	}
	
	

}
