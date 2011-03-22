package edu.usc.epigenome.uecgatk.uscec;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;

import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;

public class CpgBackedByGatk extends Cpg {

	private AlignmentContext alignmentContext = null;
	private RefMetaDataTracker metaData = null;
	private ReferenceContext refContext = null;
	
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
	
	
	

}
