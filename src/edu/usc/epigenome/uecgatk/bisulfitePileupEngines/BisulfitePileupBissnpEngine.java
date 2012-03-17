/**
 * 
 */
package edu.usc.epigenome.uecgatk.bisulfitePileupEngines;

import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;

import edu.usc.epigenome.uecgatk.pileup.BisulfitePileup;

/**
 * @author benb
 *
 */
public class BisulfitePileupBissnpEngine implements BisulfitePileupEngine {

	
	// BAC
	
	@Override
	public Boolean supportsMultithreadedMode() {
		return true;
	}

	
	/* (non-Javadoc)
	 * @see edu.usc.epigenome.uecgatk.bisulfitePileupEngines.BisulfitePileupEngine#gatkLocusToBisulfitePileup(org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker, org.broadinstitute.sting.gatk.contexts.ReferenceContext, org.broadinstitute.sting.gatk.contexts.AlignmentContext)
	 */
	@Override
	public BisulfitePileup gatkLocusToBisulfitePileup(
			RefMetaDataTracker tracker, ReferenceContext ref,
			AlignmentContext context, GenomeAnalysisEngine toolkit) {
		// TODO Auto-generated method stub
		return null;
	}



}
