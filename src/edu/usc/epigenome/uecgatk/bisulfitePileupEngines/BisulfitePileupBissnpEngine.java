/**
 * 
 */
package edu.usc.epigenome.uecgatk.bisulfitePileupEngines;

import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;

import edu.usc.epigenome.uecgatk.nmerContexts.NmerCollection;
import edu.usc.epigenome.uecgatk.pileup.BisulfitePileup;
import edu.usc.epigenome.uecgatk.pileup.BisulfitePileupContextPosteriors;

/**
 * @author benb
 *
 * Question for Yaping.  Could this class simply "extend" an existing BisSNP class, and then implment the key
 * method gatkLocusToBisulfitePileup?  I'll look how you do it in your BisSNPUtils, but I suspect this is possible.
 */
public class BisulfitePileupBissnpEngine implements BisulfitePileupEngine {

	
	/***************************
	/*** Constructors 
	/***************************/

	public BisulfitePileupBissnpEngine(NmerCollection contexts) // Pass in priors, BAC, everything else BisSNP needs to run  
	{
	}

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
			AlignmentContext context, GenomeAnalysisEngine toolkit) 
	{
		// Determine posteriors for each read group
		BisulfitePileupContextPosteriors out = new BisulfitePileupContextPosteriors(true);

		
		// Determine C/T pileups for each read group and add to pileup
		//out.add(ReadGroupPileups);
		
		return out;
	
	}



}
