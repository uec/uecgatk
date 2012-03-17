package edu.usc.epigenome.uecgatk.bisulfitePileupEngines;

import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;

import edu.usc.epigenome.uecgatk.pileup.BisulfitePileup;

/**
 * @author benb
 * 
 * Core purpose is to transform an alignmenContext and referenceContext into a bisulfitePileup object.
 * In the future, it will transform a Pileup2 object into a bisulfitePileup object, but Pileup2 spec
 * has not been written in GATK yet.
 * 
 * Engine must specify whether or not it works with multi-processor (out of order shards)
 * 
 * Should we handle context filtering here?  Or leave it up to the walker?
 * 
 */
public interface BisulfitePileupEngine {

	public Boolean supportsMultithreadedMode();
	
	
	
	/**
	 * @param tracker
	 * @param ref
	 * @param context
	 * @param toolkit 
	 * @return null if it is not a valid position of the desired context.
	 */
	public BisulfitePileup gatkLocusToBisulfitePileup(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context, GenomeAnalysisEngine toolkit);
	
}
