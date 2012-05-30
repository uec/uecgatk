package edu.usc.epigenome.uecgatk.bisulfitePileupWalkers;

import edu.usc.epigenome.uecgatk.bisulfitePileupEngines.BisulfitePileupEngine;
import edu.usc.epigenome.uecgatk.pileup.BisulfitePileup;

import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.gatk.filters.MappingQualityFilter;
import org.broadinstitute.sting.gatk.filters.NotPrimaryAlignmentFilter;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;


/**
 * Transforms a normal GATK walker into a BisulfitePileupWalker.  Uses a particular
 * bisulfite analysis engine to transform alignmentContext and refContext into 
 * meaningful bisulfite contexts (read group aware)
 * 
 * Engine must specify whether or not it works with multi-processor (out of order shards)
 * 
 * To do: Must implement tree-reducible to get parallel execution.
 */

//@ReadFilters( {MappingQualityReadFilter.class, BadMateFilter.class, NotPrimaryAlignmentReadFilter.class} ) // Filter out all reads with zero mapping quality
//@Requires( {DataSource.READS, DataSource.REFERENCE, DataSource.REFERENCE_BASES} ) // This walker requires both -I input.bam and -R reference.fasta

public abstract class LocusWalkerToBisulfitePileupWalker<MapType,ReduceType> extends LocusWalker<MapType,ReduceType> implements TreeReducible<ReduceType> {

	
	BisulfitePileupEngine engine = null;
	
	
	/**
	 * 
	 */
	public LocusWalkerToBisulfitePileupWalker() {
		super();
		
		// Check sharding strategy
		//this.getToolkit().
	}
	
	

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.Walker#initialize()
	 */
	@Override
	public void initialize() {
		// TODO Auto-generated method stub
		super.initialize();
		
	}
	
	
	
    /**
     * The map function runs once per single-base locus, and accepts a 'context', a
     * data structure consisting of the reads which overlap the locus, the sites over
     * which they fall, and the base from the reference that overlaps.
     * @param tracker The accessor for reference metadata.
     * @param ref The reference base that lines up with this locus.
     * @param context Information about reads aligning to this locus.
     * @return In this case, returns a count of how many loci were seen at this site (1).
     */
    @Override
    public MapType map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

    	MapType mapout = null;
    	
    	BisulfitePileup bsPileup = engine.gatkLocusToBisulfitePileup(tracker, ref, context, this.getToolkit());
    	

		mapout = processPosition(bsPileup);

		return mapout;
    }

   abstract protected void alertNewContig(String newContig);
   abstract protected MapType processPosition(BisulfitePileup thisC);

   
	/**
	 * Combines the result of the latest map with the accumulator.  In inductive terms,
	 * this represents the step loci[x + 1] = loci[x] + 1
	 * @param value result of the map.
	 * @param sum accumulator for the reduce.
	 * @return The total count of loci processed so far.
	 * 
	 * We implement this to ignore non-cytosines.  Subclasses should implement reduceCytosines
	 */
	@Override
	public ReduceType reduce(MapType value, ReduceType sum) {
		ReduceType out = sum;
		if (value != null)
		{
			out = reducePositions(value, sum);
		}
		return out;
	}

	abstract protected ReduceType reducePositions(MapType value, ReduceType sum);




}