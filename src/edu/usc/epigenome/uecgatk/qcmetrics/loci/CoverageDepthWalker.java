package edu.usc.epigenome.uecgatk.qcmetrics.loci;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.FailsVendorQualityCheckFilter;
import org.broadinstitute.sting.gatk.filters.NotPrimaryAlignmentFilter;
import org.broadinstitute.sting.gatk.filters.UnmappedReadFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.commandline.Output;

import java.io.PrintStream;

import org.apache.commons.math3.stat.descriptive.*;
import edu.usc.epigenome.uecgatk.filters.NonUniqueFilter;



/**
 * Zack Ramjan
 * USC Epigenome Center 
 * 06/10/2011
 */

/**
 * Bin Depths walker. calculate coverage in windows across genome.
 * report stats upon these windows
 * 
 * NOTE BadMateFilter.  Looks like this doesn't use properly paired flag, but just checks
 * that it's on the same chromosome.
 */
@By(DataSource.REFERENCE)
@ReadFilters( {NonUniqueFilter.class,NotPrimaryAlignmentFilter.class,UnmappedReadFilter.class,FailsVendorQualityCheckFilter.class} ) // Filter out all reads with zero mapping quality
public class CoverageDepthWalker extends LocusWalker<Boolean,Boolean> implements TreeReducible<Boolean>  
{
    @Output
    PrintStream out;
    protected  long count;
    protected  long current_window;
    SynchronizedSummaryStatistics stats;
   

    public void initialize() 
    {
    	
    	 stats = new SynchronizedSummaryStatistics();
    	 
    }
    
    /**
     * We pretty much bypass the whole map/reduce stuff since we use a global to keep track of counts
     * we scan the genome linearly, saving counts and resetting when ever we cross a threshold.
     * since we skipped the M/R, this will not be parallel/treereducable.
     * 
     * basically, we have a single for loop
     * 
     * @param tracker The accessor for reference metadata.
     * @param ref The reference base that lines up with this locus.
     * @param context Information about reads aligning to this locus.
     * @return return average dup count of all trials.
     */
    @Override
    public Boolean map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) 
    {
    	stats.addValue(1.0 * context.getBasePileup().getBases().length);
    	return true;
    	
    }

    
    
    /**
     * Provides an initial value for the reduce function. and array of 0's
     * @return always return true.
     */
    @Override
    public Boolean reduceInit() 
    { 
    	return true;
    }

    /**
     * Combines the result of the latest map with the accumulator.  In inductive terms,
     * this represents the step loci[x + 1] = loci[x] + 1
     * @param just a bogus param for the override. 
     * @param just a bogus param for the override. 
     * @return always truen true.
     */
    @Override
    public Boolean reduce(Boolean value, Boolean val2) {
        
    	return true;
    }

    /**
     * Retrieves the final result of the traversal.
     * @param just a bogus param for the override. 
     */
    @Override
    public void onTraversalDone(Boolean t) 
    {
    	out.println("#loci visited=" + stats.getN());
    	out.println("#mean=" + stats.getMean());
    	out.println("#max=" + stats.getMax());
    	out.println("#std dev=" + stats.getStandardDeviation());
    	out.println("#window sum=" + stats.getSum());    	    	
    }

	@Override
	public Boolean treeReduce(Boolean arg0, Boolean arg1)
	{
		return true;
	}


}