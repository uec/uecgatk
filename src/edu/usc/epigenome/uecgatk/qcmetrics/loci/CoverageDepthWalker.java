package edu.usc.epigenome.uecgatk.qcmetrics.loci;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.commandline.Output;
import java.io.PrintStream;
import org.apache.commons.math.stat.descriptive.*;



/**
 * Zack Ramjan
 * USC Epigenome Center 
 * 06/10/2011
 */

/**
 * Bin Depths walker. calculate coverage in windows across genome.
 * report stats upon these windows
 */
@By(DataSource.REFERENCE)
public class CoverageDepthWalker extends LocusWalker<Boolean,Boolean>  
{
    @Output
    PrintStream out;
    protected  long count;
    protected  long current_window;
    DescriptiveStatistics stats;
   

    public void initialize() 
    {
    	 stats = new DescriptiveStatistics();
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
    	out.println("loci visited=" + stats.getN());
    	out.println("mean=" + stats.getMean());
    	out.println("max=" + stats.getMax());
    	out.println("std dev=" + stats.getStandardDeviation());
    	for(double i=10.0; i<=100.0; i+=10.0)
    		out.println(i + " percentile=" + stats.getPercentile(i));    	
    }


}