package edu.usc.epigenome.uecgatk.qcmetrics.loci;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
//import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import java.io.PrintStream;



/**
 * Zack Ramjan
 * USC Epigenome Center 
 * 03/23/2011
 */

/**
 * downsample dups walker
 */
@By(DataSource.REFERENCE)
public class BinDepthsWalker extends LocusWalker<Boolean,Boolean>  
{
    @Output
    PrintStream out;
    
    @Argument(fullName="winSize", shortName="winsize", doc="window width", required=false)
    protected int WINSIZE = 50000;
    
    protected  long count;
    protected  long current_window;
    protected int current_contig;

    public void initialize() 
    {
    	 count=0;
    	 current_contig = -1;
    	 current_window = 0;
    }
    
    /**
     * The map function runs once per single-base locus, and accepts a 'context', a
     * data structure consisting of the reads which overlap the locus, the sites over
     * which they fall, and the base from the reference that overlaps.
     * 
     * Our map and reduce data types is an array of ints (counts) for the number or trials, and their respective 
     * dup counts. for these two lists, we use a single array and a numtrials+n offest.
     * 
     * @param tracker The accessor for reference metadata.
     * @param ref The reference base that lines up with this locus.
     * @param context Information about reads aligning to this locus.
     * @return return average dup count of all trials.
     */
    @Override
    public Boolean map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) 
    {
    	int contig = ref.getLocus().getContigIndex();
    	long pos = context.getPosition();
    	int depth = context.getBasePileup().getBases().length;
    	
    	if(current_contig != contig || pos % WINSIZE == 0)
    	{
    		out.printf("contig=%d window=%d count=%d%n", current_contig, current_window, count);
    		count=0;
    		current_window = pos;
    		current_contig=contig;
    		
    	}
    	count += depth;
    	return null;
    }

    
    
    /**
     * Provides an initial value for the reduce function. and array of 0's
     * @return 0fill array.
     */
    @Override
    public Boolean reduceInit() 
    { 
    	return true;
    }

    /**
     * Combines the result of the latest map with the accumulator.  In inductive terms,
     * this represents the step loci[x + 1] = loci[x] + 1
     * @param value result of the map.
     * @param sum accumulator for the reduce.
     * @return The total count of loci processed so far.
     */
    @Override
    public Boolean reduce(Boolean value, Boolean val2) {
        
    	return true;
    }

    /**
     * Retrieves the final result of the traversal.
     * @param result The ultimate value of the traversal, produced when map[n] is combined with reduce[n-1]
     *               by the reduce function. 
     */
    @Override
    public void onTraversalDone(Boolean t) 
    {
    	
    }


}