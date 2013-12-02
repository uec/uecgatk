package edu.usc.epigenome.uecgatk.qcmetrics.loci;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import java.io.PrintStream;
import java.util.Random;

import org.apache.commons.math.stat.descriptive.*;

import edu.usc.epigenome.uecgatk.filters.NonUniqueFilter;



/**
 * Zack Ramjan
 * USC Epigenome Center 
 * 06/10/2011
 * 08/08/2012
 */

/**
 * Bin Depths walker. calculate coverage in windows across genome.
 * report stats upon these windows
 */
@By(DataSource.REFERENCE)
@ReadFilters( {NonUniqueFilter.class} ) // Filter out all reads with zero mapping quality
public class ZPeakzWalker extends LocusWalker<Boolean,Boolean>  
{
    @Output
    PrintStream out;
    
    @Argument(fullName="winSize", shortName="winsize", doc="window width", required=false)
    protected int WINSIZE = 10;
    @Argument(fullName="dumpvals", shortName="dumpv", doc="dumps coverage for each window", required=false)
    protected Boolean dump = false;
    
    protected  long count;
    protected int current_contig;
    int current_index;
    protected String current_contigName;
    
    SummaryStatistics statsNoMem;
    protected double PROBABILITY = 1.0;
    long[] coverage = new long[300_000_000];
    Random rand = new Random();

    public void initialize() 
    {
    	 count=0;
    	 current_contig = -1;
    	 current_contigName = "";
    	 current_index = 0;
    	 statsNoMem = new SummaryStatistics();
    	 
    	 out.println("track type=wiggle_0 name=\"coverage\" description=\"winsize=" + WINSIZE + "\"");

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
    	int contig = ref.getLocus().getContigIndex();
    	String contigName = ref.getLocus().getContig();
    	long pos = context.getPosition();
    	int depth = context.getBasePileup().getBases().length;
    	
    	if(current_contig != contig)
    	{
    		if(current_contig >= 0)
    			processContig(pos);
    		
    		//reset counts for new contig
    		current_contig=contig;
    		current_contigName = contigName;
    		current_index = -1;
    		
    		    		
    	}
    	
    	if(pos % WINSIZE == 0)
    		current_index++;
    		
    	else
    		if(depth > coverage[current_index]) 
    			coverage[current_index] = depth;
    	
    	
    	if(pos % WINSIZE == 0)
    	{
    		if(current_contig != -1)
    		{
    			if(dump)
    	
    			coverage[(int) Math.round((10.0 * count / current_index))]++;
    			statsNoMem.addValue(1.0 * count / current_index);
    		}
    		count=0;
    		current_index=0;	
    	}
    	count += depth;
    	current_index ++;
    	return null;
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
    	out.println("#window size=" + WINSIZE);
    	out.println("#window count=" + statsNoMem.getN());
    	out.println("#mean=" + statsNoMem.getMean());
    	out.println("#max=" + statsNoMem.getMax());
    	out.println("#std dev=" + statsNoMem.getStandardDeviation());
    	
    	double percentile = .02;
    	long total = 0;
    	for(int i = 0; i <  coverage.length; i++)
    	{
    		total += coverage[i];
    		while(total > (statsNoMem.getN() * percentile))
    		{
    			out.println("#"+ Math.round(percentile * 100) + " percentile=" + (1.0 * i / 10.0));
    			percentile += 0.02;
    		}
    	}
    }
    
    protected void processContig(long pos)
    {
    	if(dump)
			out.printf("fixedStep  chrom=%s start=%d step=%d span=%d%n", current_contigName, pos, WINSIZE,WINSIZE);
		//out.printf("%f%n", 1.0 * count / current_index);
    }
 }