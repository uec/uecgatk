package edu.usc.epigenome.uecgatk.qcmetrics.loci;


import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.NotPrimaryAlignmentFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.Random;

import org.apache.commons.math3.stat.descriptive.*;

import java.util.concurrent.atomic.AtomicLongArray;

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
@ReadFilters( {NonUniqueFilter.class,NotPrimaryAlignmentFilter.class} ) // Filter out all reads with zero mapping quality
public class DownsampleCoverageWalker extends LocusWalker<Boolean,Boolean> implements TreeReducible<Boolean>
{
    @Output
    PrintStream out;
    
  
    @Argument(fullName="dumpvals", shortName="dumpv", doc="dumps coverage for each window", required=false)
    protected Boolean dump = false;
    @Argument(fullName="downsample", shortName="p", doc="number of reads to sample", required=false)
    protected long SAMPLESIZE = 0;
    
    protected  long count;
    protected  long current_window;
    protected int current_contig;
    int count_index;
    protected String current_contigName;
    protected int WINSIZE=1;
    
    SynchronizedSummaryStatistics statsNoMem;
    protected double PROBABILITY = 1.0;
    AtomicLongArray coverage = new AtomicLongArray(100000);
    
    Random rand = new Random();

    public void initialize() 
    {
    	 count=0;
    	 current_contig = -1;
    	 current_contigName = "";
    	 current_window = 0;
    	 count_index = 0;
    	 statsNoMem = new SynchronizedSummaryStatistics();
    	 if (dump  && this.getToolkit().getArguments().numberOfDataThreads > 1)
    	 {
    		 System.err.println("Cannot create wigs when running in parallel mode (ex: -nt 4)");
    		 dump = false;
    	 }
    	 if(dump)
    		 out.println("track type=wiggle_0 name=\"coverage\" description=\"winsize=" + WINSIZE + "\"");
    	 
    	if(SAMPLESIZE > 0)
    	{
    		Long epoch = System.currentTimeMillis();
        	String tmpFileName = "tmpLineCounterStats" + epoch.toString() + ".txt";
        	
	    	CommandLineGATK readcount = new CommandLineGATK();
	    	
	    	String region = this.getToolkit().getArguments().intervalArguments.intervals.get(0).toString();
	    	
	    	String[] countargs = {"-T", "DownsampleCoverageWalker", "-R", this.getToolkit().getArguments().referenceFile.getPath(),"-nt",new Integer(this.getToolkit().getArguments().numberOfDataThreads).toString(),"-U","ALLOW_N_CIGAR_READS","-L",region, "-I", this.getToolkit().getArguments().samFiles.get(0), "-o", tmpFileName };
	     	try
	 		{
	 			CommandLineGATK.start(readcount, countargs);
	 	
	 			// Open the file that is the first 
	 			// command line parameter
	 			FileInputStream fstream = new FileInputStream(tmpFileName);
	 			// Get the object of DataInputStream
	 			DataInputStream in = new DataInputStream(fstream);
	 			BufferedReader br = new BufferedReader(new InputStreamReader(in));
	 			String strLine;
	 			//Read File Line By Line
	 			while ((strLine = br.readLine()) != null)   
	 				if(strLine.contains("#window sum="))
	 				{
	 					System.err.println ("total bases found: " + strLine);
	 					PROBABILITY = 1.0 * SAMPLESIZE / Double.parseDouble(strLine.replace("#window sum=", ""));
	 					System.err.println("Downsampling to: " + SAMPLESIZE + " / " + Double.parseDouble(strLine.replace("#window sum=", "")) + " = " + PROBABILITY );
	 				}
	 			
	 			//Close the input stream
	 			in.close();
	 				  
	 		} catch (Exception e)
	 		{
	 			// TODO Auto-generated catch block
	 			e.printStackTrace();
	 		}
    	 
    	}
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
    	
    	
    	int depth = context.getBasePileup().getBases().length;
    	if(SAMPLESIZE > 0)
    	{   		
    		int dsDepth = 0;
    		
    		for(int i = 0; i < depth; i++)
    			if(rand.nextDouble() < PROBABILITY)
    				dsDepth++;
    		depth = dsDepth;
    	}
    	
    	coverage.incrementAndGet(depth);
    	statsNoMem.addValue(1.0 * depth);
    	
    	return null;
    }

    
    //Reduce is not used in this analysis, so we fudge them into being true;
    @Override
    public Boolean reduceInit() 
    { 
    	return true;
    }

    @Override
    public Boolean reduce(Boolean value, Boolean val2) {
        
    	return true;
    }
    
	@Override
	public Boolean treeReduce(Boolean arg0, Boolean arg1)
	{
		return true;
	}

	
	
    /**
     * Retrieves the final result of the traversal.
     * @param just a bogus param for the override. 
     */
    @Override
    public void onTraversalDone(Boolean t) 
    {
    
    	
    	double percentile = .02;
    	long total = 0;
    	for(int i = 0; i <  coverage.length(); i++)
    	{
    		total += coverage.get(i);
    		while(total > (statsNoMem.getN() * percentile))
    		{
    			out.println("#"+ Math.round(percentile * 100) + " percentile=" +  i);
    			percentile += 0.02;
    		}
    	}
    	
    	for(int i = 0; i <  coverage.length(); i++)
    	{
    		total += coverage.get(i);
    		while(total > (statsNoMem.getN() * percentile))
    		{
    			out.println("#"+ Math.round(percentile * 100) + " hist_percentile=" +  coverage.get(i));
    			percentile += 0.02;
    		}
    	}
    	
    	out.println("#window sum=" + statsNoMem.getSum());
    	out.println("#window size=" + WINSIZE);
    	out.println("#window count=" + statsNoMem.getN());
    	out.println("#mean=" + statsNoMem.getMean());
    	out.println("#max=" + statsNoMem.getMax());
    	out.println("#std dev=" + statsNoMem.getStandardDeviation());
    }
    
    

	
 }