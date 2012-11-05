package edu.usc.epigenome.uecgatk.qcmetrics.loci;


import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
//import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.Random;

import org.apache.commons.math.stat.descriptive.*;



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
public class BinDepthsWalker extends LocusWalker<Boolean,Boolean>  
{
    @Output
    PrintStream out;
    
    @Argument(fullName="winSize", shortName="winsize", doc="window width", required=false)
    protected int WINSIZE = 50000;
    @Argument(fullName="dumpvals", shortName="dumpv", doc="dumps coverage for each window", required=false)
    protected Boolean dump = false;
    @Argument(fullName="downsample", shortName="p", doc="number of reads to sample", required=false)
    protected long SAMPLESIZE = 0;
    
    protected  long count;
    protected  long current_window;
    protected int current_contig;
    int count_index;
    protected String current_contigName;
    
    SummaryStatistics statsNoMem;
    protected double PROBABILITY = 1.0;
    long[] coverage = new long[100000000];

    public void initialize() 
    {
    	 count=0;
    	 current_contig = -1;
    	 current_contigName = "";
    	 current_window = 0;
    	 count_index = 0;
    	 statsNoMem = new SummaryStatistics();
    	 
    	 out.println("track type=wiggle_0 name=\"coverage\" description=\"winsize=" + WINSIZE + "\"");
    	 
    	if(SAMPLESIZE > 0)
    	{
	    	CommandLineGATK readcount = new CommandLineGATK();
	     	String[] countargs = {"-T", "ReadCounter", "-R", this.getToolkit().getArguments().referenceFile.getPath(), "-I", this.getToolkit().getArguments().samFiles.get(0), "-o", "tmpLineCounterStats.txt" };
	     	try
	 		{
	 			CommandLineGATK.start(readcount, countargs);
	 	
	 			// Open the file that is the first 
	 			// command line parameter
	 			FileInputStream fstream = new FileInputStream("tmpLineCounterStats.txt");
	 			// Get the object of DataInputStream
	 			DataInputStream in = new DataInputStream(fstream);
	 			BufferedReader br = new BufferedReader(new InputStreamReader(in));
	 			String strLine;
	 			//Read File Line By Line
	 			while ((strLine = br.readLine()) != null)   
	 			{
	 				System.err.println ("total reads found: " + strLine);
	 				PROBABILITY = 1.0 * SAMPLESIZE / Long.parseLong(strLine) ;
	 				System.err.println("Downsampling to: " + SAMPLESIZE + " / " + strLine + " = " + PROBABILITY );
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
    	int contig = ref.getLocus().getContigIndex();
    	String contigName = ref.getLocus().getContig();
    	long pos = context.getPosition();
    	int depth = context.getBasePileup().getBases().length;
    	if(SAMPLESIZE > 0)
    	{   		
    		int dsDepth = 0;
    		Random rand = new Random();
    		for(int i = 0; i <  context.getBasePileup().getBases().length; i++)
    			if(rand.nextDouble() < PROBABILITY)
    				dsDepth++;
    		depth = dsDepth;
    	}
    	
    	
    	
    	if(current_contig != contig)
    	{
    		current_contig=contig;
    		current_contigName = contigName;
    		if(dump)
    			out.printf("fixedStep  chrom=%s start=%d step=%d%n", current_contigName, pos, WINSIZE);
    		
    	}
    	
    	if(pos % WINSIZE == 0)
    	{
    		if(current_contig != -1)
    		{
    			if(dump)
    				out.printf("%f%n", 1.0 * count / count_index);
    			coverage[(int) Math.round((1.0 * count / count_index))]++;
    			statsNoMem.addValue(1.0 * count / count_index);
    		}
    		count=0;
    		count_index=0;
    		current_window = pos;
    		
    		
    	}
    	count += depth;
    	count_index ++;
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
    	
    	double percentile = .1;
    	long total = 0;
    	for(int i = 0; i <  coverage.length; i++)
    	{
    		total += coverage[i];
    		while(total > (statsNoMem.getN() * percentile))
    		{
    			out.println("#"+ Math.round(percentile * 100) + " percentile=" + i);
    			percentile += 0.1;
    		}
    	}
    }
 }