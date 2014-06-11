package edu.usc.epigenome.uecgatk.qcmetrics.loci;


import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.FailsVendorQualityCheckFilter;
import org.broadinstitute.sting.gatk.filters.NotPrimaryAlignmentFilter;
import org.broadinstitute.sting.gatk.filters.UnmappedReadFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Random;

import org.apache.commons.math3.stat.descriptive.*;

import java.util.concurrent.atomic.AtomicLongArray;
import java.util.concurrent.atomic.AtomicReferenceArray;

import edu.usc.epigenome.uecgatk.WalkerTypes.LocusWalkerUnfiltered;
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
@ReadFilters( {NonUniqueFilter.class,NotPrimaryAlignmentFilter.class,UnmappedReadFilter.class,FailsVendorQualityCheckFilter.class} ) // Filter out all reads with zero mapping quality
public class DownsampleCoverageWalker extends LocusWalkerUnfiltered<Boolean,Boolean> implements TreeReducible<Boolean>
{
    @Output
    PrintStream out;
    
    @Argument(fullName="subsample_increment", shortName="subinc", doc="subsample the reads at increments of this length", required=false)
    protected int PROBABILITY_DISTANCE = 100_000_000;
    
    @Argument(fullName="max_subsamples", shortName="maxsubs", doc="maximum number of subsamples", required=false)
    protected int MAX_PROBABILITIES = 20;
   
    ArrayList<AtomicLongArray> coverages = new ArrayList<>();
    ArrayList<SynchronizedSummaryStatistics> statsNoMemList = new ArrayList<>();
    AtomicReferenceArray<Double> PROBABILITIES = new AtomicReferenceArray<>(MAX_PROBABILITIES+1);
    private Random RAND = new Random();
	private Double TOTALBP = 0d;
    

    public void initialize() 
    {
    
 		Long epoch = System.currentTimeMillis();
    	String tmpFileName = "tmpLineCounterStats" + epoch.toString() + ".txt";
    	
    	CommandLineGATK readcount = new CommandLineGATK();
    	String[] countargs;
    	
    	if(this.getToolkit().getArguments().intervalArguments.intervals != null && this.getToolkit().getArguments().intervalArguments.intervals.size() > 0 )
    	{
    		String[] countargsTMP = {"-T", "CoverageDepthWalker", "-R", this.getToolkit().getArguments().referenceFile.getPath(),"-nt",new Integer(this.getToolkit().getArguments().numberOfDataThreads).toString(),"-U","ALLOW_N_CIGAR_READS","-L",this.getToolkit().getArguments().intervalArguments.intervals.get(0).toString(), "-I", this.getToolkit().getArguments().samFiles.get(0), "-o", tmpFileName };
    		countargs = countargsTMP;
    	}
    	
    	else
    	{
    		String[] countargsTMP = {"-T", "CoverageDepthWalker", "-R", this.getToolkit().getArguments().referenceFile.getPath(),"-nt",new Integer(this.getToolkit().getArguments().numberOfDataThreads).toString(),"-U","ALLOW_N_CIGAR_READS", "-I", this.getToolkit().getArguments().samFiles.get(0), "-o", tmpFileName };
    		countargs = countargsTMP;
    	}
    	
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
 			{
 				if(strLine.contains("#window sum="))
 				{
 					TOTALBP = Double.parseDouble(strLine.replace("#window sum=", ""));
 					System.err.println ("total bases found: " + TOTALBP);
 				}
 			}
 			
 			for(int i = 0; i < MAX_PROBABILITIES && (PROBABILITY_DISTANCE * (i+1)) < TOTALBP; i++)
 			{
 				PROBABILITIES.set(i, (PROBABILITY_DISTANCE * (i+1)) / TOTALBP);
 				System.err.println("Will calculate downsample at: " + (PROBABILITY_DISTANCE * (i+1)) + " / " + TOTALBP + " = " + PROBABILITIES.get(i) );
 				coverages.add(new AtomicLongArray(2000));
 				statsNoMemList.add(new SynchronizedSummaryStatistics());
 			}
 			
 			//add a counter for non-subsampled
 			
 			PROBABILITIES.set(coverages.size(), 1.0);
 			coverages.add(new AtomicLongArray(2000));
			statsNoMemList.add(new SynchronizedSummaryStatistics());
			System.err.println("Will calculate downsample at: " + TOTALBP + " / " + TOTALBP + " = " + PROBABILITIES.get(coverages.size()-1) );
 			
 			//Close the input stream
 			in.close();
	 				  
	 		} catch (Exception e)
	 		{
	 			// TODO Auto-generated catch block
	 			e.printStackTrace();
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
		for(int i=0;i<coverages.size();i++)
		{
			int dsDepth=0;
			for(int j = 0; j < depth; j++)
				if(RAND.nextDouble() < PROBABILITIES.get(i))
					dsDepth++;
		
			coverages.get(i).incrementAndGet(dsDepth);
	    	statsNoMemList.get(i).addValue(1.0 * dsDepth);
		}
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
    	long total = 0;
    	for(int i=0;i<coverages.size();i++)
		{
    		long SampleSize = (long) (PROBABILITIES.get(i) *  TOTALBP);
	    	double percentile = .02;
	    	total = 0;
	    	out.print("#Percentiles " + SampleSize + "=") ;
	    	for(int j = 0; j <  coverages.get(i).length(); j++)
	    	{
	    		total += coverages.get(i).get(j);
	    		while(total > (statsNoMemList.get(i).getN() * percentile))
	    		{
	    			out.print(j + " ");
	    			percentile += 0.02;
	    		}
	    	}
	    	out.println();
		}
    }	
 }