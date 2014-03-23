package edu.usc.epigenome.uecgatk.qcmetrics.loci;

import net.sf.samtools.SAMRecord;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
//import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Random;


/**
 * Zack Ramjan
 * USC Epigenome Center 
 * 03/23/2011
 */

/**
 * downsample dups walker
 */
@By(DataSource.READS)
public class DownsampleDupsWalker extends LocusWalker<Integer[],Integer[]>  implements TreeReducible<Integer[]>
{
    @Output
    PrintStream out;
    
    @Argument(fullName="samplesize", shortName="p", doc="number of reads to sample", required=false)
    protected long SAMPLESIZE = 0L;

    @Argument(fullName="numberTrials", shortName="trials", doc="the number of trials", required=true)
    protected int NUMTRIALS = 10;

    protected double PROBABILITY = 1;
    
    public void initialize() 
    {
    	if(SAMPLESIZE > 0)
    	{
	    	Long epoch = System.currentTimeMillis();
	    	String tmpFileName = "tmpLineCounterStats" + epoch.toString() + ".txt";
	    	System.err.println();
	    	CommandLineGATK readcount = new CommandLineGATK();
	    	String[] countargs = {"-T", "ReadCounter", "-R", this.getToolkit().getArguments().referenceFile.getPath(), "-I", this.getToolkit().getArguments().samFiles.get(0), "-o", tmpFileName };
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
    	else 
    		NUMTRIALS = 1;
		
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
    public Integer[] map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) 
    {
    	Random rand = new Random();
    	
    	Integer[] result = new Integer[NUMTRIALS * 2];
    	for(int i = 0; i < result.length; i++)
    		result[i] = 0;
    	List<GATKSAMRecord> coveringReads = context.getReads();
    	ArrayList<SAMRecord> reads = new ArrayList<SAMRecord>();
    	
    	//Get only +,reads that start here
    	for(SAMRecord read : coveringReads)
    	{
    		if(read.getAlignmentStart() == context.getPosition())
    		{
    			//PE
    			if(read.getReadPairedFlag() && !read.getReadNegativeStrandFlag())
    				reads.add(read);
    			//SR
    			else if (!read.getReadPairedFlag())
    				reads.add(read);
    		}
    	}
    	//randomly sample these reads many times and calc dups
    	for(int i = 0; i < NUMTRIALS; i++)
    	{
	    	ArrayList<DupSamRecord> sampledPositions = new ArrayList<DupSamRecord>();
	    	HashMap<DupSamRecord,Integer> dupPositions = new HashMap<DupSamRecord,Integer>();
	    	//build the random subset and add to hashmap for dup checking
	    	for(SAMRecord read : reads)
	    	{
	    		if(PROBABILITY >= 1 || rand.nextDouble() < PROBABILITY)
	    		{
	    			DupSamRecord r = new DupSamRecord(read);
	    			sampledPositions.add(r);
	    			if(dupPositions.containsKey(r))
	    				dupPositions.put(r,dupPositions.get(r) + 1);
	    			else
	    				dupPositions.put(r,1);
	    		}
	    	}
	    	
	    	//count size of subset and how many dups we had
	    	result[i] = sampledPositions.size();
	    	for(int v : dupPositions.values())
	    		result[i+NUMTRIALS] += (v - 1); 	    	
    	}
    	return result;
    }

    
    
    /**
     * Provides an initial value for the reduce function. and array of 0's
     * @return 0fill array.
     */
    @Override
    public Integer[] reduceInit() 
    { 
       	Integer[] result = new Integer[NUMTRIALS * 2];
    	for(int i = 0; i < result.length; i++)
    		result[i] = 0;
    	return result;
    }

    /**
     * Combines the result of the latest map with the accumulator.  In inductive terms,
     * this represents the step loci[x + 1] = loci[x] + 1
     * @param value result of the map.
     * @param sum accumulator for the reduce.
     * @return The total count of loci processed so far.
     */
    @Override
    public Integer[] reduce(Integer[] value, Integer sum[]) {
        
    	Integer[] result = new Integer[NUMTRIALS * 2];
    	for(int i = 0; i < result.length; i++)
    	{
    		result[i] = value[i] + sum[i];
    	}
    	return result;
    }

    /**
     * Retrieves the final result of the traversal.
     * @param result The ultimate value of the traversal, produced when map[n] is combined with reduce[n-1]
     *               by the reduce function. 
     */
    @Override
    public void onTraversalDone(Integer[] result) 
    {
    	DescriptiveStatistics stats = new DescriptiveStatistics();
    	for(int i = 0; i < NUMTRIALS; i++)
    	{
    		out.printf("sampled=%d dups=%d => d/s = %f%n", result[i], result[i+NUMTRIALS], (1.0 * result[i+NUMTRIALS] / result[i]) );
    		stats.addValue((1.0 * result[i+NUMTRIALS] / result[i]));
    	}
    	out.printf("avg d/s = %f%n", stats.getMean());
    }

	@Override
	public Integer[] treeReduce(Integer[] lval, Integer[] rval)
	{
    	return reduce(lval,rval);
	}
}