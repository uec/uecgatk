package edu.usc.epigenome.uecgatk.qcmetrics.loci;
import net.sf.samtools.SAMReadGroupRecord;

import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.NotPrimaryAlignmentFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.pileup.PileupElement;

import org.broadinstitute.sting.commandline.Output;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;






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
@By(DataSource.READS)
@ReadFilters( {NotPrimaryAlignmentFilter.class} ) // Filter out all reads with zero mapping quality
public class CompareReadGroupWalker extends LocusWalker<Boolean,Boolean>   
{
    @Output
    PrintStream out;
    HashMap<String,Double> readgroups = new HashMap<String,Double>();
    HashMap<String,DescriptiveStatistics> stats = new  HashMap<String,DescriptiveStatistics>();
    String currentContig = "chr0";
   
    public void initialize() 
    {
    	ArrayList<String> rgs = new ArrayList<String>();
    	for (SAMReadGroupRecord s : this.getToolkit().getSAMFileHeader().getReadGroups())
    	{
    		System.err.println("found RG: " + s.getId());
    		rgs.add(s.getId());
    	}  
    	
    	String desc ="";
    	for(int i = 0; i< rgs.size(); i++)
    	{
    		readgroups.put(rgs.get(i), Math.pow(-1,(i+1)));
    		stats.put(rgs.get(i), new DescriptiveStatistics());
    		desc += rgs.get(i) + "=" + readgroups.get(rgs.get(i)) + ", ";
    	}
    	out.println("#type=wiggle_0 description=\"" + desc + "\"");
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
    	HashSet<String> foundReadgroups = new HashSet<String>();
    	
	    for (PileupElement s : context.getBasePileup())
	    {
	    	foundReadgroups.add(s.getRead().getReadGroup().getId());
	    	if(foundReadgroups.size() == readgroups.size())
	    		return true;
	    }
	    
	    if(!context.getContig().equals(currentContig))
    	{
    		out.println("variableStep chrom=" + context.getContig());
    		currentContig = context.getContig();
    	}
    
	    if(foundReadgroups.iterator().hasNext())
	    {
	    	String foundGroup = foundReadgroups.iterator().next();
	    	stats.get(foundGroup).addValue(context.getBasePileup().getNumberOfElements());
	    	out.println(context.getPosition() + "\t" + (readgroups.get(foundGroup) * context.getBasePileup().getNumberOfElements()));
	    }
	    
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
    public Boolean reduce(Boolean value, Boolean val2) 
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
    	for (String key : stats.keySet())
    	{
    		System.out.println(key + "\t exclusively covered loci: " + stats.get(key).getN() + "\t avg coverage on those exclusive loci: " + stats.get(key).getMean());
    	}
    }

	

}
