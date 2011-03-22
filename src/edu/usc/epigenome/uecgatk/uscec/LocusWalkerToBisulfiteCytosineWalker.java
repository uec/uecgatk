package edu.usc.epigenome.uecgatk.uscec;

import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.commandline.Output;

import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;

import java.io.PrintStream;
import java.util.Iterator;

/**
 * Translates GATK loci into cytosine objects (along with original GATK data structures). 
 * Keeps a window of cytosines upstream and downstram of current CpG (note that this
 * is not guaranteed to work well with out-of-order sharding strategies.
 */
public class LocusWalkerToBisulfiteCytosineWalker extends LocusWalker<Integer,Long> implements Iterator<Cpg> {

	
	private CpgBackedByGatk prevC = null;
	
	
	/**
	 * 
	 */
	public LocusWalkerToBisulfiteCytosineWalker() {
		super();
		
		// Check sharding strategy
		//this.getToolkit().
	}

	/*** Cpg iterator implementation *****/

	public boolean hasNext() {
		// TODO Auto-generated method stub
		return false;
	}

	public Cpg next() {
		// TODO Auto-generated method stub
		return null;
	}

	public void remove() {
		// TODO Auto-generated method stub
		
	}	
	
	
	
	/**** GATK Walker implementation ******/
    @Output
    PrintStream out;

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
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

    	// Are we on a new chrom?
    	GenomeLoc thisLoc = ref.getLocus();
    	
		if ( (prevC==null) || !ref.getLocus().onSameContig( prevC.getRefContext().getLocus()))
    	{
    		logger.info(String.format("On new contig: %s",ref.getLocus().getContig()));
    	}
    	
    	boolean isC = false;
    	boolean negStrand = false;
    	if (ref.getBase() == BaseUtils.C)
    	{
    		isC = true;
    		negStrand = false;
    	}
    	else if (ref.getBase() == BaseUtils.G)
    	{
    		isC = true;
    		negStrand = true;
    	}
    	// Check strand.
    	
    	
    	if (isC)
    	{
    		CpgBackedByGatk thisC = new CpgBackedByGatk(thisLoc.getStart(),negStrand,context, tracker, ref);
    		logger.info(String.format("Found cytosine: %d: %s", thisC.chromPos, context.getBasePileup().getPileupString(null)));

    		// Increment
    		prevC = thisC;

    	}
    	
        return 1;
    }

    /**
     * Provides an initial value for the reduce function.  Hello walker counts loci,
     * so the base case for the inductive step is 0, indicating that the walker has seen 0 loci.
     * @return 0.
     */
    @Override
    public Long reduceInit() { return 0L; }

    /**
     * Combines the result of the latest map with the accumulator.  In inductive terms,
     * this represents the step loci[x + 1] = loci[x] + 1
     * @param value result of the map.
     * @param sum accumulator for the reduce.
     * @return The total count of loci processed so far.
     */
    @Override
    public Long reduce(Integer value, Long sum) {
        return sum + value;
    }

    /**
     * Retrieves the final result of the traversal.
     * @param result The ultimate value of the traversal, produced when map[n] is combined with reduce[n-1]
     *               by the reduce function. 
     */
    @Override
    public void onTraversalDone(Long result) {
        out.println("Number of loci viewed is: " + result);
    }

}