package edu.usc.epigenome.uecgatk.benWalkers;

import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.commandline.Output;

import java.io.PrintStream;

/**
 * An example walker for illustrative purposes.
 */
public class HelloWalker extends LocusWalker<Integer,Long> {
	
	static int counter = 0;
	
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
    	
    	if ((counter % 50000)==0)
    	{
    		out.printf("Hello %d  locus %s; your ref base is %c and you have %d reads%n", counter, context.getLocation(), ref.getBase(), context.getBasePileup().getBases().length );
    	}
        counter++;
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