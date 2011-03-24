package edu.usc.epigenome.uecgatk.benWalkers;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.SAMRecord;

import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgRead;
import edu.usc.epigenome.uecgatk.BaseUtilsMore;

import org.broadinstitute.sting.gatk.filters.MappingQualityReadFilter;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.commandline.Output;

import java.io.PrintStream;
import java.util.Iterator;


/**
 * Translates GATK loci into cytosine objects (along with original GATK data structures). 
 * Keeps a window of cytosines upstream and downstram of current CpG (note that this
 * is not guaranteed to work well with out-of-order sharding strategies.
 */
@ReadFilters( {MappingQualityReadFilter.class} ) // Filter out all reads with zero mapping quality
@Requires( {DataSource.READS, DataSource.REFERENCE, DataSource.REFERENCE_BASES} ) // This walker requires both -I input.bam and -R reference.fasta
public class LocusWalkerToBisulfiteCytosineWalker extends LocusWalker<Integer,Long> {

	/**** GATK Walker implementation ******/
    @Output
    PrintStream out;

	private CpgBackedByGatk prevC = null;
	
	
	/**
	 * 
	 */
	public LocusWalkerToBisulfiteCytosineWalker() {
		super();
		
		// Check sharding strategy
		//this.getToolkit().
	}
	
	

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.Walker#initialize()
	 */
	@Override
	public void initialize() {
		// TODO Auto-generated method stub
		super.initialize();
		
	}
	
	
	
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
    	int centerCoord = thisLoc.getStart();
    	String thisContig = ref.getLocus().getContig();
    	
//		if ( (prevC==null) || !ref.getLocus().onSameContig( prevC.getRefContext().getLocus()))
//    	{
//    		logger.info(String.format("On new contig: %s",thisContig));
//    	}
 
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

       	if (isC)
    	{

       		byte[] contextSeq = 
       			this.getToolkit().getReferenceDataSource().getReference().getSubsequenceAt(thisContig, centerCoord-1, centerCoord+1).getBases();
       		contextSeq = BaseUtilsMore.toUpperCase(contextSeq);
       		byte[] contextSeqStranded = contextSeq;
       		if (negStrand) contextSeqStranded = BaseUtils.simpleReverseComplement(contextSeq);

       		byte[] contextSeqStrandedIupac = new byte[contextSeqStranded.length];
       		for (int i = 0; i < contextSeqStranded.length; i++)
       		{
       			if (i == 1) // Exclude central base, always a C
       			{
       				contextSeqStrandedIupac[i] = contextSeqStranded[i];
       			}
       			else
       			{

       				switch (contextSeqStranded[i])
       				{
       				case BaseUtils.A:
       				case BaseUtils.C:
       				case BaseUtils.T:
       					contextSeqStrandedIupac[i] = (byte)'H';
       					break;
       				default:
       					contextSeqStrandedIupac[i] = contextSeqStranded[i];
       					break;		
       				}
       			}
       		}

 
    	
    	
 
//    		logger.info(String.format("pos=%d\tcontext(%d,%d)=%s (strand=%d)",centerCoord,0,0,new String(windowBases),
//    				(negStrand?-1:1)));
//    		

//    		CpgBackedByGatk thisC = new CpgBackedByGatk(thisLoc.getStart(),negStrand,context, tracker, ref);
//    		logger.info(String.format("Found cytosine: %d: %s", thisC.chromPos, context.getBasePileup().getPileupString(null)));
//
//    		// Increment
//    		prevC = thisC;
    		
    		// Make the Cytosine
    		Cpg thisC = makeCytosine(thisLoc, ref, contextSeqStrandedIupac,negStrand,context);
    		
    		out.printf("%d\t%s\t%s\t%s\t%d\t%s\n", centerCoord,new String(ref.getBases()),
    				new String(contextSeqStranded),new String(contextSeqStrandedIupac),(negStrand?-1:1),thisC.toStringExpanded());


    		// And process it
    		processCytosine(thisC);

    	}
    	
        return 1;
    }

    protected void processCytosine(Cpg thisC) {
		// TODO Auto-generated method stub
		
	}



	/**
	 * @param thisLoc
	 * @param ref
	 * @param contextSeqStrandedIupac This has already been reverse complemented so that middle base is a cytosine
	 * @param cytosineNegStrand
	 * @param context This has all reads relative to the reference genome, so reverse strand cytosines will have to be revcomped
	 * @return
	 */
    protected Cpg makeCytosine(GenomeLoc thisLoc, ReferenceContext ref,
    		byte[] contextSeqStrandedIupac, boolean cytosineNegStrand,
    		AlignmentContext context) 
    {

    	Cpg cOut = new Cpg(thisLoc.getStart(),cytosineNegStrand);
    	short totalReadsOpposite = 0;
    	short aReadOpposite = 0;

    	//**************************************************************
    	//*** These would ideally be set based on the reads rather than the reference
    	boolean nextBaseG = BaseUtils.basesAreEqual(contextSeqStrandedIupac[2],BaseUtils.G);
    	char nextBase = (char)contextSeqStrandedIupac[2];
    	char prevBase = (char)contextSeqStrandedIupac[0];
    	//**************************************************************


    	ReadBackedPileup pileup = context.getBasePileup();
    	Iterator<PileupElement> pileupIt = pileup.iterator();

    	//		int onRead = 0;
    	//		boolean finished = false;
    	//		int curChrPos = -1;
    	while (pileupIt.hasNext())
    	{
    		PileupElement pe = pileupIt.next();

    		byte qual = pe.getQual();

    		// Make sure to handle the strand correctly
    		SAMRecord read = pe.getRead();
    		boolean readOnCytosineStrand = (read.getReadNegativeStrandFlag() == cytosineNegStrand);
    		byte base = pe.getBase();
    		
    		// We change the base to the cytosine-strand
    		byte baseCstrand = (cytosineNegStrand) ? BaseUtils.simpleComplement(base) : base;
    		
    		
    		if (!readOnCytosineStrand)
    		{
        		byte baseGstrand = BaseUtils.simpleComplement(baseCstrand);

//        		out.printf("Got base on NON-C strand: %c\n", (char)baseGstrand);
    			totalReadsOpposite++;
    			if (BaseUtils.basesAreEqual(baseGstrand, BaseUtils.A)) aReadOpposite++;
    		}
    		else
    		{
//        		out.printf("Got base on C strand: %c\n", (char)baseCstrand);

    			boolean isC = BaseUtils.basesAreEqual(baseCstrand, BaseUtils.C);
    			boolean isT = BaseUtils.basesAreEqual(baseCstrand, BaseUtils.T);
    			boolean isAG = BaseUtils.basesAreEqual(baseCstrand, BaseUtils.A) || BaseUtils.basesAreEqual(baseCstrand, BaseUtils.G);

    			//**** THIS IS NOT GUARANTEED TO BE SAFE IF TWO READ NAMES HASH 
    			//**** TO THE SAME INT
    			//**************************************************************
    			String readName = read.getReadName();
    			int readCode = readName.hashCode();
    			//**************************************************************


    			CpgRead cRead = new CpgRead(
    					readCode,
    					(short)(isC?1:0),
    					(short)0,
    					(short)(isT?1:0),
    					(short)((isC||isT)?0:1),
    					(short)(nextBaseG?1:0),
    					nextBase
    			);
    			cOut.addRead(cRead);
//    			out.printf("\tAdding read (BASEQ %d): %s\n", pe.getQual(), cRead.toString());
    		}
    	}

    	// Now finish it off
    	cOut.totalReadsOpposite = totalReadsOpposite;
    	cOut.aReadsOpposite = aReadOpposite;
    	cOut.setNextBaseRef(nextBase);
    	cOut.setPrevBaseRef(prevBase);

    	return cOut;
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