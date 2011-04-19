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
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;

import java.io.PrintStream;
import java.util.Iterator;


/**
 * Translates GATK loci into cytosine objects (along with original GATK data structures). 
 * Keeps a window of cytosines upstream and downstram of current CpG (note that this
 * is not guaranteed to work well with out-of-order sharding strategies.
 * 
 * Must implement tree-reducible to get parallel execution.
 */
@ReadFilters( {MappingQualityReadFilter.class} ) // Filter out all reads with zero mapping quality
@Requires( {DataSource.READS, DataSource.REFERENCE, DataSource.REFERENCE_BASES} ) // This walker requires both -I input.bam and -R reference.fasta
public abstract class LocusWalkerToBisulfiteCytosineWalker<MapType,ReduceType> extends LocusWalker<MapType,ReduceType> implements TreeReducible<ReduceType> {

    @Argument(fullName = "outputCph", shortName = "cph", doc = "Output CpHs in addition to Cpgs", required = false)
    public boolean outputCph = false;
	
    @Argument(fullName = "minCT", shortName = "mct", doc = "Minimum number of C/T bases to process cytosine (default=1)", required = false)
    public int minCT = 1;
    
    @Argument(fullName = "maxOppAfrac", shortName = "maxa", doc = "Maximum fraction of opposite strand reads being A (default=0.101)", required = false)
    public double maxOppAfrac = 0.101;

    @Argument(fullName = "minConv", shortName = "minc", doc = "minimum number of converted cytosines required for 5' conversion filter (default=0)", required = false)
    public int minConv = 0;
    
    
    final static protected String END1_SUFFIX = String.format("%c1", '/');
    final static protected String END2_SUFFIX = String.format("%c2", '/');
    
    
	/**** GATK Walker implementation ******/
    @Output
    protected PrintStream out;

	protected String prevContig = null;
	
	
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
    public MapType map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

    	MapType mapout = null;
    	
    	
    	// Are we on a new chrom?
    	GenomeLoc thisLoc = ref.getLocus();
    	int centerCoord = thisLoc.getStart();
    	String thisContig = ref.getLocus().getContig();
    	
		if ( (prevContig==null) || !thisContig.equalsIgnoreCase(prevContig) )
    	{
    		logger.info(String.format("On new contig: (%s, %s, %s)",prevContig,thisContig,this));
    		this.alertNewContig(thisContig);
    	}
 
		// *** NOTE :: JUST USING POSITIONS THAT ARE C IN REFERENCE :: NEED TO CHECK SAMPLE GENOME ****
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
    	
 
//    		logger.info(String.format("pos=%d\tcontext(%d,%d)=%s (strand=%d)",centerCoord,0,0,new String(contextSeqStrandedIupac),
//    				(negStrand?-1:1)));
    		
    		// Make the Cytosine
       		CpgBackedByGatk thisC = makeCytosine(thisLoc, ref, contextSeqStrandedIupac,negStrand,context,tracker);
    		
//			out.printf("%d\t%s\t%s\t%s\t%d\t%s\n", centerCoord,new String(ref.getBases()),
//					new String(contextSeqStranded),new String(contextSeqStrandedIupac),(negStrand?-1:1),thisC.toStringExpanded());

			if (this.outputCph || !thisC.isCph(false, 0.101))
    		{


    			// And process it
    			boolean minCTpasses = ((thisC.cReads + thisC.tReads) >= this.minCT);
    			double fracOppA = thisC.fracOppositeA();
    			if (Double.isNaN(fracOppA)) fracOppA = 0.0;
    			boolean maxOppAfracPasses = (fracOppA <= maxOppAfrac);
//				logger.info(String.format("Processing cytosine=%s\tCT=%d\toppA=%.2f\n",minCTpasses&&maxOppAfracPasses,(thisC.cReads+thisC.tReads),fracOppA));
   			
    			if (minCTpasses & maxOppAfracPasses)
    			{
    				mapout = processCytosine(thisC);
    			}
    		}
			
			
    	}
    	
		// Increment
		prevContig = thisContig;

		return mapout;
    }

   abstract protected void alertNewContig(String newContig);
   abstract protected MapType processCytosine(CpgBackedByGatk thisC);

   
	/**
	 * Combines the result of the latest map with the accumulator.  In inductive terms,
	 * this represents the step loci[x + 1] = loci[x] + 1
	 * @param value result of the map.
	 * @param sum accumulator for the reduce.
	 * @return The total count of loci processed so far.
	 * 
	 * We implement this to ignore non-cytosines.  Subclasses should implement reduceCytosines
	 */
	@Override
	public ReduceType reduce(MapType value, ReduceType sum) {
		ReduceType out = sum;
		if (value != null)
		{
			out = reduceCytosines(value, sum);
		}
		return out;
	}

	abstract protected ReduceType reduceCytosines(MapType value, ReduceType sum);



	/**
	 * @param thisLoc
	 * @param ref
	 * @param contextSeqStrandedIupac This has already been reverse complemented so that middle base is a cytosine
	 * @param cytosineNegStrand
	 * @param context This has all reads relative to the reference genome, so reverse strand cytosines will have to be revcomped
	 * @return
	 */
    protected CpgBackedByGatk makeCytosine(GenomeLoc thisLoc, ReferenceContext ref,
    		byte[] contextSeqStrandedIupac, boolean cytosineNegStrand,
    		AlignmentContext context, RefMetaDataTracker tracker) 
    {

    	CpgBackedByGatk cOut = new CpgBackedByGatk(thisLoc.getStart(),cytosineNegStrand, context, tracker, ref);
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

    		// Make sure to handle the strand correctly. Remember to treat second end reads as the opposite strand!
    		SAMRecord read = pe.getRead();
     		
    		//!read.getFirstOfPairFlag();	 // This is unreliable in BSMAP, they always make the forward strand read 1
    		String readName = read.getReadName();
    		boolean secondOfPair = getSecondOfPair(read);
    		
    		boolean readStrandNegativeStrand = read.getReadNegativeStrandFlag();
//    		out.printf("secondOfPair=%s\treverseStrand=%s\n", secondOfPair, readStrandNegativeStrand);
    		if (secondOfPair) readStrandNegativeStrand = !readStrandNegativeStrand;
    		boolean readOnCytosineStrand = (readStrandNegativeStrand == cytosineNegStrand);
//    		out.printf("\tsecondOfPair=%s\tadjustedReverseStrand=%s\treadOnCytosineStrand=%s\n", secondOfPair, readStrandNegativeStrand, readOnCytosineStrand);    		
    		byte base = pe.getBase();
    		
    		
    		// We don't complement it, because forward and reverse strands aren't strand relative in pileup format. //if (secondOfPair) base = BaseUtils.simpleComplement(base);
    		

    		// We change the base to the cytosine-strand
    		byte baseCstrand = (cytosineNegStrand) ? BaseUtils.simpleComplement(base) : base;


    		if (!readOnCytosineStrand)
    		{
    			byte baseGstrand = BaseUtils.simpleComplement(baseCstrand);

//    			out.printf("Got base on NON-C strand: %c\n", (char)baseGstrand);
    			totalReadsOpposite++;
    			if (BaseUtils.basesAreEqual(baseGstrand, BaseUtils.A)) aReadOpposite++;
    		}
    		else
    		{
    			boolean passesFiveprimeFilter = PassesFiveprimeFilter(pe, ref);

    			//    				out.printf("\t\tGot base on C strand: %c\n", (char)baseCstrand);

    			boolean isC = BaseUtils.basesAreEqual(baseCstrand, BaseUtils.C);
    			boolean isT = BaseUtils.basesAreEqual(baseCstrand, BaseUtils.T);
    			boolean isAG = BaseUtils.basesAreEqual(baseCstrand, BaseUtils.A) || BaseUtils.basesAreEqual(baseCstrand, BaseUtils.G);

    			//**** THIS IS NOT GUARANTEED TO BE SAFE IF TWO READ NAMES HASH 
    			//**** TO THE SAME INT
    			//**************************************************************
    			int readCode = readName.hashCode();
    			//**************************************************************


    			CpgRead cRead = new CpgRead(
    					readCode,
    					(short)((isC&passesFiveprimeFilter)?1:0),
    					(short)((isC&!passesFiveprimeFilter)?1:0),
    					(short)(isT?1:0),
    					(short)((isC||isT)?0:1),
    					(short)(nextBaseG?1:0),
    					nextBase
    			);
    			cOut.addRead(cRead);
 //   			out.printf("\tAdding read (BASEQ %d): %s\n", pe.getQual(), cRead.toString());
    		}
    	}

    	// Now finish it off
    	cOut.totalReadsOpposite = totalReadsOpposite;
    	cOut.aReadsOpposite = aReadOpposite;
    	cOut.setNextBaseRef(nextBase);
    	cOut.setPrevBaseRef(prevBase);

//		out.printf("Adding cytosine: %s\n", cOut.toString());

    	return cOut;
    }



	protected static boolean getSecondOfPair(SAMRecord read) {
		boolean secondOfPair = false;
		String readName = read.getReadName();
		if (readName.endsWith(END1_SUFFIX))
		{
			secondOfPair = false;
		}
		else if (readName.endsWith(END2_SUFFIX))
		{
			secondOfPair = true;   			
		}
		else
		{
			System.err.println("Got a read that doesn't end with /1 or /2: " + readName + ".  Can't tell which end it is.");
			System.exit(1);
		}	
		
		return secondOfPair;
	}



	/**
	 * ******* You should only pass in PileupElements on the cytosine strand. *******
	 * 
	 * @param pe
	 * @param ref
	 * @return
	 */
	protected boolean PassesFiveprimeFilter(PileupElement pe, ReferenceContext ref) 
	{
		boolean passes = false;
		if (this.minConv == 0)
		{
			passes = true;
		}
		else
		{
			SAMRecord read = pe.getRead();
			int readLen = read.getReadLength();
	    	String thisContig = ref.getLocus().getContig();
	    	boolean revStrand = read.getReadNegativeStrandFlag();
	    	
	    	byte[] readBases = read.getReadBases();
			if (revStrand) readBases = BaseUtils.simpleReverseComplement(readBases);

	    	int thisOffset = (revStrand) ? ((readLen-1)-pe.getOffset()) : pe.getOffset(); // It's unexpected that Offset is relative to genomic coords
	    	int currentLoc = ref.getLocus().getStart();
	    	
	    	int readStartPos = (revStrand) ? read.getAlignmentEnd() : read.getAlignmentStart();
	    	int increment = (revStrand) ? -1 : 1;

	    	int contigCoord = 0;
	    	int numConv = 0;
    		boolean secondEnd = getSecondOfPair(read);
	    	for (int offset = 0; (offset<=thisOffset)&& (numConv<this.minConv) ; offset++) // 
	    	{
	    		contigCoord = readStartPos + (increment * offset);
	    		byte[] refBases = BaseUtilsMore.toUpperCase(this.getToolkit().getReferenceDataSource().getReference().getSubsequenceAt(thisContig, contigCoord, contigCoord).getBases());
	    		byte refBase = (revStrand) ? BaseUtils.simpleComplement(refBases[0]) : refBases[0];
	    		byte readBase = readBases[offset];
	    		
	    		boolean conv = false;
	    		if (secondEnd)
	    		{
	    			conv = ((refBase==BaseUtils.G) && (readBase==BaseUtils.A));
	    		}
	    		else
	    		{
	    			conv = ((refBase==BaseUtils.C) && (readBase==BaseUtils.T));
	    		}
	    		
	    		
	    		if (conv) numConv++;
	    		
    		// logger.info(String.format("\t\tChecking ref(%d),index(%d) = %c,%c (%s)",contigCoord,offset,refBase,readBase, conv));
	    		
	    		
	    	}

    		if (numConv>=this.minConv) passes = true;
	    	
	    	
//	    	int lastCoordDiff = Math.abs(contigCoord-currentLoc); // Debugging
//	    	if (lastCoordDiff>0) logger.info(String.format("\t%s -- Passes(%d) got pe: %c (%d) diff=%d [read %s:%d-%d]",  passes, currentLoc, pe.getBase(), thisOffset, lastCoordDiff, thisContig, read.getAlignmentStart(), read.getAlignmentEnd()));			
			

		}

		return passes;
	}






}