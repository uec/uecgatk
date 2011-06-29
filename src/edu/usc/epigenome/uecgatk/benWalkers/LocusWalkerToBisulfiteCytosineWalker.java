package edu.usc.epigenome.uecgatk.benWalkers;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.SAMRecord;

import edu.usc.epigenome.genomeLibs.PicardUtils;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgRead;
import edu.usc.epigenome.uecgatk.BaseUtilsMore;

import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.gatk.filters.MappingQualityReadFilter;
import org.broadinstitute.sting.gatk.filters.NotPrimaryAlignmentReadFilter;
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
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;


/**
 * Translates GATK loci into cytosine objects (along with original GATK data structures). 
 * Keeps a window of cytosines upstream and downstram of current CpG (note that this
 * is not guaranteed to work well with out-of-order sharding strategies.
 * 
 * Must implement tree-reducible to get parallel execution.
 */
@ReadFilters( {MappingQualityReadFilter.class, BadMateFilter.class, NotPrimaryAlignmentReadFilter.class} ) // Filter out all reads with zero mapping quality
@Requires( {DataSource.READS, DataSource.REFERENCE, DataSource.REFERENCE_BASES} ) // This walker requires both -I input.bam and -R reference.fasta
public abstract class LocusWalkerToBisulfiteCytosineWalker<MapType,ReduceType> extends LocusWalker<MapType,ReduceType> implements TreeReducible<ReduceType> {

    @Argument(fullName = "outputCph", shortName = "cph", doc = "Output CpHs in addition to Cpgs", required = false)
    public boolean outputCph = false;
	
    @Argument(fullName = "minCT", shortName = "mct", doc = "Minimum number of C/T bases to process cytosine (default=1)", required = false)
    public int minCT = 1;
    
    @Argument(fullName = "maxOppAfrac", shortName = "maxa", doc = "Maximum fraction of opposite strand reads being A (default=0.101)", required = false)
    public double maxOppAfrac = 0.101;

    @Argument(fullName = "minNextGfrac", shortName = "minnextg", doc = "Minimum fraction of next position reads being G, only if outputCph is false (default=0.899)", required = false)
    public double minNextGfrac = 0.899;

    @Argument(fullName = "minContextFracReadsMatching", shortName = "contextminmatch", doc = "Minimum fraction of reads matching to call CpH or CpG context (default=0.899)", required = false)
    public double minContextFracReadsMatching = 0.899;

    @Argument(fullName = "minContextPhredQuality", shortName = "contextminphred", doc = "Minimum phred quality for flanking bases to determine context (default=0)", required = false)
    public int minContextPhredQuality = 0;
    
    @Argument(fullName = "minConv", shortName = "minc", doc = "minimum number of converted cytosines required for 5' conversion filter (default=0)", required = false)
    public int minConv = 0;
    
    
    //ZR final static protected String END1_SUFFIX = String.format("%c1", '/');
    //ZR final static protected String END2_SUFFIX = String.format("%c2", '/');
    
    
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
 
		// Increment
		prevContig = thisContig;
		
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

       		byte[] contextSeq = null;
       		try
       		{
       			// I've had this throw an out of range exception from time to time
       			contextSeq = 
       				this.getToolkit().getReferenceDataSource().getReference().getSubsequenceAt(thisContig, centerCoord-1, centerCoord+1).getBases();
       		}
       		catch (Exception e)
       		{
				System.err.printf("Non-fatal error, could not getToolkit().getReferenceDataSource().getReference().getSubsequenceAt(%s,%d,%d)\n%s\n",
						thisContig, centerCoord-1, centerCoord+1,e.toString());
				e.printStackTrace();
				return mapout;
       		}
       		
       		
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
       		Pair<CpgBackedByGatk,MapType> makeCytPair = makeCytosine(thisLoc, ref, contextSeqStrandedIupac,negStrand,context,tracker);
       		CpgBackedByGatk thisC = makeCytPair.first;
       		MapType reduction = makeCytPair.second;
    		
//			out.printf("%d\t%s\t%s\t%s\t%d\t%s\n", centerCoord,new String(ref.getBases()),
//					new String(contextSeqStranded),new String(contextSeqStrandedIupac),(negStrand?-1:1),thisC.toStringExpanded());

       		if (reduction != null)
       		{
       			mapout = reduction;

       		}
       		else if (this.outputCph || !thisC.isCph(false, 0.101))
    		{


    			// And process it
    			boolean minCTpasses = ((thisC.cReads + thisC.tReads) >= this.minCT);
    			double fracOppA = thisC.fracOppositeA();
    			if (Double.isNaN(fracOppA)) fracOppA = 0.0;
    			double fracNextG = thisC.fracNextBaseG();
    			if (Double.isNaN(fracNextG)) fracNextG = 1.0;
    			boolean maxOppAfracPasses = (fracOppA <= maxOppAfrac);
    			boolean minNextGfracPasses = (this.outputCph || (fracNextG>=this.minNextGfrac));
    			
//				logger.info(String.format("Processing cytosine=%s\tCT=%d\toppA=%.2f\n",minCTpasses&&maxOppAfracPasses,(thisC.cReads+thisC.tReads),fracOppA));
   			
    			if (minCTpasses && maxOppAfracPasses && minNextGfracPasses)
    			{
    				mapout = processCytosine(thisC);
    			}
    		}
			
			
    	}
    	


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
	 * @param contextSeqRefStrandedIupac This has already been reverse complemented so that middle base is a cytosine. Reference genome
	 * @param cytosineNegStrand
	 * @param context This has all reads relative to the reference genome, so reverse strand cytosines will have to be revcomped
	 * @return
	 */
    private Pair<CpgBackedByGatk,MapType> makeCytosine(GenomeLoc thisLoc, ReferenceContext ref,
    		byte[] contextSeqRefStrandedIupac, boolean cytosineNegStrand,
    		AlignmentContext context, RefMetaDataTracker tracker) 
    {

    	CpgBackedByGatk cOut = new CpgBackedByGatk(thisLoc.getStart(),cytosineNegStrand, context, tracker, ref);
    	short totalReadsOpposite = 0;
    	short aReadOpposite = 0;
    	
    	List<MapType> baseElementMapList = new ArrayList<MapType>();

    	//**************************************************************
    	//*** These would ideally be set based on the reads rather than the reference
    	boolean nextBaseGref = BaseUtils.basesAreEqual(contextSeqRefStrandedIupac[2],BaseUtils.G);
    	char nextBaseRef = (char)contextSeqRefStrandedIupac[2];
    	char prevBaseRef = (char)contextSeqRefStrandedIupac[0];
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

    		// Make sure to handle the strand correctly. Remember to treat second end reads as the opposite strand from what they actually are!
    		SAMRecord read = pe.getRead();
    		
     		
    		//!read.getFirstOfPairFlag();	 // This is unreliable in early versions of BSMAP, they always make the forward strand read 1
    		String readName = read.getReadName();
    		boolean secondOfPair = getSecondOfPair(read); // EVENTUALLY GET RID OF THIS
//    		boolean secondOfPair = read.getSecondOfPairFlag(); // EVENTUALLY USE THIS
    		
    		
    		// This is super tricky and confusing.  I basically reverse complement the second end for the purposes of methylation analysis.
    		boolean bisulfiteStrandNegativeStrand = read.getReadNegativeStrandFlag();
//    		System.err.printf("secondOfPair=%s\treverseStrand=%s\n", secondOfPair, bisulfiteStrandNegativeStrand);
    		if (secondOfPair) bisulfiteStrandNegativeStrand = !bisulfiteStrandNegativeStrand;
    		boolean readOnCytosineStrand = (bisulfiteStrandNegativeStrand == cytosineNegStrand);
//    		System.err.printf("\tsecondOfPair=%s\tbisulfiteReverseStrand=%s\treadOnCytosineStrand=%s\tcytosineNegStrand=%s\n", secondOfPair, bisulfiteStrandNegativeStrand, readOnCytosineStrand,cytosineNegStrand);    		
    		byte base = pe.getBase();

    		// Next base stuff
    		int peOffset = pe.getOffset();
    		if (read.getReadNegativeStrandFlag()) peOffset = (read.getReadLength()) - peOffset - 1; // nextBaseSeq and prevBaseSeq take the offset relative to the readString 
    	   	byte nextBaseSeqReadStrand = (byte)PicardUtils.nextBaseSeq(read, peOffset, this.minContextPhredQuality);
    	   	byte prevBaseSeqReadStrand = (byte)PicardUtils.preBaseSeq(read, peOffset,  this.minContextPhredQuality);
//    	   	System.err.printf("\tInitial nextSeq=%c, prevSeq=%c (readStart=%d, peOffset=%d)\n", nextBaseSeqReadStrand, prevBaseSeqReadStrand, read.getAlignmentStart(), peOffset);

    		
    		// We don't complement it, because forward and reverse strands aren't strand relative in pileup format. //if (secondOfPair) base = BaseUtils.simpleComplement(base);
    		

    		// We change the bases to the cytosine-strand.  This is complicated because "base" is right now relative to the
    	   	// genome assembly, while nextBaseSeq and prevBaseSeq are relative to the read strand.  So we do them differently,
    	   	// the "base" needs to get adjusted if the cytosine is on the reverse strand, while the next/prev only need
    	   	// to get adjusted when we're on paired-end 2, which means that the bisulfite strand is opposite to the actual 
    	   	// strand of the read.
    	   	byte baseCstrand = base;
    	   	if (cytosineNegStrand)
    	   	{
    			baseCstrand = BaseUtils.simpleComplement(base);
     	   	}

    	   	byte nextBaseSeqCstrand = nextBaseSeqReadStrand;
    	   	byte prevBaseSeqCstrand = prevBaseSeqReadStrand;
    	   	if (secondOfPair)
    	   	{
    	   		// Notice that these will generally not be used, because context is only relevant on the C strand.
       			nextBaseSeqCstrand = BaseUtils.simpleComplement(prevBaseSeqReadStrand); // Note that we swap with prev Base to change orientation
    			prevBaseSeqCstrand = BaseUtils.simpleComplement(nextBaseSeqReadStrand); // Note that we swap with prev Base to change orientation
//        	   	System.err.printf("\tCytosine rev strand, now nextSeq=%c, prevSeq=%c\n", nextBaseSeqCstrand, prevBaseSeqCstrand);
    	   	}
    	   	boolean nextBaseSeqG = BaseUtils.basesAreEqual(nextBaseSeqCstrand,BaseUtils.G); 
    	   	boolean prevBaseSeqG = BaseUtils.basesAreEqual(prevBaseSeqCstrand,BaseUtils.G); 


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
    			
    			int cycle = pe.getOffset();
    			int readCycleQual = pe.getQual();
    			
    			MapType mapElement = mapReferenceCytosineBase(thisLoc, cycle, qual, baseCstrand, prevBaseRef, nextBaseRef); 
    			if (mapElement != null)
    			{
    				baseElementMapList.add(mapElement);
    			}
    			

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
    					(short)(nextBaseSeqG?1:0),
    					(char)nextBaseSeqCstrand,
    					(short)(prevBaseSeqG?1:0),
    					(char)prevBaseSeqCstrand
    			);
    			cOut.addRead(cRead);
//    			System.err.printf("\tAdding read (%d, BASEQ %d): %s\n", thisLoc.getStart(), pe.getQual(), cRead.toString());
    		}
    	}

    	// Now finish it off
    	cOut.totalReadsOpposite = totalReadsOpposite;
    	cOut.aReadsOpposite = aReadOpposite;
    	cOut.setNextBaseRef(nextBaseRef);
    	cOut.setPrevBaseRef(prevBaseRef);

//		out.printf("Adding cytosine: %s\n", cOut.toString());

    	ReduceType outReduce = null;
    	for (MapType map : baseElementMapList)
    	{
    		if (outReduce == null) outReduce = this.reduceInit();
    		outReduce = this.reduce(map, outReduce);

    	}
    	
    	MapType outMap = null;
    	try
    	{
    		outMap = (MapType)outReduce;
    	}
    	catch (Exception e)
    	{
    		System.err.printf("Fatal error: If you use LocusWalkerToBisulfiteCytosineWalker::" + 
    				"mapReferenceCytosineBase, ReduceType must be castable to MapType\n%s",e.toString());
    		e.printStackTrace();
    		System.exit(1);

    	}
    	
    	return new Pair<CpgBackedByGatk,MapType>(cOut,outMap);
    }



	protected MapType mapReferenceCytosineBase(GenomeLoc thisLoc,
			int cycle, int qual, byte baseCstrand, char prevBaseRef, char nextBaseRef) {
		// Doesn't do anything.  Override for counters
		MapType out = null;
		return out;
	}



	protected static boolean getSecondOfPair(SAMRecord read) {
		boolean out = false;
		
		if (read.getReadPairedFlag())
		{
			out = read.getSecondOfPairFlag(); 
		}
		return out;
		
//		boolean secondOfPair = false;
//		String readName = read.getReadName();
//
//		if (read.getReadPairedFlag())
//		{
//			if (readName.endsWith(END1_SUFFIX))
//			{
//				secondOfPair = false;
//			}
//			else if (readName.endsWith(END2_SUFFIX))
//			{
//				secondOfPair = true;   			
//			}
//			else
//			{
//				System.err.println("LocusWalkerToBisulfiteCytosineWalker::getSecondOfPair() Got a read that doesn't end with /1 or /2: " + readName + ".  Can't tell which end it is.");
//				System.exit(1);
//			}	
//		}
//		
//		return secondOfPair;
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
	    		byte[] refBases = null;
	    		try
	    		{
	    			refBases = BaseUtilsMore.toUpperCase(this.getToolkit().getReferenceDataSource().getReference().getSubsequenceAt(thisContig, contigCoord, contigCoord).getBases());
	    		}
	      		catch (Exception e)
	       		{
					System.err.printf("Non-fatal error, could not getToolkit()getReferenceDataSource().getReference().getSubsequenceAt(%s,%d,%d)\n%s\n",
							thisContig, contigCoord, contigCoord,e.toString());
					e.printStackTrace();
					return passes;
	       		}
	      		
	    		byte refBase = (revStrand) ? BaseUtils.simpleComplement(refBases[0]) : refBases[0];
	    		byte readBase = readBases[offset];

	    		// Conversion filter is still relative to the 5' end relative to the read, whether it's the
	    		// bisulfite or bisulfite-rc strand
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
	    		
    		//logger.info(String.format("\t\tChecking ref(%d),index(%d) = %c,%c (%s)",contigCoord,offset,refBase,readBase, conv));
	    		
	    		
	    	}

    		if (numConv>=this.minConv) passes = true;
	    	
	    	
//	    	int lastCoordDiff = Math.abs(contigCoord-currentLoc); // Debugging
//	    	if (lastCoordDiff>0) logger.info(String.format("\t%s -- Passes(%d) got pe: %c (%d) diff=%d [read %s:%d-%d]",  passes, currentLoc, pe.getBase(), thisOffset, lastCoordDiff, thisContig, read.getAlignmentStart(), read.getAlignmentEnd()));			
			

		}

		return passes;
	}






}