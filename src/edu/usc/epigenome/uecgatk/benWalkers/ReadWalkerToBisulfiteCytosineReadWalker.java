package edu.usc.epigenome.uecgatk.benWalkers;

import net.sf.samtools.SAMRecord;

import edu.usc.epigenome.genomeLibs.PicardUtils;

import org.broadinstitute.sting.gatk.filters.MappingQualityReadFilter;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;

import java.io.PrintStream;


/**
 * Translates GATK loci into cytosine objects (along with original GATK data structures). 
 * Keeps a window of cytosines upstream and downstram of current CpG (note that this
 * is not guaranteed to work well with out-of-order sharding strategies.
 * 
 * Must implement tree-reducible to get parallel execution.
 * 
 * We only have one type here, because we have to do internal reduce steps
 */
@ReadFilters( {MappingQualityReadFilter.class} ) // Filter out all reads with zero mapping quality
@Requires( {DataSource.READS, DataSource.REFERENCE, DataSource.REFERENCE_BASES} ) // This walker requires both -I input.bam and -R reference.fasta
public abstract class ReadWalkerToBisulfiteCytosineReadWalker<ReduceType> extends ReadWalker<ReduceType,ReduceType> implements TreeReducible<ReduceType> {

    @Argument(fullName = "outputCph", shortName = "cph", doc = "Output CpHs in addition to Cpgs", required = false)
    public boolean outputCph = true;
	

	/**** GATK Walker implementation ******/
    @Output
    protected PrintStream out;

	protected String prevContig = null;
	
	
	/**
	 * 
	 */
	public ReadWalkerToBisulfiteCytosineReadWalker() {
		super();
		
		// Check sharding strategy
		//this.getToolkit().
	}
	
	

	abstract protected void alertNewContig(String newContig);
	abstract protected ReduceType processReadCytosine(int positionInRead, String cContext, boolean isMethylated);

//	/* (non-Javadoc)
//	 * @see org.broadinstitute.sting.gatk.walkers.Walker#initialize()
//	 */
//	@Override
//	public void initialize() {
//		super.initialize();
//		
//	}
//	
	
	   
//	@Override
//	public ReduceType treeReduce(ReduceType lhs, ReduceType rhs) {
//		return null;
//	}
//
//
//
//	@Override
//	public boolean filter(ReferenceContext ref, SAMRecord read) {
//		return super.filter(ref, read);
//	}



	@Override
	public ReduceType map(ReferenceContext ref, SAMRecord read,
			ReadMetaDataTracker metaDataTracker) 
	{
    	GenomeLoc thisLoc = ref.getLocus();
    	String thisContig = thisLoc.getContig();
		if ( (prevContig==null) || !thisContig.equalsIgnoreCase(prevContig) )
    	{
    		logger.info(String.format("On new contig: (%s, %s, %s)",prevContig,thisContig,this));
    		this.alertNewContig(thisContig);
    	}

		int mapQual = read.getMappingQuality();
		boolean unmapped = read.getReadUnmappedFlag();
    	boolean revStrand = read.getReadNegativeStrandFlag();

    	// Get the sequences
		String seq = PicardUtils.getReadString(read, true);
		byte[] readSeq = read.getReadBases();
    	byte[] refSeq = ref.getBases();
    	if (revStrand)
    	{
    		refSeq = BaseUtils.simpleReverseComplement(refSeq);
    		readSeq = BaseUtils.simpleReverseComplement(readSeq);
    	}
    	
    	
    	// Don't do first and last one because they don't have context.
		ReduceType out = this.reduceInit();
    	for (int i = 1; i < (readSeq.length-1); i++)
    	{
    		byte refBase = refSeq[i];
    		byte readBase = readSeq[i];
    		
			boolean isC = BaseUtils.basesAreEqual(refBase, BaseUtils.C) && 
			(BaseUtils.basesAreEqual(readBase, BaseUtils.C) || BaseUtils.basesAreEqual(readBase, BaseUtils.T));
			
			if (isC)
			{
				if (this.outputCph || BaseUtils.basesAreEqual(readSeq[i+1], BaseUtils.G))
				{
					String cContext = getCytosineContext(i, readSeq, refSeq);
					boolean isMethylated = BaseUtils.basesAreEqual(readBase, BaseUtils.C);
					ReduceType thisOut = processReadCytosine(i, cContext, isMethylated);
					
					out = this.reduce(out, thisOut);
				}
			}
    	}
    	
   // 	logger.info(String.format("Got read:\n\tref =%s\n\tread=%s", new String(refSeq), new String(readSeq)))
		
		// Increment
		prevContig = thisContig;
		
		return out;
	}    	


	private String getCytosineContext(int offset, byte[] readSeq, byte[] refSeq) 
	{
   		byte[] contextSeqIupac = new byte[3];
   		for (int i = (offset-1); i <= (offset+1); i++)
   		{
   			int contextPos = (i-offset+1);
   			if (i == offset) // Exclude central base, always a C
   			{
   				contextSeqIupac[contextPos] = refSeq[i]; 
   			}
   			else
   			{

   				switch (readSeq[i]) // Ok to use bisulfite read since C/T go to same IUPAC code
   				{
   				case BaseUtils.A:
   				case BaseUtils.C:
   				case BaseUtils.T:
   					contextSeqIupac[contextPos] = (byte)'H';
   					break;
   				default:
   					contextSeqIupac[contextPos] = readSeq[i];
   					break;		
   				}
   			}
   		}
   		return new String(contextSeqIupac);
	}






}