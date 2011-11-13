package edu.usc.epigenome.uecgatk.benWalkers.cytosineWalkers;

import java.io.File;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;

import org.biojava.bio.seq.StrandedFeature;
import org.broadinstitute.sting.commandline.Argument;
import org.usckeck.genome.ChromFeatures;

import edu.usc.epigenome.uecgatk.benWalkers.CpgBackedByGatk;
import edu.usc.epigenome.uecgatk.benWalkers.LocusWalkerToBisulfiteCytosineWalker;


/**
 * @author benb
 * 
 */
public class BisulfiteSeqToBareWigWalker extends LocusWalkerToBisulfiteCytosineWalker<Integer, Long> {


    @Argument(fullName = "outPrefix", shortName = "pre", doc = "Output prefix for all output files", required = true)
    public String outPrefix = null;
 
    @Argument(fullName = "gnomeMode", shortName = "gnome", doc = "GNOMe mode - outputs only GCH datapoints (default false)", required = false)
    public boolean gnomeMode = false;
    
    @Argument(fullName = "outputAllContexts", shortName = "all", doc = "Output optional contexts HCH, GCG (default false)", required = false)
    public boolean outputAllContexts = false;
    
    @Argument(fullName = "bedMode", shortName = "bed", doc = "Output with chr number field, not guaranteed to be ordered, but can run in multi-thread mode", required = false)
    public boolean bedMode = false;
   
    @Argument(fullName = "csvMode", shortName = "csv", doc = "Output with chr number field, not guaranteed to be ordered, but can run in multi-thread mode", required = false)
    public boolean csvMode = false;
   
    Map<String,PrintWriter> wigByContext = new HashMap<String,PrintWriter>();
    
    int curChrNum = 0;
    String curChr = "";
    
	/**
	 * locus walker overrides
	 */


	/**
	 * Provides an initial value for the reduce function.  Hello walker counts loci,
	 * so the base case for the inductive step is 0, indicating that the walker has seen 0 loci.
	 * @return 0.
	 */
	@Override
	public Long reduceInit()
	{ 
		Long out = 0L;
		return out;
	}

	@Override
	public void initialize() {
		super.initialize();

		if (!bedMode && !csvMode)
		{	
			if (this.getToolkit().getArguments().numberOfThreads>1)
			{
				System.err.println("GnomeSeqToBareWigWalker does not yet implement multi-threaded mode. Use -nt 1");
				System.exit(1);
			}
		}
		
//		this.outputCph = true; // Because GNOME-seq used GCH
	}


	

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.TreeReducible#treeReduce(java.lang.Object, java.lang.Object)
	 * 
	 * Must override this to get multithreaded capability
	 */
	@Override
	public Long treeReduce(Long lhs, Long rhs) 
	{
		if (!bedMode && !csvMode)
		{	
			System.err.println("BisulfiteSeqToBareWigWalker does not yet implement multi-threaded mode. Use -nt 1");
			System.exit(1);
		}
		System.err.printf("treeReduce(%d,%d)\n",lhs,rhs);
		return lhs + rhs;
	}

	/**
	 * Retrieves the final result of the traversal.
	 * @param result The ultimate value of the traversal, produced when map[n] is combined with reduce[n-1]
	 *               by the reduce function. 
	 */
	@Override
	public void onTraversalDone(Long result) 
	{
		out.printf("Saw %d cytosines\n", result);
		for (String context : this.wigByContext.keySet())
		{
			PrintWriter pw = this.wigByContext.get(context);
			pw.close();
		}
	}

	
	/***************************************************
	 * cytosine walker overrides
	 ***************************************************/
	
	
	@Override
	protected void alertNewContig(String newContig) 
	{
		System.err.printf("New chrom: %s\n",newContig);
		this.curChr = newContig;
		this.curChrNum = (new ChromFeatures()).chrom_from_public_str(newContig);
		
		for (String context : this.wigByContext.keySet())
		{
			PrintWriter pw = this.wigByContext.get(context);
			pw.printf("variableStep chrom=%s\n", newContig);
		}
	}

	

	@Override
	protected Integer processCytosine(CpgBackedByGatk thisC)
	{
		Integer out = 0;
		//writeCoverage(thisC);
		
		String context = thisC.context();
		
		if (!this.gnomeMode)
		{
			context = context.substring(1);
		}

		// Skip optional contexts
		if (!this.outputAllContexts)
		{
			if (context.equalsIgnoreCase("HCH") || context.equalsIgnoreCase("GCG")) return out;
		}
		
		double meth = thisC.fracMeth(false);
		int nreads = (int)thisC.totalReadsCorT(true);
		
		PrintWriter wigpw = null;
		if (this.wigByContext.containsKey(context))
		{
			wigpw = this.wigByContext.get(context);
			//logger.info("Found wig " + wig.toString());
		}
		else
		{
			// We use the wigByContext for merging
//			String outfn = String.format("%s.%s.wig.%s",this.outPrefix,context,this.wigByContext.hashCode());
			String name = String.format("%s.%s-minct%d-minconv%d", this.outPrefix, context, this.minCT, this.minConv);
			String suffix = (this.bedMode) ? ".bed" : ((this.csvMode) ? ".csv" : ".wig"); 
			String outfn = String.format("%s%s",name,suffix);
			logger.info("NEW wig " + outfn);
			try
			{
				wigpw = new PrintWriter(new File(outfn));
			}
			catch (Exception e)
			{
				System.err.printf("Fatal error, BisulfiteSeqToBareWigWalker can't create file: %s\n",outfn);
				e.printStackTrace();
				System.exit(1);
			}
			
			if (!csvMode)
			{
				String desc = String.format("%s-wig", name);
				wigpw.printf("track type=wiggle_0 name=%s description=%s color=204,102,0 visibility=full " +
						" graphType=points autoScale=off alwaysZero=off maxHeightPixels=64:32:10 viewLimits=0:100\n", desc, desc);	
				wigpw.printf("variableStep chrom=%s\n", this.curChr);
			}

			this.wigByContext.put(context, wigpw);
		}
		
		if (this.bedMode)
		{
			wigpw.printf("%s\t%d\t%d\t%d\n",this.curChr, thisC.getGenomeLoc().getStart(),
					thisC.getGenomeLoc().getStart(), (int)Math.round(100.0*meth));
		}
		else if (this.csvMode)
		{
			int strandint = (thisC.getStrand() == StrandedFeature.NEGATIVE) ? -1 : 1; // Unstranded does not exist for Cpg objects.
			wigpw.printf("%d,%d,%d,%d,%d\n",this.curChrNum, thisC.getGenomeLoc().getStart(), strandint, 
					(int)Math.round(100.0*meth), nreads);
		}
		else
		{
			wigpw.printf("%d\t%d\n",thisC.getGenomeLoc().getStart(), (int)Math.round(100.0*meth));
		}
		//logger.info(String.format("\tWriting pos %d to wig \"%s\" (this=%s)\n", thisC.chromPos, context,this));
		out++;
		
		return out;
	}



	@Override
	protected Long reduceCytosines(Integer value, Long sum) 
	{
		return value.longValue() + sum;
	}

	
	//////////////////////
	/// Private funcs
	//////////////////////

//	private void writeCoverage(CpgBackedByGatk thisC)
//	{
//		String context = "readcvg";
//		WiggleWriterReducible wig = null;
//		if (this.wigByContext.containsKey(context))
//		{
//			wig = this.wigByContext.get(context);
//			//logger.info("Found wig " + wig.toString());
//		}
//		else
//		{
//			// We use the wigByContext for merging
////			String outfn = String.format("%s.%s.wig.%s",this.outPrefix,context,this.wigByContext.hashCode());
//			String name = String.format("%s.%s", this.outPrefix, context);
//			String outfn = String.format("%s.wig",name,context);
//			logger.info("NEW wig " + outfn);
//			wig = new WiggleWriterReducible(new File(outfn));
//			
//			WiggleHeader head = new WiggleHeaderCytosines(name, name);
//			wig.writeHeader(head);
//			this.wigByContext.put(context, wig);
//		}
//		
//		wig.writeData(thisC.getGenomeLoc(), thisC.totalReads);
//	}

}



