package edu.usc.epigenome.uecgatk.benWalkers.cytosineWalkers;

import java.io.File;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;

import org.broadinstitute.sting.commandline.Argument;

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
    
    
    Map<String,PrintWriter> wigByContext = new HashMap<String,PrintWriter>();
    
    
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
		
		if (this.getToolkit().getArguments().numberOfThreads>1)
		{
			System.err.println("GnomeSeqToBareWigWalker does not yet implement multi-threaded mode. Use -nt 1");
			System.exit(1);
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
		System.err.println("GnomeSeqToBareWigWalker does not yet implement multi-threaded mode. Use -nt 1");
		System.exit(1);
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
			String outfn = String.format("%s.wig",name,context);
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
			
			String desc = String.format("%s-wig", name);
			wigpw.printf("track type=wiggle_0 name=%s description=%s color=204,102,0 visibility=full " +
					" graphType=points autoScale=off alwaysZero=off maxHeightPixels=64:32:10 viewLimits=0:100\n", desc, desc);	
			wigpw.printf("variableStep chrom=%s\n", thisC.getGenomeLoc().getContig());

			this.wigByContext.put(context, wigpw);
		}
		
		wigpw.printf("%d\t%d\n",thisC.getGenomeLoc().getStart(), (int)Math.round(100.0*meth));
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



