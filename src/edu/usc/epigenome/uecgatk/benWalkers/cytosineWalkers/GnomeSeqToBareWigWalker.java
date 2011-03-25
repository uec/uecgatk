package edu.usc.epigenome.uecgatk.benWalkers.cytosineWalkers;

import java.io.File;
import java.io.PrintStream;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.utils.wiggle.WiggleHeader;

import edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers.CpgMethLevelSummarizer;
import edu.usc.epigenome.uecgatk.benWalkers.CpgBackedByGatk;
import edu.usc.epigenome.uecgatk.benWalkers.LocusWalkerToBisulfiteCytosineWalker;
import edu.usc.epigenome.uecgatk.benWalkers.WiggleHeaderCytosines;
import edu.usc.epigenome.uecgatk.WiggleWriterReducible;

/**
 * @author benb
 * 
 */
public class GnomeSeqToBareWigWalker extends LocusWalkerToBisulfiteCytosineWalker<Map<String,WiggleWriterReducible>, Map<String,WiggleWriterReducible>> {

//    @Output
//    @Multiplex(value=OutputMultiplexer.class,arguments={"outPrefix"})
//    Map<String,PrintStream> outStreams;


    @Argument(fullName = "outPrefix", shortName = "pre", doc = "Output prefix for all output files", required = true)
    public String outPrefix = null;
    
    @Argument(fullName = "outputAllContexts", shortName = "all", doc = "Output optional contexts HCH, GCG (default)", required = false)
    public boolean outputAllContexts = false;
    
    Map<String,WiggleWriterReducible> wigByContext = new HashMap<String,WiggleWriterReducible>();


	/**
	 * locus walker overrides
	 */


	/**
	 * Provides an initial value for the reduce function.  Hello walker counts loci,
	 * so the base case for the inductive step is 0, indicating that the walker has seen 0 loci.
	 * @return 0.
	 */
	@Override
	public Map<String,WiggleWriterReducible> reduceInit()
	{
		//Long out = 0L;
		Map<String,WiggleWriterReducible> out = new HashMap<String,WiggleWriterReducible>();
		return out;
	}


	@Override
	public void initialize() {
		// TODO Auto-generated method stub
		super.initialize();
		
		if (this.getToolkit().getArguments().numberOfThreads>1)
		{
			System.err.println("GnomeSeqToBareWigWalker does not yet implement multi-threaded mode. Use -nt 1");
			System.exit(1);
		}
		
		this.outputCph = true; // Because GNOME-seq used GCH
	}



	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.TreeReducible#treeReduce(java.lang.Object, java.lang.Object)
	 * 
	 * Must override this to get multithreaded capability
	 */
	@Override
	public Map<String,WiggleWriterReducible> treeReduce(Map<String,WiggleWriterReducible> a, Map<String,WiggleWriterReducible> b)
	{
		System.err.println("GnomeSeqToBareWigWalker does not yet implement multi-threaded mode. Use -nt 1");
		System.exit(1);
		
		return this.reduceCytosines(a, b);
	}

	/**
	 * Retrieves the final result of the traversal.
	 * @param result The ultimate value of the traversal, produced when map[n] is combined with reduce[n-1]
	 *               by the reduce function. 
	 */
	@Override
	public void onTraversalDone(Map<String,WiggleWriterReducible> result) 
	{
//		for (String context : result.keySet())
//		{
//			WiggleWriterReducible wig = result.get(context);
//			
//		}
	}

	
	/***************************************************
	 * cytosine walker overrides
	 ***************************************************/
	
	@Override
	protected void alertNewContig(String newContig) 
	{
	}

	@Override
	protected Map<String,WiggleWriterReducible> processCytosine(CpgBackedByGatk thisC)
	{
		writeCoverage(thisC);
		
		String context = thisC.context();

		// Skip optional contexts
		if (!this.outputAllContexts)
		{
			if (context.equalsIgnoreCase("HCH") || context.equalsIgnoreCase("GCG")) return this.wigByContext;
		}
		
		double meth = thisC.fracMeth(false);
		
		WiggleWriterReducible wig = null;
		if (this.wigByContext.containsKey(context))
		{
			wig = this.wigByContext.get(context);
			//logger.info("Found wig " + wig.toString());
		}
		else
		{
			// We use the wigByContext for merging
//			String outfn = String.format("%s.%s.wig.%s",this.outPrefix,context,this.wigByContext.hashCode());
			String name = String.format("%s.%s", this.outPrefix, context);
			String outfn = String.format("%s.wig",name,context);
			logger.info("NEW wig " + outfn);
			wig = new WiggleWriterReducible(new File(outfn));
			
			WiggleHeader head = new WiggleHeaderCytosines(name, name);
			wig.writeHeader(head);
			this.wigByContext.put(context, wig);
		}
		
		wig.writeData(thisC.getGenomeLoc(), Math.round(100.0*meth));
		//logger.info(String.format("\tWriting pos %d to wig \"%s\" (this=%s)\n", thisC.chromPos, context,this));
		
		return this.wigByContext;
	}


	private void writeCoverage(CpgBackedByGatk thisC)
	{
		String context = "readcvg";
		WiggleWriterReducible wig = null;
		if (this.wigByContext.containsKey(context))
		{
			wig = this.wigByContext.get(context);
			//logger.info("Found wig " + wig.toString());
		}
		else
		{
			// We use the wigByContext for merging
//			String outfn = String.format("%s.%s.wig.%s",this.outPrefix,context,this.wigByContext.hashCode());
			String name = String.format("%s.%s", this.outPrefix, context);
			String outfn = String.format("%s.wig",name,context);
			logger.info("NEW wig " + outfn);
			wig = new WiggleWriterReducible(new File(outfn));
			
			WiggleHeader head = new WiggleHeaderCytosines(name, name);
			wig.writeHeader(head);
			this.wigByContext.put(context, wig);
		}
		
		wig.writeData(thisC.getGenomeLoc(), thisC.totalReads);
	}


	@Override
	protected Map<String,WiggleWriterReducible> reduceCytosines(Map<String,WiggleWriterReducible> newMap, Map<String,WiggleWriterReducible> oldMap) 
	{
		Map<String,WiggleWriterReducible> outMap = null;
		
		if (newMap.hashCode() == oldMap.hashCode())
		{
			// Do nothing, same object
			outMap = oldMap;
		}
		else
		{
			logger.info(String.format("Actually combining maps: (%s,%s)\n",newMap.hashCode(),oldMap.hashCode()));
			outMap = this.actuallyReduceCytosines(newMap, oldMap);
		}
		
		return outMap;
	}
	
	private Map<String,WiggleWriterReducible> actuallyReduceCytosines(Map<String,WiggleWriterReducible> newMap, Map<String,WiggleWriterReducible> oldMap) 
	{
		Map<String,WiggleWriterReducible> outmap = oldMap; //new HashMap<String,WiggleWriterReducible>();
		
		Set<String> keys = new HashSet<String>();
		keys.addAll(newMap.keySet());
		keys.addAll(oldMap.keySet());
		for (String key : keys)
		{
			if (newMap.containsKey(key) && !oldMap.containsKey(key))
			{
				outmap.put(key, newMap.get(key));
			}
			else if (!newMap.containsKey(key) && oldMap.containsKey(key))
			{
				outmap.put(key, oldMap.get(key));
			}
			else if (!newMap.containsKey(key) && !oldMap.containsKey(key))
			{
				// SHOULD NOT BE HERE
				logger.error(String.format("GnomeSeqToBareWigWalker::actuallyReduceCytosines() - How did we get key \"%s\" when neither lhs and rhs contain it??",key));
				System.exit(1);
			}
			else
			{
				// They must both contain the key.  Merge
				//logger.info(String.format("About to combine a and b WigWriters for key \"%s\"",key)); 
				String outfn = String.format("%s.%s.wig.%s",this.outPrefix,key,this.hashCode());
				WiggleWriterReducible combined = WiggleWriterReducible.merge(newMap.get(key), oldMap.get(key), new File(outfn));
				outmap.put(key, combined);
			}		
		}
			
		return outmap;
	}
	
}


//class OutputMultiplexer implements Multiplexer<String>
//{
//
//	private String prefix = "prefix";
//	
//	public OutputMultiplexer(String inPrefix) {
//		super();
//		this.prefix = inPrefix;
//	}
//
//	@Override
//	public Collection<String> multiplex() {
//		return null;
//	}
//
//	@Override
//	public String transformArgument(String multiplexedEntry, String argument) {
//		return null;
//	}
//}



