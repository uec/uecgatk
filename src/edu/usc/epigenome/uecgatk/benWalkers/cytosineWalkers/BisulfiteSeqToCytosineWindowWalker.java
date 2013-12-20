package edu.usc.epigenome.uecgatk.benWalkers.cytosineWalkers;

import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.lang.reflect.Constructor;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.broadinstitute.sting.commandline.Argument;
import org.kohsuke.args4j.Option;

import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerParams;
import edu.usc.epigenome.uecgatk.benWalkers.CpgBackedByGatk;
import edu.usc.epigenome.uecgatk.benWalkers.CpgWalker;
import edu.usc.epigenome.uecgatk.benWalkers.LocusWalkerToBisulfiteCytosineWalker;


/**
 * @author benb
 * 
 * Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalker>
 */
public abstract class BisulfiteSeqToCytosineWindowWalker extends LocusWalkerToBisulfiteCytosineWalker<Map<BisulfiteSeqToCytosineWindowWalker.ContextConditions,CpgWalker>,Map<BisulfiteSeqToCytosineWindowWalker.ContextConditions,CpgWalker>>
{

    @Argument(fullName = "outPrefix", shortName = "pre", doc = "Output prefix for all output files", required = true)
    public String outPrefix = null;
     
    @Argument(fullName="minCpgs",shortName = "mincs", doc="minimum Cpgs (each strand counted separately) in window (8)", required = true)
    protected int minCpgs = 8;
    
    @Argument(fullName="maxWindStretch",shortName = "maxwind", doc="maximum amount to stretch window to find minCpgs cpgs (10000)", required = false)
    protected int maxWindStretch = 10000;

    @Argument(fullName = "gnomeMode", shortName = "gnome", doc = "GNOMe mode - outputs only GCH datapoints (default false)", required = false)
    public boolean gnomeMode = false;
    
    @Argument(fullName = "useFixedWind", shortName = "fixed", doc = "Use a fixed window (must set fixedWindMinLength)", required = false)
    public boolean useFixedWind = false;

    @Argument(fullName="fixedWindMinLength",shortName = "fixedlen", doc="Minimum fixed window length (5000)", required = false)
    protected int fixedWindMinLength = 5000;
   
    //  @Argument(fullName = "windSize", shortName = "w", doc = "minimum window size (500)", required = true)
//  public int windSize = 500;

    ////////////////////////////////
    // Private vars
    ////////////////////////////////
    protected Map<BisulfiteSeqToCytosineWindowWalker.ContextConditions,CpgWalker> walkerByCondition = null;
    protected String cpgWalkerType = null; // Must set this to something like CpgWalkerWindowWigWriter
    
    
    
    /**
	 * @param cpgWalkerType
	 */
	public BisulfiteSeqToCytosineWindowWalker(String cpgWalkerType) {
		super();
		this.cpgWalkerType = cpgWalkerType;
	}


	public Map<BisulfiteSeqToCytosineWindowWalker.ContextConditions,CpgWalker> emptyMap()
    {
		Map<BisulfiteSeqToCytosineWindowWalker.ContextConditions,CpgWalker> out = 
			new HashMap<BisulfiteSeqToCytosineWindowWalker.ContextConditions,CpgWalker>();
		for (ContextConditions cond : this.validContexts())
		{
			CpgWalker walker = cond.createWalker(this.minCpgs, this.maxWindStretch, this.cpgWalkerType, this.useFixedWind, this.fixedWindMinLength);
			out.put(cond, walker);
		}
		
		return out;
	}

    
	public ContextConditions[] validContexts()
	{
		return ContextConditions.validValues(this.gnomeMode, this.outputCph);
	}
	
	/**
	 * locus walker overrides
	 */

    @Override
	public void initialize() {
		// TODO Auto-generated method stub
		super.initialize();
		walkerByCondition = this.emptyMap();
		//logger.info(String.format("Initializing map (%d)",walkerByCondition.hashCode()));
		this.outputCph = this.gnomeMode; // Because GNOME-seq used GCH
	}   

    /**
     * This can get called tons of times before even the first cytosine gets processed.
     * We can avoid the merges if we set it to null here and then make the merger know
     * how to deal with nulls.
     * @return
     */
    @Override
	public Map<BisulfiteSeqToCytosineWindowWalker.ContextConditions,CpgWalker> reduceInit()
	{
    	Map<BisulfiteSeqToCytosineWindowWalker.ContextConditions,CpgWalker> out = null; // this.emptyMap();
		
		return out;
	}
    
 


	@Override
	public Map<BisulfiteSeqToCytosineWindowWalker.ContextConditions,CpgWalker>
	treeReduce(Map<BisulfiteSeqToCytosineWindowWalker.ContextConditions,CpgWalker> a, Map<BisulfiteSeqToCytosineWindowWalker.ContextConditions,CpgWalker> b)
	{
		// I am still getting errors from out-of-order cytosines given to the CpgWalker.
		// **** FIX LATER ***
		System.err.println("GnomeSeqAutocorrByReadWalker does not yet implement multi-threaded mode. Use -nt 1");
		System.exit(1);
		
		return this.reduceCytosines(a, b);
	}

	/**
	 * Retrieves the final result of the traversal.
	 * @param result The ultimate value of the traversal, produced when map[n] is combined with reduce[n-1]
	 *               by the reduce function. 
	 */
	@Override
	public void onTraversalDone(Map<BisulfiteSeqToCytosineWindowWalker.ContextConditions,CpgWalker> result) 
	{
		int count = 0;
		for (ContextConditions cond : this.validContexts())
		{
			CpgWalker walker = this.walkerByCondition.get(cond);
			walker.finishChr();
			count++;
		}
	}

	
	/***************************************************
	 * cytosine walker overrides
	 ***************************************************/
	
	@Override
	protected void alertNewContig(String newContig) 
	{
		for (ContextConditions cond : this.validContexts())
		{
			CpgWalker walker = this.walkerByCondition.get(cond);
			if (this.prevContig != null) walker.finishChr();
			walker.setCurChr(newContig);
		}
	}

	
	@Override
	protected Map<BisulfiteSeqToCytosineWindowWalker.ContextConditions,CpgWalker> 
	processCytosine(CpgBackedByGatk thisC)
	{
		int grandTotal = 0;
		for (ContextConditions cond : this.validContexts())
		{
			CpgWalker walker = this.walkerByCondition.get(cond);
			walker.streamCpg(thisC);
		}
		
//		logger.info(String.format("processCytosine(%s,%d) returning walker with %d pairs (this=%s, map=%s)",
//				thisC.getChrom(),thisC.chromPos, grandTotal, this.hashCode(), this.walkerByCondition.hashCode()));
		return this.walkerByCondition;
	}


	@Override
	protected Map<BisulfiteSeqToCytosineWindowWalker.ContextConditions,CpgWalker> 
	reduceCytosines(Map<BisulfiteSeqToCytosineWindowWalker.ContextConditions,CpgWalker> newMap, 
			Map<BisulfiteSeqToCytosineWindowWalker.ContextConditions,CpgWalker> oldMap) 
	{
		Map<BisulfiteSeqToCytosineWindowWalker.ContextConditions,CpgWalker> outMap = null;
		
		if (newMap == null)
		{
			outMap = oldMap;
		}
		else if (oldMap == null)
		{
			outMap = newMap;
		}
		else if (newMap.hashCode() == oldMap.hashCode())
		{
			// Do nothing, same object
			outMap = oldMap;
		}
		else
		{
			logger.info(String.format("Actually combining maps: (%s,%d pairs) + (%s, %d pairs)",
					newMap.hashCode(),mapTotal(newMap),oldMap.hashCode(),mapTotal(oldMap)));
			
			try
			{
				throw new Exception("FAIL");
			}
			catch (Exception e)
			{
//				e.printStackTrace();
//				System.exit(1);
			}
			
			outMap = this.actuallyReduceCytosines(newMap, oldMap);
		}
		
		return outMap;
	}
	
	private Map<BisulfiteSeqToCytosineWindowWalker.ContextConditions,CpgWalker>
	actuallyReduceCytosines(Map<BisulfiteSeqToCytosineWindowWalker.ContextConditions,CpgWalker> newMap, 
			Map<BisulfiteSeqToCytosineWindowWalker.ContextConditions,CpgWalker> oldMap) 
	{
		Map<BisulfiteSeqToCytosineWindowWalker.ContextConditions,CpgWalker> outmap = newMap; //new HashMap<String,WiggleWriterReducible>();
		
//		Set<BisulfiteSeqToCytosineWindowWalker.ContextConditions> keys = new HashSet<BisulfiteSeqToCytosineWindowWalker.ContextConditions>();
//		keys.addAll(newMap.keySet());
//		keys.addAll(oldMap.keySet());
//		for (BisulfiteSeqToCytosineWindowWalker.ContextConditions key : keys)
//		{
//			if (newMap.containsKey(key) && !oldMap.containsKey(key))
//			{
//				outmap.put(key, newMap.get(key));
//			}
//			else if (!newMap.containsKey(key) && oldMap.containsKey(key))
//			{
//				outmap.put(key, oldMap.get(key));
//			}
//			else if (!newMap.containsKey(key) && !oldMap.containsKey(key))
//			{
//				// SHOULD NOT BE HERE
//				logger.error(String.format("GnomeSeqToBareWigWalker::actuallyReduceCytosines() - How did we get key \"%s\" when neither lhs and rhs contain it??",key));
//				System.exit(1);
//			}
//			else
//			{
//				// They must both contain the key.  Merge
//				//logger.info(String.format("About to combine a and b WigWriters for key \"%s\"",key)); 
//				String outfn = String.format("%s.%s.wig.%s",this.outPrefix,key,this.hashCode());
//				CpgWalker combined = CpgWalker.merge(newMap.get(key), oldMap.get(key));
//				outmap.put(key, combined);
//			}		
//		}
			
		return outmap;
	}
	


	public  int mapTotal(Map<BisulfiteSeqToCytosineWindowWalker.ContextConditions,CpgWalker> in)
	{
		int grandTotal = 0;
		for (ContextConditions cond : this.validContexts())
		{
			CpgWalker walker = in.get(cond);
		}
		return grandTotal;
	}

	public enum ContextConditions
	{
		HCG ("HCG"),
		GCH ("GCH"),
		GCG ("GCG"),
		HCH ("HCH"),
		CG ("CG"),
		CH ("CH");
		
		private String context;

		private ContextConditions(String inContext) 
		{
			this.context = inContext;
		}
		
		public static ContextConditions[] validValues(boolean gnomeMode, boolean useCph)
		{
			ContextConditions[] out = null;
			
			if (gnomeMode)
			{
				out = new ContextConditions[]{HCG,GCH};
			}
			else if (useCph)
			{
				out = new ContextConditions[]{CG,CH};
			}
			else
			{
				out = new ContextConditions[]{CG};
			}
			return out;
		}
		
		public CpgWalker createWalker(int minCs, int maxWindStretch, String cpgWalkerClassName)
		{
			// Default to variable window
			return createWalker(minCs, maxWindStretch, cpgWalkerClassName, false, 0);
		}
		
		public CpgWalker createWalker(int minCs, int maxWindStretch, String cpgWalkerClassName, boolean useFixedWind, int fixedWindMinSize)
		{

			CpgWalkerParams params = new CpgWalkerParams();
			params.useVariableWindow = !useFixedWind;
			params.maxScanningWindSize = maxWindStretch;
			params.minScanningWindCpgs = minCs;
			if (useFixedWind)
			{
				params.minScanningWindSize = fixedWindMinSize;
			}
			
			CpgWalker out = null;
			
			try
			{
				Class<?> cpgWalkerClass = Class.forName(cpgWalkerClassName);
				Constructor constructor = cpgWalkerClass.getConstructor(CpgWalkerParams.class);
				out = (CpgWalker)constructor.newInstance(params);
			}
			catch (Exception e)
			{
				System.err.printf("Fatal error, could not make class of type %s\n%s\n",cpgWalkerClassName,e.toString());
				e.printStackTrace();
				System.exit(1);
			}
			//out.walkParams = params;

			return out;
		}
		
	}
}