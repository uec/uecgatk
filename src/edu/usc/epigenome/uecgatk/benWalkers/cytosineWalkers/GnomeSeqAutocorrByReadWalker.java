package edu.usc.epigenome.uecgatk.benWalkers.cytosineWalkers;

import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.broadinstitute.sting.commandline.Argument;

import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerAllpairsAutocorrByread;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerAllpairsAutocorrByreadWcontext;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerParams;
import edu.usc.epigenome.uecgatk.benWalkers.CpgBackedByGatk;
import edu.usc.epigenome.uecgatk.benWalkers.LocusWalkerToBisulfiteCytosineWalker;


/**
 * @author benb
 * 
 * Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread>
 */
public class GnomeSeqAutocorrByReadWalker extends LocusWalkerToBisulfiteCytosineWalker<Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread>,Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread>>
{

    @Argument(fullName = "windSize", shortName = "ws", doc = "Output window size (default 2kb)", required = false)
    public int windSize = 2000;

    @Argument(fullName = "outPrefix", shortName = "pre", doc = "Output prefix for all output files", required = true)
    public String outPrefix = null;


    Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread> walkerByCondition = null;
    
    
    public Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread> emptyMap()
    {
		Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread> out = 
			new HashMap<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread>();
		for (AutocorrConditions cond : AutocorrConditions.values())
		{
			CpgWalkerAllpairsAutocorrByread walker = cond.createWalker(this.windSize);
			out.put(cond, walker);
		}
		
		return out;
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
		this.outputCph = true; // Because GNOME-seq used GCH
	}   

    /**
     * This can get called tons of times before even the first cytosine gets processed.
     * We can avoid the merges if we set it to null here and then make the merger know
     * how to deal with nulls.
     * @return
     */
    @Override
	public Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread> reduceInit()
	{
    	Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread> out = null; // this.emptyMap();
		
		return out;
	}
    
 


	@Override
	public Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread>
	treeReduce(Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread> a, Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread> b)
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
	public void onTraversalDone(Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread> result) 
	{
		int count = 0;
		for (AutocorrConditions cond : AutocorrConditions.values())
		{
			CpgWalkerAllpairsAutocorrByread walker = this.walkerByCondition.get(cond);
			
			try {

				String outfn = String.format("%s-%s.csv", this.outPrefix, cond.toString());
				PrintWriter pw = new PrintWriter(new FileOutputStream(outfn));
				pw.println(walker.headerStr());

				pw.println(walker.toCsvStr());
				pw.close();
			} catch (Exception e) {
				logger.error("Error CpgWalkerAllpairsAutocorrByread::toCsvStr\n" + e.toString());
				e.printStackTrace();
			}
			

			count++;
		}
	}

	
	/***************************************************
	 * cytosine walker overrides
	 ***************************************************/
	
	@Override
	protected void alertNewContig(String newContig) 
	{
		for (AutocorrConditions cond : AutocorrConditions.values())
		{
			CpgWalkerAllpairsAutocorrByread walker = this.walkerByCondition.get(cond);
			walker.setCurChr(newContig);
		}
	}

	
	@Override
	protected Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread> 
	processCytosine(CpgBackedByGatk thisC)
	{
		int grandTotal = 0;
		for (AutocorrConditions cond : AutocorrConditions.values())
		{
			CpgWalkerAllpairsAutocorrByread walker = this.walkerByCondition.get(cond);
			walker.streamCpg(thisC);
			grandTotal += walker.totalCount();
		}
		
//		logger.info(String.format("processCytosine(%s,%d) returning walker with %d pairs (this=%s, map=%s)",
//				thisC.getChrom(),thisC.chromPos, grandTotal, this.hashCode(), this.walkerByCondition.hashCode()));
		return this.walkerByCondition;
	}


	@Override
	protected Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread> 
	reduceCytosines(Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread> newMap, 
			Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread> oldMap) 
	{
		Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread> outMap = null;
		
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
	
	private Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread>
	actuallyReduceCytosines(Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread> newMap, 
			Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread> oldMap) 
	{
		Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread> outmap = newMap; //new HashMap<String,WiggleWriterReducible>();
		
		Set<GnomeSeqAutocorrByReadWalker.AutocorrConditions> keys = new HashSet<GnomeSeqAutocorrByReadWalker.AutocorrConditions>();
		keys.addAll(newMap.keySet());
		keys.addAll(oldMap.keySet());
		for (GnomeSeqAutocorrByReadWalker.AutocorrConditions key : keys)
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
				CpgWalkerAllpairsAutocorrByread combined = CpgWalkerAllpairsAutocorrByread.merge(newMap.get(key), oldMap.get(key));
				outmap.put(key, combined);
			}		
		}
			
		return outmap;
	}
	


	public static int mapTotal(Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread> in)
	{
		int grandTotal = 0;
		for (AutocorrConditions cond : AutocorrConditions.values())
		{
			CpgWalkerAllpairsAutocorrByread walker = in.get(cond);
			grandTotal += walker.totalCount();
		}
		return grandTotal;
	}

	public enum AutocorrConditions
	{
		HCG_HCG_SAMEREAD ("HCG", "HCG", true, false),
		HCG_GCH_SAMEREAD ("HCG", "GCH", true, false),
		GCH_HCG_SAMEREAD ("GCH", "HCG", true, false),
		GCH_GCH_SAMEREAD ("GCH", "GCH", true, false),
		HCG_HCG_ANYREAD ("HCG", "HCG", false, false),
		HCG_GCH_ANYREAD ("HCG", "GCH", false, false),
		GCH_HCG_ANYREAD ("GCH", "HCG", false, false),
		GCH_GCH_ANYREAD ("GCH", "GCH", false, false);
		
		private final String fromContext;
		private final String toContext;
		private final boolean sameRead;
		private final boolean sameStrand;

		private AutocorrConditions(String inFromContext, String inToContext, boolean inSameRead, boolean inSameStrand) 
		{
			this.sameRead = inSameRead;
			this.fromContext = inFromContext;
			this.toContext = inToContext;
			this.sameStrand = inSameStrand;
		}
		
		public CpgWalkerAllpairsAutocorrByread createWalker(int windSize)
		{
			CpgWalkerParams params = new CpgWalkerParams();
			params.maxScanningWindSize = windSize;
			params.minScanningWindSize = windSize;
			params.minScanningWindCpgs = 2;

			CpgWalkerAllpairsAutocorrByread out = new CpgWalkerAllpairsAutocorrByreadWcontext(
					params,
					this.sameStrand,
					false,
					this.sameRead,
					false,
					fromContext,
					toContext);
			out.useSummarizers(false); // These summarizers don't work with merging
			
			out.useOnlyCG(false);
			return out;
		}
		
	}
}