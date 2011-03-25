package edu.usc.epigenome.uecgatk.benWalkers.cytosineWalkers;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.broadinstitute.sting.commandline.Argument;
import org.kohsuke.args4j.Option;

import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerAllpairsAutocorrByread;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerParams;
import edu.usc.epigenome.uecgatk.WiggleWriterReducible;
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


    Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread> walkerByCondition = 
    	new HashMap<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread>();
    
    
	/**
	 * locus walker overrides
	 */


    @Override
	public Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread> reduceInit()
	{
		//Long out = 0L;
		Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread> out = 
			new HashMap<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread>();
		
		for (AutocorrConditions cond : AutocorrConditions.values())
		{
			CpgWalkerAllpairsAutocorrByread walker = cond.createWalker(this.windSize);
			out.put(cond, walker);
		}
		
		return out;
	}



	@Override
	public Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread>
	treeReduce(Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread> a, Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread> b)
	{
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
	}

	
	/***************************************************
	 * cytosine walker overrides
	 ***************************************************/
	
	@Override
	protected Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread> processCytosine(CpgBackedByGatk thisC)
	{
		String context = thisC.context();
		double meth = thisC.fracMeth(false);

		
		return this.walkerByCondition;
	}


	@Override
	protected Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread> 
	reduceCytosines(Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread> newMap, 
			Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread> oldMap) 
	{
		Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread> outMap = null;
		
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
	
	private Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread>
	actuallyReduceCytosines(Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread> newMap, 
			Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread> oldMap) 
	{
		Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread> outmap = oldMap; //new HashMap<String,WiggleWriterReducible>();
		
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
	




	public enum AutocorrConditions
	{
		HCG_HCG_SAMEREAD ("HCG", "HCG", true, false),
		HCG_GCH_SAMEREAD ("HCG", "GCH", true, false),
		GCH_GCH_SAMEREAD ("GCH", "GCH", true, false),
		HCG_HCG_ANYREAD ("HCG", "HCG", false, false),
		HCG_GCH_ANYREAD ("HCG", "GCH", false, false),
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

			CpgWalkerAllpairsAutocorrByread out = new CpgWalkerAllpairsAutocorrByread(
					params,
					this.sameStrand,
					false,
					this.sameRead,
					false);
			return out;
		}
		
	}
}