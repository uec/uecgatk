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

    @Argument(fullName = "elementGff", shortName = "gff", doc = "A gff containing elements to interrogate", required = true)
    public String elementGff = null;

    @Argument(fullName = "downscaleCols", shortName = "cols", doc = "Number of columns in output alignment", required = false)
    public Integer downscaleCols = 0;

    @Argument(fullName = "featWindSize", shortName = "wind", doc = "window size", required = true)
    public Integer featWindSize = 0;
    
    @Argument(fullName = "censor", shortName = "censor", doc = "Only count cytosines within the element itself", required = false)
    public boolean censor = false;

    @Argument(fullName = "useOnlyDifferentReads", shortName = "diff", doc = "only use those on different reads.  Otherwise, any is any read", required = false)
    public static boolean useOnlyDifferentReads = false;
    
    Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread> walkerByCondition = null;
    
    
    public Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread> emptyMap()
    {
		Map<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread> out = 
			new HashMap<GnomeSeqAutocorrByReadWalker.AutocorrConditions,CpgWalkerAllpairsAutocorrByread>();
		for (AutocorrConditions cond : AutocorrConditions.values())
		{
			CpgWalkerAllpairsAutocorrByread walker = cond.createWalker(this.windSize, this.elementGff, this.featWindSize, this.censor, this.downscaleCols);
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
				
				if (this.elementGff != null)
				{
					outfn = String.format("%s-featAligner-%s.csv", this.outPrefix, cond.toString());
					pw = new PrintWriter(new FileOutputStream(outfn));
					walker.getAligner().matlabCsv(pw, false, true);
					pw.close();
				}
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
//		HCG_HCG_SAMEREAD ("HCG", "HCG", true, false, true, true),
//		HCG_GCH_SAMEREAD ("HCG", "GCH", true, false, true, true),
//		GCH_HCG_SAMEREAD ("GCH", "HCG", true, false, true, true),
//		GCH_GCH_SAMEREAD ("GCH", "GCH", true, false, true, true),
//		HCG_HCG_ANYREAD ("HCG", "HCG", false, false, true, true),
//		HCG_GCH_ANYREAD ("HCG", "GCH", false, false, true, true),
//		GCH_HCG_ANYREAD ("GCH", "HCG", false, false, true, true),
//		GCH_GCH_ANYREAD ("GCH", "GCH", false, false, true, true),
		
		HCG_GCH_ANYREAD_00 ("HCG", "GCH", false, false, false, false),
		HCG_GCH_ANYREAD_01 ("HCG", "GCH", false, false, false, true),
		HCG_GCH_ANYREAD_10 ("HCG", "GCH", false, false, true, false),
		HCG_GCH_ANYREAD_11 ("HCG", "GCH", false, false, true, true),
		
		HCG_GCH_SAMEREAD_00 ("HCG", "GCH", true, false, false, false),
		HCG_GCH_SAMEREAD_01 ("HCG", "GCH", true, false, false, true),
		HCG_GCH_SAMEREAD_10 ("HCG", "GCH", true, false, true, false),
		HCG_GCH_SAMEREAD_11 ("HCG", "GCH", true, false, true, true),
		
		HCG_HCG_ANYREAD_00 ("HCG", "HCG", false, false, false, false),
		HCG_HCG_ANYREAD_01 ("HCG", "HCG", false, false, false, true),
		HCG_HCG_ANYREAD_10 ("HCG", "HCG", false, false, true, false),
		HCG_HCG_ANYREAD_11 ("HCG", "HCG", false, false, true, true),
		
		HCG_HCG_SAMEREAD_00 ("HCG", "HCG", true, false, false, false),
		HCG_HCG_SAMEREAD_01 ("HCG", "HCG", true, false, false, true),
		HCG_HCG_SAMEREAD_10 ("HCG", "HCG", true, false, true, false),
		HCG_HCG_SAMEREAD_11 ("HCG", "HCG", true, false, true, true);
		
		private final String fromContext;
		private final String toContext;
		private final boolean sameRead;
		private final boolean sameStrand;
		private final boolean firstMeth;
		private final boolean secondMeth;

		private AutocorrConditions(String inFromContext, String inToContext, boolean inSameRead, boolean inSameStrand, boolean inFirstMeth, boolean inSecondMeth) 
		{
			this.sameRead = inSameRead;
			this.fromContext = inFromContext;
			this.toContext = inToContext;
			this.sameStrand = inSameStrand;
			this.firstMeth = inFirstMeth;
			this.secondMeth = inSecondMeth;
		}
		
		public CpgWalkerAllpairsAutocorrByread createWalker(int windSize)
		{
			return createWalker(windSize, null, 0, false, 0);
		}
			
			public CpgWalkerAllpairsAutocorrByread createWalker(int windSize, String inGff, int inFeatWindSize, boolean inCensor, int inDownscaleCols)
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
					(useOnlyDifferentReads) ? !this.sameRead : false, // in different read
					fromContext,
					toContext);
			out.useSummarizers(false); // These summarizers don't work with merging
			out.useOnlyCG(false);
			
			if (inGff != null)
			{
				// Add gff
				out.enableFeatAlignment(inGff, inFeatWindSize, inCensor, inDownscaleCols, firstMeth, secondMeth);
			}
				
			
			return out;
		}
		
	}
}