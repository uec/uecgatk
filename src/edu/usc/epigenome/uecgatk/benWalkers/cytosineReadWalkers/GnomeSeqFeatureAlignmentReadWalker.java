package edu.usc.epigenome.uecgatk.benWalkers.cytosineReadWalkers;

import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.TreeMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

import org.apache.commons.math.fraction.Fraction;
import org.biojava.bio.seq.StrandedFeature;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.utils.collections.Pair;
import org.usckeck.genome.ChromFeatures;

import edu.usc.epigenome.genomeLibs.IupacPatterns;
import edu.usc.epigenome.genomeLibs.FeatAligners.FeatAlignerEachfeat;
import edu.usc.epigenome.genomeLibs.GenomicRange.GenomicRangeWithRefpoint;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.uecgatk.FractionNonidentical;
import edu.usc.epigenome.uecgatk.benWalkers.CpgBackedByGatkWithAlignmentRelCoords;
import edu.usc.epigenome.uecgatk.benWalkers.ReadWalkerToBisulfiteCytosineReadWalker;
import edu.usc.epigenome.uecgatk.benWalkers.ReadWithCpgMeths;
import edu.usc.epigenome.uecgatk.benWalkers.ReadWithCpgMethsQuadrants;


/**
 * @author benb
 * 
 * HOW DOES THE READ WALKER WORK WITH PAIRED-END READS?  WE NEED TO GET THEM IN ONE MAP STEP.
 * 
 */
public class GnomeSeqFeatureAlignmentReadWalker extends
		ReadWalkerToBisulfiteCytosineReadWalker<ReadWithCpgMeths, Map<String,FeatAlignerEachfeat>> {

    @Argument(fullName = "outPrefix", shortName = "pre", doc = "Output prefix for all output files", required = true)
    public String outPrefix = null;

    @Argument(fullName = "elementGff", shortName = "gff", doc = "A gff containing elements to interrogate", required = true)
    public String elementGff = null;

    @Argument(fullName = "iupacPatterns", shortName = "pats", doc = "A list of IUPAC contexts to interrogate", required = true)
    public String[] iupacPatterns = new String[]{"WCG","GCH"};

    @Argument(fullName = "downscaleCols", shortName = "cols", doc = "Number of columns in output alignment", required = false)
    public Integer downscaleCols = 0;

    @Argument(fullName = "windSize", shortName = "wind", doc = "window size", required = true)
    public Integer windSize = 0;

    @Argument(fullName = "threadedNfeatsStart", shortName = "as", doc = "Advanced. For multithreaded mode, the number of feats to initially start the array with (default 10)", required = false)
    public Integer threadedNfeatsStart = 10;
    
    @Argument(fullName = "threadedGrowsize", shortName = "ags", doc = "Advanced. For multithreaded mode, the number of feats to grow the arrays by when they overflow (default 100)", required = false)
    public Integer threadedGrowsize = 100;
    
    @Argument(fullName = "censor", shortName = "censor", doc = "Only count cytosines within the element itself", required = false)
    public boolean censor = false;

//    @Argument(fullName = "hcgMin", shortName = "hcgmin", doc = "Minimum methylation across HCGs in a read (0-1, default 0)", required = false)
//    public double hcgMin = 0.0;
//    @Argument(fullName = "hcgMax", shortName = "hcgmax", doc = "Maximum methylation across HCGs in a read (0-1, default 1)", required = false)
//    public double hcgMax = 1.0;
//    @Argument(fullName = "gchMin", shortName = "gchmin", doc = "Minimum methylation across GCHs in a read (0-1, default 0)", required = false)
//    public double gchMin = 0.0;
//    @Argument(fullName = "gchMax", shortName = "gchmax", doc = "Maximum methylation across GCHs in a read (0-1, default 1)", required = false)
//    public double gchMax = 1.0;
//

    @Argument(fullName = "contextCombos", shortName = "combos", doc = "Output combinations of contexts (default false)", required = false)
    public boolean contextCombos = false;
    
    
    @Argument(fullName = "minNumWcg", shortName = "minwcg", doc = "Minimum number of WCGs in the read to count (default 1)", required = false)
    public int minNumWcg = 1;
    @Argument(fullName = "minNumGch", shortName = "mingch", doc = "Minimum number of GCHs in the read to count (default 1)", required = false)
    public int minNumGch = 1;    
    
    
    // object variables

    // class variables
    protected static ChromFeatures feats = null;
    protected static boolean featsCompletelyPreloaded = false;
    protected static IupacPatterns patternMap = null;


    public GnomeSeqFeatureAlignmentReadWalker() {
		super();
	}


    @Override
    protected void alertNewContig(String newContig)
    {
    	System.err.printf("Number of threads = %d\n",this.getToolkit().getArguments().numberOfDataThreads);
    	if (this.getToolkit().getArguments().numberOfDataThreads == 1)
    	{
    		if (this.prevContig != null) this.feats.unload_coord_filtered_features();
    		if (this.feats != null)
    		{
    			int chrNum = feats.chrom_from_public_str(newContig);
    			this.feats.preload_coord_filtered_features(chrNum);
    		}
    	}
    	else
    	{
    		// Multi threads, load it up all at once so they won't have to fight over it.
    		synchronized(feats)
    		{
    			if (!featsCompletelyPreloaded)
    			{
    				System.err.println("Multi-threads, preloading coord filt featured (got lock)");
    				feats.preload_coord_filtered_features();
    				featsCompletelyPreloaded = true;
    			}
    		}
    	}
    }

    /**
	 * @return the nFeats
	 */
	public int getnFeats() {
		//System.err.printf("getting n feats: (%s) chr=%s\n", this.feats, this.prevContig);
		return feats.num_features(); // Total number of features across chroms
	}
    
	/***************************************************
	 * cytosine read walker overrides
	 ***************************************************/
	
//	public Map<String,FeatAlignerEachfeat> emptyMap()
//    {
//    	
//    	Map<String,FeatAlignerEachfeat> out = 
//			new HashMap<String,FeatAlignerEachfeat>();
//		for (String cond : this.iupacPatterns)
//		{
//			int nFeats = this.getnFeats();
//			//System.err.printf("Initializing aligner %s with %d elements\n", cond.context, nFeats);
//			
//			// **** CRITICAL SECTION.  In multithreaded mode, it spawns a new emptyMap() for almost every single locus,
//			// so creating a large datastructure is very inefficient.  FeatAlignerEachfeat is now selfgrowing,
//			// so we can create a small one to start out with
//			//if (nThreads<0) nThreads = this.getToolkit().getArguments().numberOfThreads;
//			if (this.getToolkit().getArguments().numberOfThreads>1)
//			{
//				nFeats = threadedNfeatsStart;
//			}
//			// *** END CRITICAL SECTION
//			
//			FeatAlignerEachfeat walker = new FeatAlignerEachfeat(windSize,this.censor, nFeats,this.downscaleCols);
//			walker.setArrayGrowsize(threadedGrowsize);
//			out.put(cond, walker);
//		}
//		
//		return out;
//	}

	public FeatAlignerEachfeat emptyAligner()
	{
		int nFeats = this.getnFeats();

		if (this.getToolkit().getArguments().numberOfDataThreads>1)
		{
			nFeats = threadedNfeatsStart;
		}

		FeatAlignerEachfeat walker = new FeatAlignerEachfeat(windSize,this.censor, nFeats,this.downscaleCols);
		walker.setArrayGrowsize(threadedGrowsize);

		return walker;
	}

	
	@Override
	public Map<String,FeatAlignerEachfeat>
	treeReduce(Map<String,FeatAlignerEachfeat> a, Map<String,FeatAlignerEachfeat> b)
	{
		
		if (a==null)
		{
			//System.err.println("______a==null_______");
			return b;
		}
		
		if (b==null)
		{
			//System.err.println("______b==null_______");
			return a;
		}
		
		// Go through each context in each one.
		Map<String,FeatAlignerEachfeat> out = new HashMap<String,FeatAlignerEachfeat>();
		Set<String> contexts = new HashSet(a.keySet());
		contexts.addAll(b.keySet());
		for (String context : contexts)
		{
			FeatAlignerEachfeat aAligner = a.get(context);
			FeatAlignerEachfeat bAligner = b.get(context);
			FeatAlignerEachfeat newAligner = null;
			if (aAligner==null)
			{
				newAligner = bAligner;
				//System.err.println("______aAligner==null_______");
			}
			else if (bAligner==null)
			{
				newAligner = aAligner;
				//System.err.println("______bAligner==null_______");
			}
			else
			{
				//System.err.printf("Merging FeatAligners with %d and %d elements\n",aAligner.numFeats(), bAligner.numFeats());
				newAligner = aAligner.mergeInAlignments(bAligner, false);//*** Not sure if there is a rule against modifying the pre-reduce objects.  It's much more efficient without.
				if (newAligner==null)
				{
					//System.err.println("______newAligner==null_______");
				}
			}
			
			if (newAligner != null)
			{
				out.put(context, newAligner);	
			}
		}
		
//		
//		// I am still getting errors from out-of-order cytosines given to the CpgWalker.
//		// **** FIX LATER ***
//		System.err.println("GnomeSeqAutocorrByReadWalker does not yet implement multi-threaded mode. Use -nt 1");
//		System.exit(1);
		
		//System.err.printf("__________Out=%s\n", out);
		return out;
	}

	
	/* (non-Javadoc)
	 * @see edu.usc.epigenome.uecgatk.benWalkers.ReadWalkerToBisulfiteCytosineReadWalker#processReadCytosines(java.util.List)
	 * 
	 * Called during Map phase.  The output is a map with keys context/chromPos pairs, with the average methylation within the 
	 * read as the map value
	 */
	@Override
	protected ReadWithCpgMeths processReadCytosines(ReadWithCpgMeths read) 
	{
		ReadWithCpgMeths newread = ReadWithCpgMeths.copy(read);
		if (read.size() == 0)
		{
//			System.err.printf("Found a read with 0 cytosines\n");	
		}
		else
		{
			newread.addAlignmentPoints(feats, this.windSize);
		}
		
		return newread;
		
	}

    /**
     * This can get called tons of times before even the first cytosine gets processed.
     * We can avoid the merges if we set it to null here and then make the merger know
     * how to deal with nulls.
     * @return
     */
    @Override
	public Map<String,FeatAlignerEachfeat> reduceInit()
	{
    	Map<String,FeatAlignerEachfeat> out = null; // this.emptyMap();
		
		return out;
	}
    
    
    @Override
    public Map<String,FeatAlignerEachfeat> reduce(ReadWithCpgMeths read, Map<String,FeatAlignerEachfeat> oldMap) {

    	if (oldMap == null)
    	{
    		//System.err.printf("oldMap==null, calling emptyMap()\n");
    		oldMap = new HashMap<String,FeatAlignerEachfeat>(); // this.emptyMap();
    	}
		int chromPos = read.midpointChromPos();

		Map<String,FeatAlignerEachfeat> outMap = new HashMap<String,FeatAlignerEachfeat>(oldMap);

		Map<String,FractionNonidentical> methVals = (this.contextCombos) ? ReadWithCpgMethsQuadrants.methLevelsFractions(read, patternMap) : 
				read.methLevelsFractions(patternMap);

		// Check if we have sufficient cytosines in different categories
		int nGch = 0;
		if (methVals.containsKey("GCH")) nGch = methVals.get("GCH").getDenominator();
		int nWcg = 0;
		if (methVals.containsKey("HCG")) nWcg = methVals.get("HCG").getDenominator();
		else if (methVals.containsKey("WCG")) nWcg = methVals.get("WCG").getDenominator();
		else if (methVals.containsKey("CCG")) nWcg = methVals.get("CCG").getDenominator();

		if ((nGch < this.minNumGch) || (nWcg < this.minNumWcg))
		{
			// System.err.printf("Found reads with <%d GCH or <%d WCG\n", this.minNumGch, this.minNumWcg);
		}
		else
		{
			for (String context : methVals.keySet())
			{
				double mLevel = methVals.get(context).doubleValue();
				FeatAlignerEachfeat aligner = oldMap.get(context);
				if (aligner == null) aligner = this.emptyAligner(); // Make it on the fly if we haven't seen it yet
				
				//    		System.err.printf("reduce(): pos=%d, context=%s, mLevel=%.2f\n",chromPos, context, mLevel);

				// We should only see each position in each feature once, so it should always
				// be empty when we encounter it here.
				// Meth val
				Iterator<GenomicRangeWithRefpoint> it = read.getAlignmentPoints();
				if (it != null)
				{
					while (it.hasNext())
					{
						GenomicRangeWithRefpoint ap = it.next();
						//    	    		System.err.printf("\treduce(ap=%d): pos=%d, context=%s, mLevel=%.2f\n", ap.getStart(), chromPos, context, mLevel);

						aligner.addAlignmentPos(
								chromPos,
								(read.getStrand() == StrandedFeature.NEGATIVE) ? Double.NaN : mLevel,
										(read.getStrand() == StrandedFeature.NEGATIVE) ? mLevel: Double.NaN,
												"", ap.getChrom(), ap.getRefPoint(), ap.getStrand(), 0);
					}
				}

				outMap.put(context, aligner);
			}
		}

		return outMap;
    }

	@Override
	public void initialize() {
		super.initialize();
		
		try
		{
			feats = new ChromFeatures(elementGff, true);
		}
		catch (Exception e)
		{
			System.err.printf("Could not read element file %s. Quitting. \n%s\n", elementGff, e.toString());
			e.printStackTrace();
			System.exit(1);
		}
		
		//logger.info(String.format("Initializing map (%d)",walkerByCondition.hashCode()));
		patternMap = new IupacPatterns();
		for (String pat : this.iupacPatterns)
		{
			patternMap.register(pat);
		}
		
		this.outputCph = true; // Because GNOME-seq used GCH
	}

	/**
	 * Retrieves the final result of the traversal.
	 * @param result The ultimate value of the traversal, produced when map[n] is combined with reduce[n-1]
	 *               by the reduce function. 
	 */
	@Override
	public void onTraversalDone(Map<String,FeatAlignerEachfeat> result) 
	{
		int count = 0;
		for (String cond : result.keySet())
		{
			if (result != null)
			{
				FeatAlignerEachfeat aligner = result.get(cond);

				try {

					String outfn = String.format("%s-%s.csv", this.outPrefix, cond.toString());
					PrintWriter pw = new PrintWriter(new FileOutputStream(outfn));

					aligner.matlabCsv(pw, false, true);
					pw.close();

				} catch (Exception e) {
					logger.error("Error FeatAlignerEachfeat::toCsvStr\n" + e.toString());
					e.printStackTrace();
				}


				count++;
			}
		}
	}

	public int mapTotal(Map<String,FeatAlignerEachfeat> in)
	{
		int grandTotal = 0;
		for (String cond : this.iupacPatterns)
		{
			FeatAlignerEachfeat feats = in.get(cond);
			grandTotal += feats.numFeats();
		}
		return grandTotal;
	}
	

	
	
}
