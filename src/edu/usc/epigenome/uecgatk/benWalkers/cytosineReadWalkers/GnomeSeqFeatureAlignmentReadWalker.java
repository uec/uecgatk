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

import edu.usc.epigenome.genomeLibs.FeatAligners.FeatAlignerEachfeat;
import edu.usc.epigenome.genomeLibs.GenomicRange.GenomicRangeWithRefpoint;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.uecgatk.FractionNonidentical;
import edu.usc.epigenome.uecgatk.IupacPatterns;
import edu.usc.epigenome.uecgatk.benWalkers.CpgBackedByGatkWithAlignmentRelCoords;
import edu.usc.epigenome.uecgatk.benWalkers.ReadWalkerToBisulfiteCytosineReadWalker;
import edu.usc.epigenome.uecgatk.benWalkers.ReadWithCpgMeths;


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
    	System.err.printf("Number of threads = %d\n",this.getToolkit().getArguments().numberOfThreads);
    	if (this.getToolkit().getArguments().numberOfThreads == 1)
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
	
	public Map<String,FeatAlignerEachfeat> emptyMap()
    {
    	
    	Map<String,FeatAlignerEachfeat> out = 
			new HashMap<String,FeatAlignerEachfeat>();
		for (String cond : this.iupacPatterns)
		{
			int nFeats = this.getnFeats();
			//System.err.printf("Initializing aligner %s with %d elements\n", cond.context, nFeats);
			
			// **** CRITICAL SECTION.  In multithreaded mode, it spawns a new emptyMap() for almost every single locus,
			// so creating a large datastructure is very inefficient.  FeatAlignerEachfeat is now selfgrowing,
			// so we can create a small one to start out with
			//if (nThreads<0) nThreads = this.getToolkit().getArguments().numberOfThreads;
			if (this.getToolkit().getArguments().numberOfThreads>1)
			{
				nFeats = threadedNfeatsStart;
			}
			// *** END CRITICAL SECTION
			
			FeatAlignerEachfeat walker = new FeatAlignerEachfeat(windSize,this.censor, nFeats,this.downscaleCols);
			walker.setArrayGrowsize(threadedGrowsize);
			out.put(cond, walker);
		}
		
		return out;
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
		
//		FracPair outPair = new FracPair(false);
//		for (Cpg c : cs)
//		{
//			int positionInRead = c.chromPos;
//			String cContext = c.context();
//			boolean isMethylated = (c.cReads>0);
//			
//			if (cContext.equalsIgnoreCase("GCH"))
//			{
//				outPair.incrementGCH(isMethylated);
//			}
//			else if (cContext.equalsIgnoreCase("HCG"))
//			{
//				outPair.incrementHCG(isMethylated);
//			}
//			else
//			{
//				// Ignore other context
//				//System.err.printf("Ignoring cytosine: %d, %s, %s\n", positionInRead, cContext, isMethylated);
//			}
//		}

//		// Check if we're in range and output. HCG is first
//		if ( (outPair.hcg().getDenominator() >= this.minNumHcg) && (outPair.gch().getDenominator() >= this.minNumGch) )
//		{
//			if ( (outPair.hcg().doubleValue() >= this.hcgMin) && (outPair.hcg().doubleValue() <= this.hcgMax ) &&
//					(outPair.gch().doubleValue() >= this.gchMin) && (outPair.gch().doubleValue() <= this.gchMax ) )
//			{
//				// Bonanza.  Just use the middle of the cytosines as the point.  They should be ordered.
//				int midIndex = (int)Math.round((double)cs.size()/2.0);
//				int point = cs.get(midIndex).chromPos;
//				String chr = this.prevContig;
//				
//				// Use the methyldb schema
//				out.printf("%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
//						chr, point, "+", 1, 1, 0, 0, 0, 0, 0, 1, 1, 0);
//			}
//		}
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
    		oldMap = this.emptyMap();
    	}
		int chromPos = read.midpointChromPos();
    	Map<String,Double> methVals = read.methLevels(patternMap);
    	
    	Map<String,FeatAlignerEachfeat> outMap = new HashMap<String,FeatAlignerEachfeat>(oldMap);
    	for (String context : methVals.keySet())
    	{
    		double mLevel = methVals.get(context);
    		FeatAlignerEachfeat aligner = oldMap.get(context);
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
		for (String cond : this.iupacPatterns)
		{
			if (result != null)
			{
				FeatAlignerEachfeat aligner = result.get(cond);

				try {

					String outfn = String.format("%s-%s.csv", this.outPrefix, cond.toString());
					PrintWriter pw = new PrintWriter(new FileOutputStream(outfn));

					aligner.matlabCsv(pw, false);
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
	
	// -- // -- // -- // -- // -- 
	
	
//	public class FracPair extends Pair<FractionNonidentical,FractionNonidentical> implements Comparable<FracPair>
//	{
//		public FracPair(boolean mergeIdentical) {
//			super(new FractionNonidentical(), new FractionNonidentical());
//			this.first.setUseIdentical(mergeIdentical);
//			this.second.setUseIdentical(mergeIdentical);
//		}
//
//		public FracPair(FractionNonidentical hcg, FractionNonidentical gch, boolean mergeIdentical) {
//			super(hcg, gch);
//			this.first.setUseIdentical(mergeIdentical);
//			this.second.setUseIdentical(mergeIdentical);
//		}
//		
//		public void incrementHCG(boolean methylated)
//		{
//			increment(false, methylated);
//		}
//		
//		public void incrementGCH(boolean methylated)
//		{
//			increment(true, methylated);
//		}
//		
//		public FractionNonidentical hcg()
//		{
//			return this.first;
//		}
//		
//		public FractionNonidentical gch()
//		{
//			return this.second;
//		}
//
//		private void increment(boolean gch, boolean methylated)
//		{
//			FractionNonidentical frac = (gch) ? this.second : this.first;
//			frac.incDenominator();
//			if (methylated) frac.incNumerator();
//		}
//
//		@Override
//		public String toString() 
//		{
//			String str = String.format("HCG=%s, GCH=%s", this.first.toString(), this.second.toString());
//			return str;
//		}
//
//		@Override
//		public int compareTo(FracPair arg0) {
//			int result = this.first.compareTo(arg0.first);
//			
//			if (result == 0)
//			{
//				result = this.second.compareTo(arg0.second);
//			}
//			//System.err.printf("CompareTo(%s,%s)=%d\n",this,arg0,result);
//			return result;
//		}
//		
//		
//		
//	}
	
	
}
