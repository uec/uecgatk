package edu.usc.epigenome.uecgatk.benWalkers.cytosineWalkers;

import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.biojava.bio.seq.StrandedFeature;
import org.broadinstitute.sting.commandline.Argument;
import org.usckeck.genome.ChromFeatures;

import edu.usc.epigenome.genomeLibs.FeatAligners.FeatAlignerEachfeat;
import edu.usc.epigenome.genomeLibs.GenomicRange.GenomicRangeWithRefpoint;
import edu.usc.epigenome.uecgatk.benWalkers.CpgBackedByGatk;
import edu.usc.epigenome.uecgatk.benWalkers.CpgBackedByGatkWithAlignmentRelCoords;
import edu.usc.epigenome.uecgatk.benWalkers.LocusWalkerToBisulfiteCytosineWalker;


/**
 * @author benb
 * 
 * Map<GnomeSeqFeatureAlignmentWalker.MethConditions,FeatAlignerEachfeat>
 */
public class GnomeSeqFeatureAlignmentWalker extends LocusWalkerToBisulfiteCytosineWalker<
CpgBackedByGatkWithAlignmentRelCoords,
Map<GnomeSeqFeatureAlignmentWalker.MethConditions,FeatAlignerEachfeat>>
{

    @Argument(fullName = "outPrefix", shortName = "pre", doc = "Output prefix for all output files", required = true)
    public String outPrefix = null;

    @Argument(fullName = "elementGff", shortName = "gff", doc = "A gff containing elements to interrogate", required = true)
    public String elementGff = null;

    @Argument(fullName = "downscaleCols", shortName = "cols", doc = "Number of columns in output alignment", required = false)
    public Integer downscaleCols = 0;

    @Argument(fullName = "windSize", shortName = "wind", doc = "window size", required = true)
    public Integer windSize = 0;

    @Argument(fullName = "censor", shortName = "censor", doc = "Only count cytosines within the element itself", required = false)
    public boolean censor = false;

//    Map<GnomeSeqFeatureAlignmentWalker.MethConditions,FeatAlignerEachfeat> alignerByCondition = null;
 //   protected int nFeats = 0;
    protected ChromFeatures feats = null;

    
    /**
	 * @return the nFeats
	 */
	public int getnFeats() {
		System.err.printf("getting n feats: (%s) chr=%s\n", this.feats, this.prevContig);
		return feats.num_features(this.feats.chrom_from_public_str(this.prevContig));
	}



	public Map<GnomeSeqFeatureAlignmentWalker.MethConditions,FeatAlignerEachfeat> emptyMap()
    {

    	
 //    	// Create arrays
//    	System.err.println("About to initialize aligners");
//    	fStatMats = new FeatAlignerEachfeat[nS][4];
//    	//		fDeltaMats = new FeatAligner[(nS*(nS-1))/2];
//    	//		fVarianceMat = new FeatAlignerAveraging(flankSize,false);
//    	int onDeltaMat = 0;
//    	for (int i = 0; i < nS; i++)
//    	{
//    		// 0=readCount, 1=nCpGs, 2=mLevel, 3= CpG weights
//    		if (readCounts) fStatMats[i][0] = new FeatAlignerEachfeat(flankSize,!this.censor, nFeats,this.downscaleCols); // Changed this to start with NaN.  Now we explicitly zero out features to allow censoring
//    		//			fStatMats[i][1] = new FeatAlignerEachfeat(flankSize,true, nFeats,500);
//    		if (!nometh) fStatMats[i][2] = new FeatAlignerEachfeat(flankSize,false, nFeats,this.downscaleCols);
//    		if (!nometh) fStatMats[i][3] = new FeatAlignerEachfeat(flankSize,false, nFeats,this.downscaleCols);
//    		//			for (int j = (i+1); j < nS; j++)
//    		//			{
//    		//				fDeltaMats[onDeltaMat++] = new FeatAlignerAveraging(flankSize,false);
//    		//			}
//    	}
    	
    	Map<GnomeSeqFeatureAlignmentWalker.MethConditions,FeatAlignerEachfeat> out = 
			new HashMap<GnomeSeqFeatureAlignmentWalker.MethConditions,FeatAlignerEachfeat>();
		for (MethConditions cond : MethConditions.values())
		{
			int nFeats = this.getnFeats();
			System.err.printf("Initializing aligner %s with %d elements\n", cond.context, nFeats);
			FeatAlignerEachfeat walker = cond.createAligner(this.windSize, nFeats, this.censor, this.downscaleCols);
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
		this.outputCph = true; // Because GNOME-seq used GCH
	}   

    /**
     * This can get called tons of times before even the first cytosine gets processed.
     * We can avoid the merges if we set it to null here and then make the merger know
     * how to deal with nulls.
     * @return
     */
    @Override
	public Map<GnomeSeqFeatureAlignmentWalker.MethConditions,FeatAlignerEachfeat> reduceInit()
	{
    	Map<GnomeSeqFeatureAlignmentWalker.MethConditions,FeatAlignerEachfeat> out = null; // this.emptyMap();
		
		return out;
	}
    
 


	@Override
	public Map<GnomeSeqFeatureAlignmentWalker.MethConditions,FeatAlignerEachfeat>
	treeReduce(Map<GnomeSeqFeatureAlignmentWalker.MethConditions,FeatAlignerEachfeat> a, Map<GnomeSeqFeatureAlignmentWalker.MethConditions,FeatAlignerEachfeat> b)
	{
		// I am still getting errors from out-of-order cytosines given to the CpgWalker.
		// **** FIX LATER ***
		System.err.println("GnomeSeqAutocorrByReadWalker does not yet implement multi-threaded mode. Use -nt 1");
		System.exit(1);
		
		return null;
	}

	/**
	 * Retrieves the final result of the traversal.
	 * @param result The ultimate value of the traversal, produced when map[n] is combined with reduce[n-1]
	 *               by the reduce function. 
	 */
	@Override
	public void onTraversalDone(Map<GnomeSeqFeatureAlignmentWalker.MethConditions,FeatAlignerEachfeat> result) 
	{
		int count = 0;
		for (MethConditions cond : MethConditions.values())
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

	
	/***************************************************
	 * cytosine walker overrides
	 ***************************************************/
	
	@Override
	protected void alertNewContig(String newContig) 
	{
		if (this.prevContig != null)
		{
			this.feats.unload_coord_filtered_features();
		}
		
		if (this.feats != null)
		{
			this.feats.preload_coord_filtered_features(feats.chrom_from_public_str(newContig));
		}
	}

	
	@Override
	protected CpgBackedByGatkWithAlignmentRelCoords 
	processCytosine(CpgBackedByGatk thisC)
	{
		int grandTotal = 0;
		
		CpgBackedByGatkWithAlignmentRelCoords out = 
			CpgBackedByGatkWithAlignmentRelCoords.copy(thisC);
		
		out.addAlignmentPoints(this.feats, this.windSize); // Including windSize flank will make sure we get all overlapping elements.
		
//		for (MethConditions cond : MethConditions.values())
//		{
//			FeatAlignerEachfeat walker = this.walkerByCondition.get(cond);
//		}
		
		return out;
	}


	@Override
	protected Map<GnomeSeqFeatureAlignmentWalker.MethConditions,FeatAlignerEachfeat> 
	reduceCytosines(CpgBackedByGatkWithAlignmentRelCoords cytosine, 
			Map<GnomeSeqFeatureAlignmentWalker.MethConditions,FeatAlignerEachfeat> oldMap) 
	{
		
		if (oldMap == null)
		{
			oldMap = this.emptyMap();
		}
		
		// First get the context
		GnomeSeqFeatureAlignmentWalker.MethConditions context = null;
		try 
		{
			context = GnomeSeqFeatureAlignmentWalker.MethConditions.valueOf(cytosine.context());
		}
		catch (Exception e)
		{
//			System.err.printf("Found unknown context: %s\n", cytosine.context());
		}

		Map<GnomeSeqFeatureAlignmentWalker.MethConditions,FeatAlignerEachfeat> outMap = oldMap;
		if (context != null)
		{
//			System.err.printf("Found KNOWN context: %s\n", cytosine.context());
			FeatAlignerEachfeat aligner = oldMap.get(context);

			// We should only see each position in each feature once, so it should always
			// be empty when we encounter it here.
			// Meth val
			double mLevel = cytosine.fracMeth(true);
			Iterator<GenomicRangeWithRefpoint> it = cytosine.getAlignmentPoints();
			while (it.hasNext())
			{
				GenomicRangeWithRefpoint ap = it.next();

				aligner.addAlignmentPos(
						cytosine.chromPos,
						(cytosine.getStrand() == StrandedFeature.NEGATIVE) ? Double.NaN : mLevel,
								(cytosine.getStrand() == StrandedFeature.NEGATIVE) ? mLevel: Double.NaN,
										"", ap.getChrom(), ap.getRefPoint(), ap.getStrand(), 0);
			}


			outMap = new HashMap<GnomeSeqFeatureAlignmentWalker.MethConditions,FeatAlignerEachfeat>(oldMap);
			outMap.put(context, aligner);
		}

		return outMap;
	}
	



	public static int mapTotal(Map<GnomeSeqFeatureAlignmentWalker.MethConditions,FeatAlignerEachfeat> in)
	{
		int grandTotal = 0;
		for (MethConditions cond : MethConditions.values())
		{
			FeatAlignerEachfeat feats = in.get(cond);
			grandTotal += feats.numFeats();
		}
		return grandTotal;
	}

	public enum MethConditions
	{
		HCG ("HCG"),
//		HCG_GCH_SAMEREAD ("HCG", "GCH", true, false),
		GCH ("GCH"),
//		GCH_GCH_SAMEREAD ("GCH", "GCH", true, false),
//		HCG_HCG_ANYREAD ("HCG", "HCG", false, false),
//		HCG_GCH_ANYREAD ("HCG", "GCH", false, false),
//		GCH_HCG_ANYREAD ("GCH", "HCG", false, false),
//		GCH_GCH_ANYREAD ("GCH", "GCH", false, false);
		;
	
		private final String context;
		private int count = 0;

		
		private MethConditions(String inContext) 
		{
			this.context = inContext;
		}
		
		public FeatAlignerEachfeat createAligner(int windSize, int nFeats, boolean censor, int downscaleCols)
		{


			FeatAlignerEachfeat out = new FeatAlignerEachfeat(windSize,censor, nFeats,downscaleCols);
			

			return out;
		}
		
//		static public MethConditions conditionByContext(String inContext)
//		{
//			MethConditions out = null;
//			for (MethConditions cond : MethConditions.values())
//			{
//				if (cond.context.equals(inContext))
//				{
//					out = cond;
//				}
//			}
//			return out;
//		}
		
	}
}