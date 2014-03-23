package edu.usc.epigenome.uecgatk.benWalkers;

//import java.util.ArrayList;
//import java.util.Arrays;
//import java.util.HashMap;
//import java.util.List;
//import java.util.Map;
//import java.util.TreeSet;

//import com.mallardsoft.tuple.Pair;
//import com.mallardsoft.tuple.Quadruple;
//import com.mallardsoft.tuple.Tuple;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;

import org.biojava.bio.seq.StrandedFeature;
import org.usckeck.genome.ChromFeatures;

import edu.usc.epigenome.genomeLibs.MatUtils;
import edu.usc.epigenome.genomeLibs.FeatAligners.AlignmentRelCoords;
import edu.usc.epigenome.genomeLibs.FeatAligners.FeatAlignerEachfeat;
import edu.usc.epigenome.genomeLibs.GenomicRange.GenomicRangeWithRefpoint;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerParams;
import edu.usc.epigenome.uecgatk.benWalkers.cytosineWalkers.Cpg;
import edu.usc.epigenome.uecgatk.benWalkers.cytosineWalkers.CpgRead;


public class CpgWalkerAllpairsAutocorrByread extends CpgWalkerAllpairs {

//	protected double[] binEdges = null;
//	protected TreeSet<Double> binTree = null;
//	protected Map<Object, int[]> counters = null;
	
	protected static final int SAME_READ = 1;
	protected static final int DIFFERENT_READ = 2;
	protected static final int ANY_READ = 3;
	
	protected int readType = ANY_READ;
	
	protected boolean useOnlyCG = true;
	protected final static boolean useOnlyOneDifferentRead = true; // I found that having a wild imbalance between the number of reads considered in SAME_READ and DIFFERENT_READ could skew results.
	
	
	// By dist version
	protected int[] nMM;
	protected int[] nMU;
	protected int[] nUM;
	protected int[] nUU; 
	protected int[] nM;
	protected int[] nU;
	
	// Feat alignment version
    protected FeatAlignerEachfeat aligner = null;
    protected ChromFeatures feats = null;
    protected boolean featsCompletelyPreloaded = false;
    protected String prevChr = null;
    protected boolean firstMeth = false;
    protected boolean secondMeth = false;
   
	
	public CpgWalkerAllpairsAutocorrByread(CpgWalkerParams inWalkParams, boolean inSamestrandOnly, boolean inOppstrandOnly,int inReadType) {
		super(inWalkParams, inSamestrandOnly, inOppstrandOnly);
		
		this.readType = inReadType;

		init();
	}
	
	public CpgWalkerAllpairsAutocorrByread(CpgWalkerParams inWalkParams, boolean inSamestrandOnly, boolean inOppstrandOnly, boolean inSameRead, boolean inDifferentRead) {
		super(inWalkParams, inSamestrandOnly, inOppstrandOnly);
		
		if (inSameRead && inDifferentRead)
		{
			System.err.println("CpgWalkerAllpairsAutocorrByread::constructor can not be called with both inSameRead and inDifferentRead set");
			(new Exception()).printStackTrace();
			System.exit(1);
		}
		else if (inSameRead)
		{
			readType = SAME_READ;
		}
		else if (inDifferentRead)
		{
			readType = DIFFERENT_READ;
		}
		else
		{
			readType = ANY_READ;
		}
		
		init();
	}
	
	public void enableFeatAlignment(String gfffn, int featWindSize, boolean censor, int downscaleCols, boolean inFirstMeth, boolean inSecondMeth)
	{
		try
		{
			feats = new ChromFeatures(gfffn, true);
			int nFeats = feats.num_features();
			this.firstMeth = inFirstMeth;
			this.secondMeth = inSecondMeth;
			
			// Initialize the aligner. Preset all feats to zero
			this.aligner = new FeatAlignerEachfeat(featWindSize,censor, nFeats,downscaleCols);
			System.err.println("Built aligner with " + aligner.numFeats() + " feats");
		}
		catch (Exception e)
		{
			System.err.printf("Could not read element file %s. Quitting. \n%s\n", gfffn, e.toString());
			e.printStackTrace();
			System.exit(1);
		}		
	}

	
	
    /**
	 * @return the aligner
	 */
	public FeatAlignerEachfeat getAligner() {
		return aligner;
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalker#alertNewChrom()
	 */
	@Override
	protected void alertNewChrom() {
		// TODO Auto-generated method stub
		super.alertNewChrom();
		
// ONLY WORKS SINGLE THREADED
    		if (prevChr != null) feats.unload_coord_filtered_features();
    		if (feats != null)
    		{
    			int chrNum = feats.chrom_from_public_str(this.curChr);
    			feats.preload_coord_filtered_features(chrNum);
    		}
    		prevChr = this.curChr;
    		
//    	else
//    	{
//    		// Multi threads, load it up all at once so they won't have to fight over it.
//    		synchronized(feats)
//    		{
//    			if (!featsCompletelyPreloaded)
//    			{
//    				System.err.println("Multi-threads, preloading coord filt featured (got lock)");
//    				feats.preload_coord_filtered_features();
//    				featsCompletelyPreloaded = true;
//    			}
//    		}
//    	}	
    		
	}

    
	/**
	 * This makes it reducible for map-reduce
	 * 
	 * @param cpgWalkerAllpairsAutocorrByread
	 * @param cpgWalkerAllpairsAutocorrByread2
	 * @return
	 */
	public static CpgWalkerAllpairsAutocorrByread merge(
			CpgWalkerAllpairsAutocorrByread a,
			CpgWalkerAllpairsAutocorrByread b) {
		
				System.err.printf("Merging Autcorr(%d) + Autocorr(%d)...", a.totalCount(), b.totalCount());

		if (a.useSummarizers || b.useSummarizers)
		{
			System.err.printf("CpgWalkerAllpairsAutocorr::merge() does not work with useSummarizers set. Quitting.\n");
			System.exit(1);
		}
		
		// ***** DEBUGGING **** TURN OFF ****
		if (a.readType != b.readType)
		{
			System.err.printf("a.readType(%d) != b.readType(%d)\n",a.readType,b.readType);
			System.exit(1);
		}
		if (a.samestrandOnly != b.samestrandOnly)
		{
			System.err.printf("samestrandOnly(%s) != samestrandOnly(%s)\n",a.samestrandOnly,b.samestrandOnly);
			System.exit(1);
		}
		if (a.oppstrandOnly != b.oppstrandOnly)
		{
			System.err.printf("oppstrandOnly(%s) != oppstrandOnly(%s)\n",a.oppstrandOnly,b.oppstrandOnly);
			System.exit(1);
		}
		// END ***** DEBUGGING **** TURN OFF ****
		
		CpgWalkerAllpairsAutocorrByread out = new CpgWalkerAllpairsAutocorrByread(a.walkParams, a.samestrandOnly, a.oppstrandOnly, a.readType);
		int windSize = a.walkParams.maxScanningWindSize;
		for (int i = 0; i < (windSize-1); i++)
		{
			out.nMM[i] = a.nMM[i] + b.nMM[i];
			out.nMU[i] = a.nMU[i] + b.nMU[i];
			out.nUM[i] = a.nUM[i] + b.nUM[i];
			out.nUU[i] = a.nUU[i] + b.nUU[i];
			out.nM[i] = a.nM[i] + b.nM[i];
			out.nU[i] = a.nU[i] + b.nU[i];
		}
		
		// Reduce the feats
		out.aligner = a.aligner.mergeInAlignments(b.aligner, false);
		//aAligner.mergeInAlignments(bAligner, false)

		
		System.err.printf("Result has Autocorr(%d)\n", out.totalCount());
		return out;
	}

	
	
	protected void init()
	{
		// Initalize counters
		
		int windSize = this.walkParams.maxScanningWindSize;
		nMM = new int[windSize-1]; 
		nMU = new int[windSize-1]; 
		nUM = new int[windSize-1]; 
		nUU = new int[windSize-1]; 
		nM = new int[windSize-1]; 
		nU = new int[windSize-1]; 
	}

	
	
	
	public boolean useOnlyCG() {
		return useOnlyCG;
	}



	public void useOnlyCG(boolean useOnlyCG) {
		this.useOnlyCG = useOnlyCG;
	}

	protected int getPairDist(Cpg a, Cpg b)
	{
		return Math.abs(b.chromPos - a.chromPos) - 1;
	}

	
	
	@Override
	protected void recordPair(Cpg a, Cpg b)
	{
		// Get the distance
		
		
//		Map<Integer,CpgRead> aCgReads = a.getReads();
//		Map<Integer,CpgRead> bCgReads = b.getReads();
		
		Collection<CpgRead> aCgReads = a.getReads().values();
		Collection<CpgRead> bCgReads = b.getReads().values();
	
		
		// Check if we hit any feats
		AlignmentRelCoords pairAlignments = null;
		int midpoint = Math.round( ((float)a.chromPos + (float)b.chromPos)/2);
		if (feats != null)
		{
			pairAlignments = new AlignmentRelCoords(this.curChr, midpoint, false);
			pairAlignments.addAlignmentPoints(feats, aligner.flankSize);	
			//System.err.printf("Found %d pairAlignments (%d bp from midpoint %d)\n",pairAlignments.numAlignmentPoints(),aligner.flankSize,midpoint);
		}
		

		//System.err.println("Comparing CpG (A): " + a.toString());

		for (CpgRead aCgRead : aCgReads)
		{
			boolean aCg = aCgRead.validCg( (this.walkParams.methylParams==null) || this.walkParams.methylParams.getUseNonconversionFilter());
			//if (!aCg) System.err.println("Got an uncounted cytosine (A): " + a.toString() + "\n\t" + aCgRead.toString()); // COMMENT OUT
			if (!useOnlyCG || aCg)
			{
				int aReadId = aCgRead.readId;
				boolean aMeth = aCgRead.meth(true);
				
				// I found that having a wild imbalance between the number of reads checked for SAME_READ (1) and DIFFERENT_READ (many) could
				// cause a skew in the number of pairs because SAME_READ was undersampling in regions of low CpG content.  So to equalize,
				// we can choose to only take one random DIFFERENT read.
				if (useOnlyOneDifferentRead)
				{
					// Randomly permute so we don't introduce some bias by selecting the first one.
					Collections.shuffle(new ArrayList<CpgRead>(bCgReads));
				}

				B_READS: for (CpgRead bCgRead : bCgReads)
				{
					boolean bCg = bCgRead.validCg( (this.walkParams.methylParams==null) || this.walkParams.methylParams.getUseNonconversionFilter());
					//if (!bCg) System.err.println("Got uncounted cytosine (B): " + b.toString() + "\n\t" +  bCgRead.toString());  // COMMENT OUT
					if (!useOnlyCG || bCg)
					{
						int bReadId = bCgRead.readId;
						boolean identicalRead = (aReadId == bReadId);

						boolean include = (readType==ANY_READ) || (identicalRead&&(readType==SAME_READ)) || (!identicalRead&&(readType==DIFFERENT_READ));
						if (include)
						{
							boolean bMeth = bCgRead.meth(true);

//							// OLD symmetric way
//							if (aMeth) { nM[dist]++; } else { nU[dist]++; }
//							if (bMeth) { nM[dist]++; } else { nU[dist]++; }
//
//							if (aMeth && bMeth)
//							{
//								nMM[dist]+=2;
//							}
//							else if (!aMeth && !bMeth)
//							{
//								nUU[dist]+=2;
//							}
//							else
//							{
//								nMU[dist]++;
//							}

							// Now that it's asymmetric, we can the pair in either orientation.
							int dist = this.getPairDist(a, b);

//							System.err.printf("\tRecording %d,%d dist=%d\tameth=%s\tbmeth=%s\n", a.chromPos,b.chromPos,dist,aMeth,bMeth);

							// Add to Dist arrays
							if (dist<nM.length)
							{

								if (aMeth) { nM[dist]++; } else { nU[dist]++; }

								if (aMeth && bMeth)
								{
									nMM[dist]++;
								}
								else if (aMeth && !bMeth)
								{
									nMU[dist]++;
								}
								else if (!aMeth && bMeth)
								{
									nUM[dist]++;
								}
								else
								{
									nUU[dist]++;
								}							
							}

							// Add to feat aligner.
							// Check if we're the right kind
							if ((aMeth == this.firstMeth) && (bMeth == this.secondMeth))
							{
								if (pairAlignments != null)
								{
									Iterator<GenomicRangeWithRefpoint> it = pairAlignments.getAlignmentPoints();
									double score = 1.0;
									if (it != null)
									{
										while (it.hasNext())
										{
											GenomicRangeWithRefpoint ap = it.next();
											aligner.addAlignmentPos(
													midpoint,
													//												(read.getStrand() == StrandedFeature.NEGATIVE) ? Double.NaN : mLevel,
													//														(read.getStrand() == StrandedFeature.NEGATIVE) ? mLevel: Double.NaN,
													score, Double.NaN,
													"", ap.getChrom(), ap.getRefPoint(), ap.getStrand(), 0);
										}
									}
								}
							}
							
							// If we are using only a single different read, we are done with the B reads
							if (useOnlyOneDifferentRead  && !identicalRead)
							{
								break B_READS;
							}
						}  // End: if (include)
					}
				}
			}
		}
		

	}
	
	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalker#reset()
	 */
	@Override
	public void reset() {
		super.reset();
		init();
	}

	
	@Override
	public String headerStr()
	{
		StringBuffer sb = new StringBuffer((int)1E5);
		
		//sb.append(String.format("readType = %d,\tsamestrand=%s\n", this.readType, (this.samestrandOnly)?"true":"false"  ));
		sb.append(String.format("%s,%s,%s,%s,%s,%s,%s","Distance","nM","nU","nMM","nMU","nUM","nUU"));
		
		//System.err.println("Header="+sb.length());
		return sb.toString();
	}
	
	
	public String toCsvStrDumb()
	{
		StringBuffer sb = new StringBuffer((int)1E6);
		
		
//		ListUtils.setDelim(",");
//		for (Object key : counters.keySet())
//		{
//			sb.append(key.toString().replace(", ", "/"));
//			sb.append(',');
//
//			int[] counter = counters.get(key);
//			sb.append(ListUtils.excelLine(counter));
//			sb.append('\n');
//		}
		
		
		return sb.toString();
	}
	
	public String toCsvStr(String firstCol)
	throws Exception
	{
		StringBuffer sb = new StringBuffer((int)1E6);
		//ListUtils.setDelim(",");

		String firstColSec = (firstCol==null) ? "" : (firstCol + ",");
		
		for (int i = 0; i < nMM.length; i++)
		{
			sb.append(String.format("%s%d,%d,%d,%d,%d,%d,%d\n",firstColSec,i+1,nM[i],nU[i],nMM[i],nMU[i],nUM[i],nUU[i]));
		}
		
		return sb.toString();
	}
	
	@Override
	public String toCsvStr()
	throws Exception
	{
		return toCsvStr(null);
	}

	
	public int totalCount()
	{
		int[][] mat = new int[][]{nM, nU, nMM, nMU, nUM, nUU};
		return MatUtils.sumAll(mat);
	}


	

}
