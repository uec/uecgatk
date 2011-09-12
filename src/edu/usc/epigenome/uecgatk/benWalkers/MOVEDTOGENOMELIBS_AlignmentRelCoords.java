package edu.usc.epigenome.uecgatk.benWalkers;

import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import org.biojava.bio.program.gff.GFFEntrySet;
import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.program.gff.SimpleGFFRecord;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.RangeLocation;
import org.usckeck.genome.ChromFeatures;

import edu.usc.epigenome.genomeLibs.FeatAligners.FeatAligner;
import edu.usc.epigenome.genomeLibs.GenomicRange.GenomicRangeWithRefpoint;

public class MOVEDTOGENOMELIBS_AlignmentRelCoords {

	
	protected Set<GenomicRangeWithRefpoint> alignmentPoints = new HashSet<GenomicRangeWithRefpoint>(); // A list of alignmentPoints, one per feature within range of the CpG
	public String chr;
	public int elementChromPos;
	public boolean negStrand;
	
	
	public MOVEDTOGENOMELIBS_AlignmentRelCoords(String inChrom, int inElementChromPos, boolean inNegStrand)
	{
		this.chr = inChrom;
		this.elementChromPos = inElementChromPos;
		this.negStrand = inNegStrand;
	}

	public static MOVEDTOGENOMELIBS_AlignmentRelCoords copy(MOVEDTOGENOMELIBS_AlignmentRelCoords orig)
	{
		MOVEDTOGENOMELIBS_AlignmentRelCoords out = new MOVEDTOGENOMELIBS_AlignmentRelCoords(orig.chr, orig.elementChromPos, orig.negStrand);
		out.alignmentPoints.addAll(orig.alignmentPoints);
		return out;
	}


	public void resetAlignmentPoints()
	{
		//relativeCoords = new HashMap<GenomicRange,Integer>();
		alignmentPoints = new HashSet<GenomicRangeWithRefpoint>();
	}

	
	public Iterator<GenomicRangeWithRefpoint> getAlignmentPoints()
	{
		return alignmentPoints.iterator();	
	}
	
	public void addAlignmentPoint(GenomicRangeWithRefpoint element)
	{
		this.alignmentPoints.add(element);
	}
	
	public void addAlignmentPoints(ChromFeatures feats)
	{
		this.addAlignmentPoints(feats, 0);
	}
	
	
	public void addAlignmentPoints(ChromFeatures feats, int windowSize)
	{
		int cPos = this.elementChromPos;
		Location posLoc = new RangeLocation(cPos-windowSize, cPos+windowSize); 

		int curChr = feats.chrom_from_public_str(this.chr);
		GFFEntrySet ovs = feats.coord_filtered_features(curChr, posLoc, false); // Not reentrant, can fail in threaded use
		
		int nOvs = ovs.size();
		
		if (nOvs!=1)
		{
//			System.err.printf("processCytosine(%s,%d) found %d overlapping feats\n",
//					this.getChrom(),this.chromPos,nOvs);
		}
		if (nOvs==0)
		{
//			System.err.printf("Why can't we find gff record for position %s\n",posLoc);
		}
		
		Iterator it = ovs.lineIterator();
		while (it.hasNext())
		{
			GFFRecord rec = (GFFRecord)it.next();
			this.addRelativeCoords((SimpleGFFRecord)rec, windowSize);
		}
	}
	
	
	
	public void addRelativeCoords(SimpleGFFRecord rec, int flankSize)
	{
		// Code taken from MethylDbToMultisampleFeatAlignmentsStratified::processChrom
		
		StrandedFeature.Strand featStrand = rec.getStrand();
		String featName = null; // rec.getSeqName();
		int featS = rec.getStart();
		int featE = rec.getEnd();
		
		boolean skipUnoriented = true;
		boolean alignToStart = true;
		boolean alignToEnd = false;
		boolean censor = false;
		int extendRead = 0;
		if (skipUnoriented)
		{
			// Don't use those without orientation
			if (featStrand == StrandedFeature.UNKNOWN) return;
		}
		else
		{
			if (featStrand == StrandedFeature.UNKNOWN) rec.setStrand(StrandedFeature.POSITIVE);
			featStrand = rec.getStrand();
		}
		
		
		GenomicRangeWithRefpoint flankRange = FeatAligner.getAlignmentpointAndFlank(rec, 
				flankSize, alignToStart, alignToEnd, censor);
		flankRange.setStart(flankRange.getStart() - extendRead);
		flankRange.setEnd(flankRange.getEnd() + extendRead);
		flankRange.setStrand(featStrand);
		
		this.addAlignmentPoint(flankRange);
		
		
	}
}
