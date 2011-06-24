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
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.usckeck.genome.ChromFeatures;

import edu.usc.epigenome.genomeLibs.FeatAligners.FeatAligner;
import edu.usc.epigenome.genomeLibs.GenomicRange.GenomicRangeWithRefpoint;

public class CpgBackedByGatkWithAlignmentRelCoords extends CpgBackedByGatk {

	
	protected Set<GenomicRangeWithRefpoint> alignmentPoints = new HashSet<GenomicRangeWithRefpoint>(); // A list of alignmentPoints, one per feature within range of the CpG

	
	public CpgBackedByGatkWithAlignmentRelCoords() {
	}

	public static CpgBackedByGatkWithAlignmentRelCoords copy(CpgBackedByGatk superObj) {
		
		CpgBackedByGatkWithAlignmentRelCoords out = new CpgBackedByGatkWithAlignmentRelCoords(superObj.chromPos, superObj.negStrand, 
				superObj.totalReads, superObj.cReads, superObj.cReadsNonconversionFilt,	superObj.tReads, superObj.agReads, 
				superObj.totalReadsOpposite, superObj.aReadsOpposite, (int) superObj.getCpgWeight(), superObj.nextBaseGreads, 
				superObj.nextBaseTotalReads, superObj.getNextBaseRef());
		out.alignmentContext = superObj.alignmentContext;
		out.metaData = superObj.metaData;
		out.refContext = superObj.refContext;
		out.prevBaseRefUpperCase = superObj.getPrevBaseRef();
		
		return out;
	}

	public CpgBackedByGatkWithAlignmentRelCoords(int chromPos, boolean negStrand) {
		super(chromPos, negStrand);
	}

	public CpgBackedByGatkWithAlignmentRelCoords(int chromPos,
			boolean negStrand, AlignmentContext ac, RefMetaDataTracker meta,
			ReferenceContext rc) {
		super(chromPos, negStrand, ac, meta, rc);
	}

	public CpgBackedByGatkWithAlignmentRelCoords(int chromPos,
			boolean negStrand, short totalReads, short cReads,
			short cReadsNonconversionFilt, short tReads, short agReads,
			short totalReadsOpposite, short aReadsOpposite, int cpgWeight,
			short nextBaseGreads, short nextBaseTotalReads,
			char nextBaseRefUpperCase) {
		super(chromPos, negStrand, totalReads, cReads, cReadsNonconversionFilt,
				tReads, agReads, totalReadsOpposite, aReadsOpposite, cpgWeight,
				nextBaseGreads, nextBaseTotalReads, nextBaseRefUpperCase);
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
		int cPos = this.chromPos;
		Location posLoc = new RangeLocation(cPos-windowSize, cPos+windowSize); 

		int curChr = feats.chrom_from_public_str(this.getChrom());
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
		
		
//		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).fine(String.format(
//				"Fetching coords: censor=%s\talignToStart=%s\tchr=%s\tfeatS=%d\tfeatE=%d\tfeatStrand=%s\talignmentPoint=%d\tflankS=%d\tflankEnd=%d\t\n",
//				""+ censor, ""+ alignToStart, this.getChrom(), featS, featE, ""+featStrand, alignmentPoint, flankStart, flankEnd));
		
//		int cpgChromPos = this.chromPos;
//		StrandedFeature.Strand cpgStrand = this.getStrand();
//
//		if (extendRead > 0)
//		{
//			int strandMult = 0;
//			if (cpgStrand == StrandedFeature.POSITIVE)
//			{
//				strandMult = 1;
//			}
//			else if (cpgStrand == StrandedFeature.NEGATIVE)
//			{
//				strandMult = -1;
//			}
//			else
//			{
//				System.err.println("Why did we get  astrandless CpG?");
//			}
//			cpgChromPos += (  strandMult * extendRead );
//		}
		
	}
}
