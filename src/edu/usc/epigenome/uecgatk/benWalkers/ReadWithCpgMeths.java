package edu.usc.epigenome.uecgatk.benWalkers;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.StrandedFeature.Strand;
import org.usckeck.genome.ChromFeatures;

import edu.usc.epigenome.genomeLibs.IupacPatterns;
import edu.usc.epigenome.genomeLibs.FeatAligners.AlignmentRelCoords;
import edu.usc.epigenome.genomeLibs.GenomicRange.GenomicRangeWithRefpoint;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.uecgatk.FractionNonidentical;

public class ReadWithCpgMeths extends ArrayList<Cpg>{

	// Object variables
	protected StrandedFeature.Strand strand = null;
	protected AlignmentRelCoords alignments = null;
	protected String chrom = null;

	/**
	 * @param strand
	 */
	public ReadWithCpgMeths(Strand inStrand, String inChrom) {
		super();
		this.strand = inStrand;
		this.chrom = inChrom;
	}

	public static ReadWithCpgMeths copy(ReadWithCpgMeths old)
	{
		ReadWithCpgMeths out = new ReadWithCpgMeths(old.strand, old.chrom);
		out.addAll(old);
		if (old.alignments != null) out.alignments = AlignmentRelCoords.copy(old.alignments);
		return out;
	}
	
	
	public Map<String,Double> methLevels(IupacPatterns patterns)
	{
		Map<String,Double> out = new HashMap<String,Double>();
		Map<String, FractionNonidentical> orig = this.methLevelsFractions(patterns);
		
		for (String key : orig.keySet())
		{
			FractionNonidentical frac = orig.get(key);
			out.put(key, frac.doubleValue());
		}
		
		return out;
	}
	
	public Map<String,FractionNonidentical> methLevelsFractions(IupacPatterns patterns)
	{
		Map<String,FractionNonidentical> out = new HashMap<String,FractionNonidentical>();
		
		Iterator<Cpg> it = this.iterator();
		while (it.hasNext())
		{
			Cpg cpg = it.next();
			String contextOrig = cpg.getPrevBasesRef() + "C" + cpg.getNextBasesRef();
			String context = patterns.firstMatch(contextOrig);
//			System.err.printf("Read at pos %d, cpg at pos %d: contextOrig=%s, context=%s\n", this.midpointChromPos(), cpg.chromPos, contextOrig, context);
			if (context != null)
			{
				FractionNonidentical frac = out.get(context);
				if (frac == null)
				{
					frac = new FractionNonidentical();
				}
				frac.incDenominator();
				if (cpg.cReads > 0) frac.incNumerator();
				out.put(context, frac);
			}
		}
		
		return out;
	}

	public int midpointChromPos()
	{
		double min = Double.POSITIVE_INFINITY;
		double max = Double.NEGATIVE_INFINITY;
		
		Iterator<Cpg> it = this.iterator();
		while (it.hasNext())
		{
			Cpg cpg = it.next();
			min = Math.min(min, (double)cpg.chromPos);
			max = Math.max(max, (double)cpg.chromPos);
		}
		
		int out = (int)((min+max)/2);
		return out;
	}
	
	public String getChrom()
	{
		return this.chrom;
	}

	public Iterator<GenomicRangeWithRefpoint> getAlignmentPoints() {
		return (this.alignments==null) ? null : this.alignments.getAlignmentPoints();
	}

	/**
	 * @return the strand
	 */
	public StrandedFeature.Strand getStrand() {
		return strand;
	}

	/**
	 * @param strand the strand to set
	 */
	public void setStrand(StrandedFeature.Strand strand) {
		this.strand = strand;
	}


	// ****** ALIGNMENT points ********
	public void addAlignmentPoints(ChromFeatures feats, int windowSize)
	{
		if (this.alignments == null)
		{
			if (this.size() == 0)
			{
				System.err.printf("ReadWithCpgMeths::addAligmentPoints() called before Cpgs added.  Need cpg positions to add alignments");
				(new Exception()).printStackTrace();
				System.exit(1);
			}
			this.alignments = new AlignmentRelCoords(this.getChrom(), this.midpointChromPos(), (this.strand==StrandedFeature.NEGATIVE));
		}
		this.alignments.addAlignmentPoints(feats, windowSize);
	}

	
}
