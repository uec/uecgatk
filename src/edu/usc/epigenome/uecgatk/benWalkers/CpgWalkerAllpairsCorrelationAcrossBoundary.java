package edu.usc.epigenome.uecgatk.benWalkers;

import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.PointLocation;
import org.usckeck.genome.ChromFeatures;
import org.usckeck.genome.GFFUtils;

import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerParams;
import edu.usc.epigenome.uecgatk.benWalkers.cytosineWalkers.Cpg;

public class CpgWalkerAllpairsCorrelationAcrossBoundary extends
		CpgWalkerAllpairsAutocorrByreadWcontext {

	protected ChromFeatures feats = null;
	private int chromInt = 0;
	private int boundaryBuffer = 0;
	
	public CpgWalkerAllpairsCorrelationAcrossBoundary(
			CpgWalkerParams inWalkParams, boolean inSamestrandOnly,
			boolean inOppstrandOnly, boolean inSameRead,
			boolean inDifferentRead, String inFromContext, String inToContext,
			ChromFeatures inFeats, int inBoundaryBuffer) {
		super(inWalkParams, inSamestrandOnly, inOppstrandOnly, inSameRead,
				inDifferentRead, inFromContext, inToContext);
		this.feats = inFeats;
		this.boundaryBuffer = inBoundaryBuffer;
		
		if (this.boundaryBuffer > 0)
		{
			System.err.printf("Error, CpgWalkerAllpairsCorrelationAcrossBoundary can not handle boundary buffer > 0 yet. Quitting\n");
			System.exit(1);
		}
	}
	
	
	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalker#alertNewChrom()
	 */
	@Override
	protected void alertNewChrom() {
		super.alertNewChrom();

		this.chromInt = (new ChromFeatures()).chrom_from_public_str(this.curChr);
	}


	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerAllpairsAutocorrByreadWcontext#recordPair(edu.usc.epigenome.genomeLibs.MethylDb.Cpg, edu.usc.epigenome.genomeLibs.MethylDb.Cpg)
	 */
	@Override
	protected void recordPair(Cpg a, Cpg b) {
		Location aloc = new PointLocation(a.chromPos);
		Location bloc = new PointLocation(b.chromPos);

		boolean aov = feats.overlaps_loc(this.chromInt, aloc);
		boolean bov = feats.overlaps_loc(this.chromInt, bloc);

		// One must overlap and one not
		if ((aov && !bov) || (!aov && bov))	super.recordPair(a, b);
	}


	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerAllpairsAutocorrByread#getPairDist(edu.usc.epigenome.genomeLibs.MethylDb.Cpg, edu.usc.epigenome.genomeLibs.MethylDb.Cpg)
	 */
	@Override
	protected int getPairDist(Cpg a, Cpg b) {
		
		
		Location aloc = new PointLocation(a.chromPos);
		Location bloc = new PointLocation(b.chromPos);

		// One should overlap and one not.  
		boolean aov = feats.overlaps_loc(this.chromInt, aloc);
		boolean bov = feats.overlaps_loc(this.chromInt, bloc);
		Location ovloc = null;
		Location novloc = null;
		if (aov && !bov)
		{
			ovloc = aloc;
			novloc = bloc;
		}
		else if (!aov && bov)
		{
			ovloc = bloc;
			novloc = aloc;
		}
		else
		{
			System.err.printf("Recording pair where a does not overlap features, or b does. Killing.\n");
			System.exit(1);
		}

		int elementDist = 0;
		GFFRecord ovrec = feats.closest_feature(this.chromInt, ovloc, 0, false);
		if (ovloc.getMin() > novloc.getMin())
		{
			elementDist = ovrec.getStart() - novloc.getMin();
		}
		else
		{
			elementDist = novloc.getMin() - ovrec.getEnd();
		}
		
		int realDist = Math.abs(ovloc.getMin()-novloc.getMin());
		//System.err.printf("PAIR\tovpos=%d\tnovpos=%d\trealdist=%d\telementdist=%d\trec=%s\n",ovloc.getMin(),novloc.getMin(),realDist,elementDist,GFFUtils.gffBetterString(ovrec));
		
		
		return elementDist;
	}
	
	

}
