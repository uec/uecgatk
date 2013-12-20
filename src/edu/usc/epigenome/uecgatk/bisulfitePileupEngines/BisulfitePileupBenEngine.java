package edu.usc.epigenome.uecgatk.bisulfitePileupEngines;

import java.io.PrintStream;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.collections.Pair;

import edu.usc.epigenome.uecgatk.BaseUtilsMore;
import edu.usc.epigenome.uecgatk.benWalkers.CpgBackedByGatk;
import edu.usc.epigenome.uecgatk.pileup.BisulfitePileup;

public class BisulfitePileupBenEngine implements BisulfitePileupEngine {



	@Argument(fullName = "outputCph", shortName = "cph", doc = "Output CpHs in addition to Cpgs", required = false)
	public boolean outputCph = false;

	@Argument(fullName = "minCT", shortName = "mct", doc = "Minimum number of C/T bases to process cytosine (default=1)", required = false)
	public int minCT = 1;

	@Argument(fullName = "maxOppAfrac", shortName = "maxa", doc = "Maximum fraction of opposite strand reads being A (default=0.101)", required = false)
	public double maxOppAfrac = 0.101;

	@Argument(fullName = "minNextGfrac", shortName = "minnextg", doc = "Minimum fraction of next position reads being G, only if outputCph is false (default=0.899)", required = false)
	public double minNextGfrac = 0.899;

	@Argument(fullName = "minContextFracReadsMatching", shortName = "contextminmatch", doc = "Minimum fraction of reads matching to call CpH or CpG context (default=0.899)", required = false)
	public double minContextFracReadsMatching = 0.899;

	@Argument(fullName = "minContextPhredQuality", shortName = "contextminphred", doc = "Minimum phred quality for flanking bases to determine context (default=0)", required = false)
	public int minContextPhredQuality = 0;

	@Argument(fullName = "minConv", shortName = "minc", doc = "minimum number of converted cytosines required for 5' conversion filter (default=1)", required = false)
	public int minConv = 1;

	@Argument(fullName = "upstreamNumContextBases", shortName = "upcontext", doc = "Number of bases of context upstream of the cytosine (default=1)", required = false)
	public int upstreamNumContextBases = 1;

	@Argument(fullName = "downstreamNumContextBases", shortName = "downcontext", doc = "Number of bases of context downstream of the cytosine (default=1)", required = false)
	public int downstreamNumContextBases = 1;


	// private state variables
	protected String prevContig = null;
    protected static Logger logger = Logger.getLogger(Walker.class);


	@Override
	public Boolean supportsMultithreadedMode() {
		return false;
	}


	@Override
	public BisulfitePileup gatkLocusToBisulfitePileup(
			RefMetaDataTracker tracker, ReferenceContext ref,
			AlignmentContext context, GenomeAnalysisEngine toolkit) {

		BisulfitePileup out = null;
		
		// Are we on a new chrom?
		GenomeLoc thisLoc = ref.getLocus();
		int centerCoord = thisLoc.getStart();
		String thisContig = ref.getLocus().getContig();



		// Increment
		prevContig = thisContig;

		// *** NOTE :: JUST USING POSITIONS THAT ARE C IN REFERENCE :: NEED TO CHECK SAMPLE GENOME ****
		boolean isC = false;
		boolean negStrand = false;
		if (ref.getBase() == BaseUtils.Base.C.base)
		{
			isC = true;
			negStrand = false;
		}
		else if (ref.getBase() == BaseUtils.Base.G.base)
		{
			isC = true;
			negStrand = true;
		}

		if (isC)
		{

			byte[] contextSeq = null;
			try
			{
				// I've had this throw an out of range exception from time to time
				contextSeq = 
						toolkit.getReferenceDataSource().getReference().getSubsequenceAt(thisContig, centerCoord-1, centerCoord+1).getBases();
			}
			catch (Exception e)
			{
				System.err.printf("Non-fatal error, could not getToolkit().getReferenceDataSource().getReference().getSubsequenceAt(%s,%d,%d)\n%s\n",
						thisContig, centerCoord-1, centerCoord+1,e.toString());
				e.printStackTrace();
				return out;
			}


			contextSeq = BaseUtilsMore.toUpperCase(contextSeq);
			byte[] contextSeqStranded = contextSeq;
			if (negStrand) contextSeqStranded = BaseUtils.simpleReverseComplement(contextSeq);

			byte[] contextSeqStrandedIupac = new byte[contextSeqStranded.length];
			for (int i = 0; i < contextSeqStranded.length; i++)
			{
				if (i == 1) // Exclude central base, always a C
				{
					contextSeqStrandedIupac[i] = contextSeqStranded[i];
				}
				else
				{
					if(contextSeqStranded[i] == BaseUtils.Base.A.base || contextSeqStranded[i] == BaseUtils.Base.C.base || contextSeqStranded[i] == BaseUtils.Base.T.base)
						contextSeqStrandedIupac[i] = (byte)'H';
					else
						contextSeqStrandedIupac[i] = contextSeqStranded[i];
					
//					switch (contextSeqStranded[i])
//					{
//					case BaseUtils.Base.A.base:
//					case BaseUtils.Base.C.base:
//					case BaseUtils.Base.T.base:
//						contextSeqStrandedIupac[i] = (byte)'H';
//						break;
//					default:
//						contextSeqStrandedIupac[i] = contextSeqStranded[i];
//						break;		
//					}
				}
			}


			//    		logger.info(String.format("pos=%d\tcontext(%d,%d)=%s (strand=%d)",centerCoord,0,0,new String(contextSeqStrandedIupac),
			//    				(negStrand?-1:1)));

			// Make the bisulfitePileup
			out = null;
////			Pair<CpgBackedByGatk,MapType> makeCytPair = makeCytosine(thisLoc, ref, contextSeqStrandedIupac,contextSeqStranded, negStrand,context,tracker);
////
//			
//			CpgBackedByGatk thisC = makeCytPair.first;
//			MapType reduction = makeCytPair.second;
//
//			//			out.printf("MAP(): %d\t%s\t%s\t%s\t%d\t%s\t%s\n", centerCoord,new String(ref.getBases()),
//			//					new String(contextSeqStranded),new String(contextSeqStrandedIupac),(negStrand?-1:1),thisC.toStringExpanded(), thisC.context());
//
//			if (reduction != null)
//			{
//				mapout = reduction;
//
//			}
//			else if (this.outputCph || !thisC.isCph(false, 0.101))
//			{
//
//
//				// And process it
//				boolean minCTpasses = ((thisC.cReads + thisC.tReads) >= this.minCT);
//				double fracOppA = thisC.fracOppositeA();
//				if (Double.isNaN(fracOppA)) fracOppA = 0.0;
//				double fracNextG = thisC.fracNextBaseG();
//				if (Double.isNaN(fracNextG)) fracNextG = 1.0;
//				boolean maxOppAfracPasses = (fracOppA <= maxOppAfrac);
//				boolean minNextGfracPasses = (this.outputCph || (fracNextG>=this.minNextGfrac));
//
//				//				logger.info(String.format("Processing cytosine=%s\tCT=%d\toppA=%.2f\n",minCTpasses&&maxOppAfracPasses,(thisC.cReads+thisC.tReads),fracOppA));
//
//				if (minCTpasses && maxOppAfracPasses && minNextGfracPasses)
//				{
//					// it's a good one.  return it.
//					out = null;
//				}
//			}


		}
		return out;
	}

}
