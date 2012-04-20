package edu.usc.epigenome.uecgatk.BisSNP;


import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.TreeSet;

import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMRecord;

import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import edu.usc.epigenome.uecgatk.BisSNP.BisulfiteArgumentCollection;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Mar 20, 2012 3:39:46 PM
 * 
 */
public class DownsamplingBAM {

    private  int SAMPLE_READS_MEAN_COVERAGE = 30;
    private  BisulfiteArgumentCollection BAC;
    protected  SAMFileWriter samWriter = null;
	/**
	 * 
	 */
	public DownsamplingBAM(BisulfiteArgumentCollection BAC, SAMFileWriter samWriter) {
		// TODO Auto-generated constructor stub
		this.BAC = BAC;
		this.samWriter = samWriter;
	}

	
	 /*
     * we use different downsampling strategy. e.g. when downsampling 10X, it randomly pick up s*10/r reads.(r is the mean coverage of the sample, we use 30 here for our sample, 
     * s is the total reads covered in this position) 
     * it does NOT treat reads like GATK, when reads number more than 10, cut off the reads number to 10, and keep the same when reads number lower than 10. 
     */
    public void  downsamplingBamFile(AlignmentContext rawContext){
    	if(rawContext.hasReads()){
			String tag = "Xi";
			Integer coverageMarked = 0;
			//int covergaeLimit = getToolkit().getArguments().downsampleCoverage;
			int covergaeLimit = BAC.orcad;
			covergaeLimit = Math.max((covergaeLimit * rawContext.getBasePileup().depthOfCoverage())/SAMPLE_READS_MEAN_COVERAGE,1);
			//getToolkit().getArguments().downsampleCoverage = covergaeLimit;
			//covergaeLimit = (covergaeLimit * rawContext.getBasePileup().size())/SAMPLE_READS_MEAN_COVERAGE;
			//System.err.println("loc: " + rawContext.getLocation().getStart() + "\tcovergaeLimit: " + covergaeLimit + "\trawContext.getBasePileup().size(): " + rawContext.getBasePileup().size() + "\tdownsampleCoverage: " + getToolkit().getArguments().downsampleCoverage);
			
			//ReadBackedPileup downsampledPileup = rawContext.getBasePileup().getDownsampledPileup(covergaeLimit);
			ReadBackedPileup downsampledPileup = getDownsampledPileup(rawContext.getBasePileup(), covergaeLimit);
			//if(rawContext.getBasePileup().size() < covergaeLimit){
			//	downsampledPileup = rawContext.getBasePileup();
			//}
			//else{
			//	downsampledPileup = rawContext.getBasePileup().getDownsampledPileup(covergaeLimit);
			//}
			// = rawContext.getBasePileup().getDownsampledPileup(covergaeLimit);
			for ( PileupElement p : rawContext.getBasePileup() ) {
				if(p.getRead().getIntegerAttribute(tag) != null){
					if(p.getRead().getIntegerAttribute(tag) == 2)
						//System.out.println("loc: " + rawContext.getLocation().getStart() + " tag: " + p.getRead().getIntegerAttribute(tag));
					if(p.getRead().getIntegerAttribute(tag) == 1)
						coverageMarked++;
				}
					
			}
			//System.out.println("loc: " + rawContext.getLocation().getStart() + " coverageMarked: " + coverageMarked);
			for ( PileupElement p : downsampledPileup ) {
				//System.out.println(p.toString());
				if(p.getRead().getIntegerAttribute(tag) != null){
					//System.out.println("loc: " + rawContext.getLocation().getStart() + " tag: " + p.getRead().getIntegerAttribute(tag));
				}
				if(coverageMarked >= covergaeLimit)
					break;
				if(p.getRead().getIntegerAttribute(tag) == null){
					samWriter.addAlignment(p.getRead());
    				p.getRead().setAttribute(tag, 1);
    				coverageMarked++;
				}
				
					
			}
			for ( PileupElement p : rawContext.getBasePileup() ) {
				if(p.getRead().getIntegerAttribute(tag) == null)
					p.getRead().setAttribute(tag, 2);
			}
			
		}

    }
    
    
    public void downsamplingBamFileLikeGATK(AlignmentContext rawContext, GenomeAnalysisEngine toolkit){
    	if(rawContext.hasReads()){
			String tag = "Xi";
			Integer coverageMarked = 0;
			int covergaeLimit = toolkit.getArguments().downsampleCoverage;
			ReadBackedPileup downsampledPileup = rawContext.getBasePileup().getDownsampledPileup(covergaeLimit);
			for ( PileupElement p : rawContext.getBasePileup() ) {
				if(p.getRead().getIntegerAttribute(tag) != null){
					if(p.getRead().getIntegerAttribute(tag) == 2)
						System.out.println("loc: " + rawContext.getLocation().getStart() + " tag: " + p.getRead().getIntegerAttribute(tag));
					if(p.getRead().getIntegerAttribute(tag) == 1)
						coverageMarked++;
				}
					
			}
			//System.out.println("loc: " + rawContext.getLocation().getStart() + " coverageMarked: " + coverageMarked);
			for ( PileupElement p : downsampledPileup ) {
				//System.out.println(p.toString());
				if(coverageMarked >= covergaeLimit)
					break;
				if(p.getRead().getIntegerAttribute(tag) == null){
					samWriter.addAlignment(p.getRead());
    				p.getRead().setAttribute(tag, 1);
    				coverageMarked++;
				}
				
					
			}
			for ( PileupElement p : rawContext.getBasePileup() ) {
				if(p.getRead().getIntegerAttribute(tag) == null)
					p.getRead().setAttribute(tag, 2);
			}
			
		}

    }
	
    
    public ReadBackedPileup getDownsampledPileup(ReadBackedPileup pileup, int desiredCov){
		if ( pileup.depthOfCoverage() <= desiredCov )
            return pileup;

        // randomly choose numbers corresponding to positions in the reads list
        Random generator = new Random();
        TreeSet<Integer> positions = new TreeSet<Integer>();
        for ( int i = 0; i < desiredCov; /* no update */ ) {
            if ( positions.add(generator.nextInt(pileup.depthOfCoverage())) )
                i++;
        }
		GenomeLoc loc = pileup.getLocation();
		List<GATKSAMRecord> reads =  new ArrayList<GATKSAMRecord>();;
		List<Integer> elementOffsets = new ArrayList<Integer>();
		int i = 0;
		for ( PileupElement p : pileup ) {
			if(positions.contains(i)){
				int elementOffset = p.getOffset();
				if(elementOffset < 0 || elementOffset > p.getRead().getReadLength()-1)
					continue;
				elementOffsets.add(elementOffset);
				reads.add(p.getRead());
			}
				i++;
		}
		ReadBackedPileup downsampledPileup = new ReadBackedPileupImpl(loc,reads,elementOffsets);
		
		return downsampledPileup;
		
	}
    
}
