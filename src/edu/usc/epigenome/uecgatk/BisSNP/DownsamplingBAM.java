/**
 * 
 */
package edu.usc.epigenome.uecgatk.BisSNP;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Mar 20, 2012 3:39:46 PM
 * 
 */
public class DownsamplingBAM {

	/**
	 * 
	 */
	public DownsamplingBAM() {
		// TODO Auto-generated constructor stub
	}

	
	 /*
     * we use different downsampling strategy. e.g. when downsampling 10X, it randomly pick up s*10/r reads.(r is the mean coverage of the sample, we use 30 here for our sample, 
     * s is the total reads covered in this position) 
     * it does NOT treat reads like GATK, when reads number more than 10, cut off the reads number to 10, and keep the same when reads number lower than 10. 
     */
    public void downsamplingBamFile(AlignmentContext rawContext){
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
			ReadBackedPileup downsampledPileup = BisSNPUtils.getDownsampledPileup(rawContext.getBasePileup(), covergaeLimit);
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
    
    public void downsamplingBamFileLikeGATK(AlignmentContext rawContext){
    	if(rawContext.hasReads()){
			String tag = "Xi";
			Integer coverageMarked = 0;
			int covergaeLimit = getToolkit().getArguments().downsampleCoverage;
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
	
}
