package edu.usc.epigenome.uecgatk.YapingWriter;

import net.sf.samtools.SAMRecord;

import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time 2012 Feb 10, 2012 10:48:20 PM
 * 
 */
public class GetCytosineContext {

	/**
	 * 
	 */
	public GetCytosineContext() {
		// TODO Auto-generated constructor stub
		
	}
	
	public contextStatus getContext(ReadBackedPileup pileup, boolean paired, boolean useBadMate){
		return getContext(pileup, paired, useBadMate, 0, 0);
	}
	
	public contextStatus getContext(ReadBackedPileup pileup, boolean paired, int minDepth, int minCTDepth){
		return getContext(pileup, paired, false, 0, 0);
	}
	
	public contextStatus getContext(ReadBackedPileup pileup, boolean paired){
		return getContext(pileup, paired, false, 0, 0);
	}
	
	public contextStatus getContext(ReadBackedPileup pileup, boolean paired, boolean useBadMate, int minDepth, int minCTDepth){
		
		contextStatus num = new contextStatus();
		
		int numC = 0;
		int numT = 0;
		int numO = 0;
		for( PileupElement p : pileup){
			SAMRecord samRecord = p.getRead();
			if(samRecord.getDuplicateReadFlag()){ //get rid of duplicate reads
            	continue;
            }
        	int offset = p.getOffset();
        	if(offset < 0)//is deletion
        		continue;
        	if(paired){
        		try {
					samRecord = (SAMRecord) p.getRead().clone();
				} catch (CloneNotSupportedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
            	boolean Paired = samRecord.getReadPairedFlag();	
	        	boolean secondOfPair = samRecord.getSecondOfPairFlag();

	        	if (samRecord.getNotPrimaryAlignmentFlag())
				{
					continue;
				}
				
				// Inverted dups, count only one end
				if (samRecord.getAlignmentStart() == samRecord.getMateAlignmentStart() && samRecord.getReadNegativeStrandFlag() == samRecord.getMateNegativeStrandFlag())
				{
					if (samRecord.getSecondOfPairFlag()) continue;
   				}
	        	if (Paired  && !useBadMate && !samRecord.getProperPairFlag())
				{
					continue;
				}
	        	
	        	
	        	if(secondOfPair){	        		
		        	samRecord.setReadNegativeStrandFlag(!samRecord.getReadNegativeStrandFlag());        		
	        	}
	        	
        	}
			
        	boolean negStrand = samRecord.getReadNegativeStrandFlag();
		//	int alignmentS = samRecord.getAlignmentStart();
		//	int	onRefCoord = (negStrand) ? samRecord.getUnclippedEnd() : alignmentS;
			
			
			if(((GATKSAMRecord)p.getRead()).isGoodBase(p.getOffset())){
				if(negStrand){
					if(p.getBase()==BaseUtils.G){
						numC++;
					}
					else if(p.getBase()==BaseUtils.A){
						numT++;
					}
					else{
						numO++;
					}
					
				}
				else{
					if(p.getBase()==BaseUtils.C){
						numC++;
					}
					else if(p.getBase()==BaseUtils.T){
						numT++;
					}
					else{
						numO++;
					}
				}
			}
			
		}
		if((numC + numT) >= minCTDepth && (numC + numT + numO) >= minDepth){
			num.numC = numC;
			num.numT = numT;
			num.numOther = numO;
			return num;
		}
		else{
			return null;
		}
	}

	
	public class contextStatus{
		public int numC = 0;
		public int numT = 0;
		public int numOther = 0;
		public contextStatus(){
			
		}
	}
}
