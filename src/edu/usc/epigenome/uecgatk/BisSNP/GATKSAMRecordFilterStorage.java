/**
 * 
 */
package edu.usc.epigenome.uecgatk.BisSNP;

import java.util.BitSet;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import edu.usc.epigenome.uecgatk.BisSNP.BadBaseFilterBisulfite;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time 2012 Mar 16, 2012 3:48:16 PM
 * 
 */
public class GATKSAMRecordFilterStorage {

	/**
	 * 
	 */

	private GATKSAMRecord GATKrecord = null;
	private boolean goodBase = false;
	private BisulfiteArgumentCollection BAC;
	
	public GATKSAMRecordFilterStorage(GATKSAMRecord GATKrecord, BisulfiteArgumentCollection BAC, int offset) {
		// TODO Auto-generated constructor stub
		this.GATKrecord = GATKrecord;
		this.BAC = BAC;
		setGoodBases(offset);
	}


    public boolean isGoodBase() {
        return goodBase;
    }
    
	private void setGoodBases(int offset) {
		byte[] quals = GATKrecord.getBaseQualities();
		
		if ( GATKrecord.getMappingQuality() >= BAC.MIN_MAPPING_QUALTY_SCORE && quals[offset] >= BAC.MIN_BASE_QUALTY_SCORE &&
	             (BAC.USE_BADLY_MATED_READS || (!BadMateFilter.hasBadMate(GATKrecord)) && !GATKrecord.getNotPrimaryAlignmentFlag()) ) {
	        	//System.out.println("bad mates");
				//if((GATKrecord.getReadPairedFlag() && GATKrecord.getProperPairFlag()))
				goodBase = true;
	     }

    }

	
}
