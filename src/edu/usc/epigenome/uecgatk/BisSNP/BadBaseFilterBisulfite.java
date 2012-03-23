/**
 * 
 */
package edu.usc.epigenome.uecgatk.BisSNP;

import java.util.BitSet;

import net.sf.samtools.util.StringUtil;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import edu.usc.epigenome.uecgatk.BisSNP.BisulfiteEnums.MethylSNPModel;


/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time 2012 Mar 16, 2012 4:02:12 PM
 * This is copied from old GATK toolkit!!
 */
public class BadBaseFilterBisulfite{

	/**
	 * 
	 */
	private ReferenceContext refContext;
    private final BisulfiteArgumentCollection BAC;
    private final int MISMATCH_WINDOW_SIZE = 20;

    public BadBaseFilterBisulfite(ReferenceContext refContext, BisulfiteArgumentCollection BAC) {
        this.refContext = refContext;
        this.BAC = BAC;
    }

    public BitSet getGoodBases(final GATKSAMRecord record) {
        BitSet bitset = new BitSet(record.getReadLength());

        // if the mapping quality is too low or the mate is bad, we can just zero out the whole read and continue
        if ( record.getMappingQuality() < BAC.MIN_MAPPING_QUALTY_SCORE ||
             (!BAC.USE_BADLY_MATED_READS && BadMateFilter.hasBadMate(record)) ) {
        	//System.out.println("bad mates");
        	return bitset;
        }
        
        byte[] quals = record.getBaseQualities();
        for (int i = 0; i < quals.length; i++) {
            if ( quals[i] >= BAC.MIN_BASE_QUALTY_SCORE )
                bitset.set(i);
        }

        // if a read is too long for the reference context, extend the context (being sure not to extend past the end of the chromosome)
        if ( record.getAlignmentEnd() > refContext.getWindow().getStop() ) {
            GenomeLoc window = refContext.getWindow();
            byte[] bases = refContext.getBases();
            StringUtil.toUpperCase(bases);
            refContext = new ReferenceContext(refContext.getGenomeLocParser(),refContext.getLocus(), window, bases);
        }            

        BitSet mismatches;
        boolean paired = record.getReadPairedFlag();
        if(BAC.sequencingMode == MethylSNPModel.BM || BAC.sequencingMode == MethylSNPModel.GM){
        	mismatches = BisulfiteAlignmentUtils.mismatchesInRefWindow(record, refContext, BAC.MAX_MISMATCHES, MISMATCH_WINDOW_SIZE, BAC.sequencingMode, paired);
        }
        else{
        	mismatches = AlignmentUtils.mismatchesInRefWindow(record, refContext, BAC.MAX_MISMATCHES, MISMATCH_WINDOW_SIZE);
        }
      //  if ( mismatches != null )
      //      bitset.and(mismatches);

        return bitset;
    }

}
