package edu.usc.epigenome.uecgatk.YapingWriter;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFWriter;

public class SortingVCFWriterOwn extends SortingVCFWriterOwnBase {

	// the maximum START distance between records that we'll cache
    private int maxCachingStartDistance;

    /**
     * create a local-sorting VCF writer, given an inner VCF writer to write to
     *
     * @param innerWriter        the VCFWriter to write to
     * @param maxCachingStartDistance the maximum start distance between records that we'll cache
     * @param takeOwnershipOfInner Should this Writer close innerWriter when it's done with it
     */
    public SortingVCFWriterOwn(VCFWriter innerWriter, int maxCachingStartDistance, boolean takeOwnershipOfInner) {
        super(innerWriter, takeOwnershipOfInner);
        this.maxCachingStartDistance = maxCachingStartDistance;
    }

    public SortingVCFWriterOwn(VCFWriter innerWriter, int maxCachingStartDistance) {
        this(innerWriter, maxCachingStartDistance, false); // by default, don't own inner
    }

    protected void noteCurrentRecord(VariantContext vc) {
        super.noteCurrentRecord(vc); // first, check for errors

        // then, update mostUpstreamWritableLoc:
        int mostUpstreamWritableIndex = vc.getStart() - maxCachingStartDistance;
        this.mostUpstreamWritableLoc = Math.max(BEFORE_MOST_UPSTREAM_LOC, mostUpstreamWritableIndex);
    }
    
    public void enableDiscreteLoci(boolean enableDiscreteLoci) {
    	this.enableDiscreteLoci = enableDiscreteLoci;
    }

}
