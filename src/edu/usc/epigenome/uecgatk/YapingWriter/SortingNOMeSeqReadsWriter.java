/**
 * 
 */
package edu.usc.epigenome.uecgatk.YapingWriter;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time 2012 Mar 14, 2012 6:20:58 PM
 * 
 */
public class SortingNOMeSeqReadsWriter extends SortingFormatWriterBase {

	private int maxCachingStartDistance;
	/**
	 * @param readsWriter
	 * @param mAXIMUM_CACHE_FOR_OUTPUT_VCF
	 */
	public SortingNOMeSeqReadsWriter(FormatWriterBase innerWriter, int maxCachingStartDistance,
			boolean takeOwnershipOfInner) {
		super(innerWriter, takeOwnershipOfInner);
		this.maxCachingStartDistance = maxCachingStartDistance;
		// TODO Auto-generated constructor stub
	}

	public SortingNOMeSeqReadsWriter(FormatWriterBase readsWriter, int maxCachingStartDistance) {
		this((NOMeSeqReadsWriterImp)readsWriter, maxCachingStartDistance, false);
		// TODO Auto-generated constructor stub
	}

	
	protected void noteCurrentRecord(genomeObject obj) {
        super.noteCurrentRecord(obj); // first, check for errors

        // then, update mostUpstreamWritableLoc:
        int mostUpstreamWritableIndex = obj.getStart() - maxCachingStartDistance;
        this.mostUpstreamWritableLoc = Math.max(BEFORE_MOST_UPSTREAM_LOC, mostUpstreamWritableIndex);
    }

}
