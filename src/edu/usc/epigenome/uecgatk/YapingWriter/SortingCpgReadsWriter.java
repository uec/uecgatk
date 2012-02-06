package edu.usc.epigenome.uecgatk.YapingWriter;


public class SortingCpgReadsWriter extends SortingFormatWriterBase {

	private int maxCachingStartDistance;
	
	public SortingCpgReadsWriter(cpgReadsWriterImp innerWriter, int maxCachingStartDistance,
			boolean takeOwnershipOfInner) {
		super(innerWriter, takeOwnershipOfInner);
		this.maxCachingStartDistance = maxCachingStartDistance;
		// TODO Auto-generated constructor stub
	}

	public SortingCpgReadsWriter(cpgReadsWriterImp innerWriter, int maxCachingStartDistance) {
		this(innerWriter, maxCachingStartDistance, false);
		// TODO Auto-generated constructor stub
	}


    protected void noteCurrentRecord(genomeObject obj) {
        super.noteCurrentRecord(obj); // first, check for errors

        // then, update mostUpstreamWritableLoc:
        int mostUpstreamWritableIndex = obj.getStart() - maxCachingStartDistance;
        this.mostUpstreamWritableLoc = Math.max(BEFORE_MOST_UPSTREAM_LOC, mostUpstreamWritableIndex);
    }
	
}
