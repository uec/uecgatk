package edu.usc.epigenome.uecgatk.bisulfitesnpmodel;

import org.broad.tribble.vcf.SortingVCFWriter;
import org.broad.tribble.vcf.VCFHeader;
import org.broad.tribble.vcf.VCFWriter;

public class SortingTcgaVCFWriter extends SortingVCFWriter {

	protected TcgaVCFWriter tcgaInnerWriter = null;
	public SortingTcgaVCFWriter(TcgaVCFWriter innerWriter,
			int maxCachingStartDistance, boolean takeOwnershipOfInner) {
		super(innerWriter, maxCachingStartDistance, takeOwnershipOfInner);
		tcgaInnerWriter = innerWriter;
		// TODO Auto-generated constructor stub
	}

	public SortingTcgaVCFWriter(TcgaVCFWriter innerWriter,
			int maxCachingStartDistance) {
		super(innerWriter, maxCachingStartDistance);
		tcgaInnerWriter = innerWriter;
		// TODO Auto-generated constructor stub
	}

	public TcgaVCFWriter getInnerWriter(){
		return this.tcgaInnerWriter;
	}
	
}
