package edu.usc.epigenome.uecgatk.qcmetrics.loci;
import  net.sf.samtools.SAMRecord;

public class DupSamRecord extends SAMRecord
{
	private SAMRecord read;
	public SAMRecord getRead()
	{
		return read;
	}

	public DupSamRecord(SAMRecord readIn)
	{
		super(readIn.getHeader());
		read = readIn;	
	}
	
	@Override
	public boolean equals(Object obj)
	{
		if(obj instanceof DupSamRecord)
		{
			DupSamRecord other = (DupSamRecord) obj;
			//SR
			if(this.getRead().getReadPairedFlag() == other.getRead().getReadPairedFlag() == false)
				return (
						this.getRead().getAlignmentStart() == other.getRead().getAlignmentStart() &&
						this.getRead().getReferenceIndex() == other.getRead().getReferenceIndex() 
						);
			
			//PE
			if(this.getRead().getReadPairedFlag() == other.getRead().getReadPairedFlag() == true)
				return (
						this.getRead().getAlignmentStart() == other.getRead().getAlignmentStart() &&
						this.getRead().getMateAlignmentStart() == other.getRead().getMateAlignmentStart() &&
						this.getRead().getReferenceIndex() == other.getRead().getReferenceIndex() &&
						this.getRead().getMateReferenceIndex() == other.getRead().getMateReferenceIndex() 
						);
					
		}
		return false;
		
	}
	
	@Override
	public int hashCode()
	{
		//PE
		if(this.getRead().getReadPairedFlag() == true)
			return (read.getAlignmentStart() + read.getMateAlignmentStart());
		//SR
		else
			return (read.getAlignmentStart());
	}
}
