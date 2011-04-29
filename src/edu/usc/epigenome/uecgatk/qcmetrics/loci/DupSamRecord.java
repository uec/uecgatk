package edu.usc.epigenome.uecgatk.qcmetrics.loci;
import  net.sf.samtools.SAMRecord;

public class DupSamRecord extends SAMRecord
{
	SAMRecord read1;
	public DupSamRecord(SAMRecord read)
	{
		super(read.getHeader());
		read1 = read;
	}
	
	@Override
	public boolean equals(Object obj)
	{
		if(obj instanceof DupSamRecord)
		{
			DupSamRecord other = (DupSamRecord) obj;
			return (
					this.read1.getAlignmentStart() == other.getAlignmentStart() &&
					this.read1.getAlignmentEnd() == other.getAlignmentEnd() &&
					this.read1.getMateAlignmentStart() == other.getMateAlignmentStart() &&
					this.read1.getReferenceIndex() == other.getReferenceIndex() &&
					this.read1.getMateReferenceIndex() == other.getMateReferenceIndex() 
					);
				
		}
		return false;
		
	}
	
	@Override
	public int hashCode()
	{
		return (read1.getAlignmentStart() + read1.getAlignmentEnd() + read1.getMateAlignmentStart());		
	}
}
