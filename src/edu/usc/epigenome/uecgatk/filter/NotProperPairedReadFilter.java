/**
 * get rid of Not proper paired reads
 */
package edu.usc.epigenome.uecgatk.filter;

import net.sf.samtools.SAMRecord;

import org.broadinstitute.sting.gatk.filters.ReadFilter;


/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Dec 11, 2012 11:14:32 AM
 * 
 */
public class NotProperPairedReadFilter extends ReadFilter {
	public boolean filterOut( final SAMRecord read ) {
        return read.getReadPairedFlag() && !read.getProperPairFlag();
    }

	/* (non-Javadoc)
	 * @see net.sf.picard.filter.SamRecordFilter#filterOut(net.sf.samtools.SAMRecord, net.sf.samtools.SAMRecord)
	 */
	@Override
	public boolean filterOut(SAMRecord arg0, SAMRecord arg1) {
		// TODO Auto-generated method stub
		return false;
	}

}
