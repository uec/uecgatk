/**
 * This is to implement Robert Shoemaker 2010 Genome Research paper's algorithm's comparison with BisSNP
 * discard reads in multiple location, bad mates were discarded.
 * T/C ratio bigger than A/G ratio is identified as bisulfite-conversion reads, then non-bisulfite-conversion read were converted to its reverse complement reads.
 * Cs in reads are converted to Ts, this is called in silico demethylation. 
 */
package edu.usc.epigenome.uecgatk.YapingWalker;

import java.util.ArrayList;
import java.util.List;
import java.util.MissingResourceException;
import java.util.ResourceBundle;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMProgramRecord;
import net.sf.samtools.SAMRecord;

import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.gatk.filters.DuplicateReadFilter;
import org.broadinstitute.sting.gatk.filters.MappingQualityUnavailableFilter;
import org.broadinstitute.sting.gatk.filters.NotPrimaryAlignmentFilter;
import org.broadinstitute.sting.gatk.filters.UnmappedReadFilter;
import org.broadinstitute.sting.gatk.io.StingSAMFileWriter;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.WalkerName;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.gatk.walkers.recalibration.Covariate;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.text.TextFormattingUtils;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Apr 29, 2012 2:48:50 PM
 * 
 */
@WalkerName("ShoemakerReadsConversion")
@Requires({DataSource.READS})
@ReadFilters( {UnmappedReadFilter.class, BadMateFilter.class, NotPrimaryAlignmentFilter.class} )
public class ShoemakerReadsConversionWalker extends ReadWalker<SAMRecord, SAMFileWriter> {

	public static final String PROGRAM_RECORD_NAME = "ShoemakerReadsConversionWalker";
	/**
	 * 
	 */
	
	 @Output(doc = "The output converted BAM file", required = true)
	    private StingSAMFileWriter OUTPUT_BAM = null;
	
	 public void initialize() {
		 final SAMFileHeader header = getToolkit().getSAMFileHeader().clone();
		 final SAMProgramRecord programRecord = new SAMProgramRecord(PROGRAM_RECORD_NAME);
         final ResourceBundle headerInfo = TextFormattingUtils.loadResourceBundle("StingText");
         try {
             final String version = headerInfo.getString("org.broadinstitute.sting.gatk.version");
             programRecord.setProgramVersion(version);
         } catch (MissingResourceException e) {
         }

         StringBuffer sb = new StringBuffer();
         sb.append(getToolkit().createApproximateCommandLineArgumentString(getToolkit(), this));
         
         programRecord.setCommandLine(sb.toString());

         List<SAMProgramRecord> oldRecords = header.getProgramRecords();
         List<SAMProgramRecord> newRecords = new ArrayList<SAMProgramRecord>(oldRecords.size() + 1);
         for (SAMProgramRecord record : oldRecords) {
             if (!record.getId().startsWith(PROGRAM_RECORD_NAME))
                 newRecords.add(record);
         }
         newRecords.add(programRecord);
         header.setProgramRecords(newRecords);

         // Write out the new header
         OUTPUT_BAM.writeHeader(header);
	 }

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.ReadWalker#map(org.broadinstitute.sting.gatk.contexts.ReferenceContext, org.broadinstitute.sting.utils.sam.GATKSAMRecord, org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker)
	 */
	@Override
	public SAMRecord map(ReferenceContext ref, GATKSAMRecord read,
			ReadMetaDataTracker metaDataTracker) {
		boolean isBsConvertedStrand = isBisulfiteConversionStrand(read);
		if(isBsConvertedStrand){
			read.setReadNegativeStrandFlag(!read.getReadNegativeStrandFlag());
		}
		GATKSAMRecord newRead = demethylatedReads(read);
		
		return newRead;
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.Walker#reduceInit()
	 */
	@Override
	public SAMFileWriter reduceInit() {
		// TODO Auto-generated method stub
		 return OUTPUT_BAM;
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.Walker#reduce(java.lang.Object, java.lang.Object)
	 */
	@Override
	public SAMFileWriter reduce(SAMRecord read, SAMFileWriter output) {
		if (output != null) {
            output.addAlignment(read);
        }
        return output;
	}

	public void onTraversalDone(SAMFileWriter output) {
	//	output.close();
		
	}
	
	private boolean isBisulfiteConversionStrand(GATKSAMRecord read){
		double tcRatio = basesRatioInReads(read, BaseUtils.T, BaseUtils.C);
		double agRatio = basesRatioInReads(read, BaseUtils.A, BaseUtils.G);
		return tcRatio > agRatio ? true: false;
		
	}
	
	private double basesRatioInReads(GATKSAMRecord read, byte a, byte b){
		byte[] bases = read.getReadNegativeStrandFlag() ? BaseUtils.simpleReverseComplement(read.getReadBases()) : read.getReadBases();
		int countBaseA = 0;
		int countBaseB = 0;
		for(byte base : bases){
			if(BaseUtils.basesAreEqual(a, base)){
				countBaseA++;
			}
			if(BaseUtils.basesAreEqual(b, base)){
				countBaseB++;
			}
		}
		return (double)countBaseA/(double)countBaseB;
	}
	
	private GATKSAMRecord demethylatedReads(GATKSAMRecord read){
		byte[] bases = read.getReadBases();
		if(read.getReadNegativeStrandFlag()){
			for(int i=0; i<bases.length; i++){
				if(BaseUtils.basesAreEqual(bases[i], BaseUtils.G)){ // in negative strand, gatk return 'c' as 'g'...
					bases[i] = BaseUtils.A;
				}
			}
		}
		else{
			for(int i=0; i<bases.length; i++){
				if(BaseUtils.basesAreEqual(bases[i], BaseUtils.C)){
					bases[i] = BaseUtils.T;
				}
			}
		}
		
		read.setReadBases(bases);
		return read;
	}
}
