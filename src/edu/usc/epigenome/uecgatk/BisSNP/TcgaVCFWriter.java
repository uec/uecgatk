package edu.usc.epigenome.uecgatk.BisSNP;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.text.SimpleDateFormat;
import java.util.Calendar;

import net.sf.samtools.SAMSequenceDictionary;

import org.broad.tribble.TribbleException;
import org.broadinstitute.sting.utils.codecs.vcf.StandardVCFWriter;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.codecs.vcf.VCFFilterHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderVersion;

/*
 * Bis-SNP/BisSNP: It is a genotyping and methylation calling in bisulfite treated 
 * massively parallel sequencing (Bisulfite-seq and NOMe-seq) on Illumina platform
 * Copyright (C) <2011>  <Yaping Liu: lyping1986@gmail.com>

 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
 * VCF Writer to generate TCGA specific VCF file, it is only for sorted coordinate now, so in multi-thread mode, this could not be used
 */
public class TcgaVCFWriter extends StandardVCFWriter {
	//store reference file path and name
	protected String ref = null;
	
	public TcgaVCFWriter(File location,SAMSequenceDictionary refDict) {
		super(location, refDict);
		// TODO Auto-generated constructor stub
	}

	public TcgaVCFWriter(File location,SAMSequenceDictionary refDict, boolean enableOnTheFlyIndexing) {
		super(location,refDict, enableOnTheFlyIndexing);
		// TODO Auto-generated constructor stub
	}

	public TcgaVCFWriter(OutputStream output,SAMSequenceDictionary refDict) {
		super(output,refDict,false);
		// TODO Auto-generated constructor stub
	}

	public TcgaVCFWriter(OutputStream output,SAMSequenceDictionary refDict, boolean doNotWriteGenotypes) {
		super(output,refDict, doNotWriteGenotypes);
		// TODO Auto-generated constructor stub
	}

	public TcgaVCFWriter(File location, OutputStream output,SAMSequenceDictionary refDict,
			boolean enableOnTheFlyIndexing, boolean doNotWriteGenotypes) {
		super(location, output,refDict, enableOnTheFlyIndexing, doNotWriteGenotypes);
		// TODO Auto-generated constructor stub
	}

	@Override
	public void writeHeader(VCFHeader header) {
        mHeader = doNotWriteGenotypes ? new VCFHeader(header.getMetaData()) : header;
        String refGenomeversion = null;
        //System.err.println(ref);
        if(ref.contains("assembly18")){
        	refGenomeversion = "hg18";
        }
        else if(ref.contains("hg18")){
        	refGenomeversion = "hg18";
        }
        else if(ref.contains("hg19")){
        	refGenomeversion = "hg19";
        }
        else{
        	refGenomeversion = BisSNPUtils.getRefGenomeVersion();
        }
        try {
            // the file format field needs to be written first, specially for TCGA VCF header
            mWriter.write(VCFHeader.METADATA_INDICATOR + new VCFHeaderLine(BisulfiteVCFConstants.VCF_HEADER_VERSION_FORMAT,"VCFv4.1").toString() + "\n");
            mWriter.write(VCFHeader.METADATA_INDICATOR + new VCFHeaderLine(BisulfiteVCFConstants.VCF_HEADER_VERSION_DATE,now("yyyyMMdd")).toString() + "\n");
            mWriter.write(VCFHeader.METADATA_INDICATOR + new VCFHeaderLine(BisulfiteVCFConstants.VCF_HEADER_VERSION_TCGA_VERSION,"1.1").toString() + "\n");
            mWriter.write(VCFHeader.METADATA_INDICATOR + new VCFHeaderLine(BisulfiteVCFConstants.VCF_HEADER_VERSION_LOG, "<InputVCF=<>, InputVCFSource=<" + refGenomeversion + ">, InputVCFVer=<1.0>, InputVCFParam=<> InputVCFgeneAnno=<>>").toString() + "\n");
            mWriter.write(VCFHeader.METADATA_INDICATOR + new VCFHeaderLine(BisulfiteVCFConstants.VCF_HEADER_VERSION_REF,"<ID=" + refGenomeversion +",Source=" + ref + ">").toString() + "\n");
            mWriter.write(VCFHeader.METADATA_INDICATOR + new VCFHeaderLine(BisulfiteVCFConstants.VCF_HEADER_VERSION_ASSEMBLY,refGenomeversion).toString() + "\n");
            mWriter.write(VCFHeader.METADATA_INDICATOR + new VCFHeaderLine(BisulfiteVCFConstants.VCF_HEADER_VERSION_CENTER,"USC Epigenome Center").toString() + "\n");
            mWriter.write(VCFHeader.METADATA_INDICATOR + new VCFHeaderLine(BisulfiteVCFConstants.VCF_HEADER_VERSION_PHASE,"none").toString() + "\n");
            mWriter.write(VCFHeader.METADATA_INDICATOR + new VCFHeaderLine(BisulfiteVCFConstants.VCF_HEADER_VERSION_GAF,"none").toString() + "\n");
            for ( VCFHeaderLine line : mHeader.getMetaData() ) {
                if ( line.getKey().equals(VCFHeaderVersion.VCF4_0.getFormatString()) ||
                        line.getKey().equals(VCFHeaderVersion.VCF3_3.getFormatString()) ||
                        line.getKey().equals(VCFHeaderVersion.VCF3_2.getFormatString()))
                    continue;

                // are the records filtered (so we know what to put in the FILTER column of passing records) ?
                if ( line instanceof VCFFilterHeaderLine)
                    filtersWereAppliedToContext = true;

                mWriter.write(VCFHeader.METADATA_INDICATOR);
                mWriter.write(line.toString());
                mWriter.write("\n");
            }

            // write out the column line
            mWriter.write(VCFHeader.HEADER_INDICATOR);
            for ( VCFHeader.HEADER_FIELDS field : mHeader.getHeaderFields() ) {
                mWriter.write(field.toString());
                mWriter.write(VCFConstants.FIELD_SEPARATOR);
            }

            if ( mHeader.hasGenotypingData() ) {
                mWriter.write("FORMAT");
                for ( String sample : mHeader.getGenotypeSamples() ) {
                    mWriter.write(VCFConstants.FIELD_SEPARATOR);
                    mWriter.write(sample);
                }
            }

            mWriter.write("\n");
            mWriter.flush();  // necessary so that writing to an output stream will work
        }
        catch (IOException e) {
            throw new TribbleException("IOException writing the VCF header to " + e);
        }
    }
	
	//get the system time that this VCF file generated
	public static String now(String dateFormat) {
	    Calendar cal = Calendar.getInstance();
	    SimpleDateFormat sdf = new SimpleDateFormat(dateFormat);
	    return sdf.format(cal.getTime());

	  }
	
	public void setRefSource(String ref){
    	this.ref = ref;
    }
	
	public void writeFlush(){
		
		try {
			mWriter.flush();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			throw new TribbleException("IOException writing the VCF flush to " + e);
		}
		
	}
}
