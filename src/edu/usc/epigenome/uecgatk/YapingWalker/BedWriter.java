package edu.usc.epigenome.uecgatk.YapingWalker;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;

import org.broad.tribble.Tribble;
import org.broad.tribble.TribbleException;
import org.broad.tribble.index.DynamicIndexCreator;
import org.broad.tribble.index.Index;
import org.broad.tribble.index.IndexFactory;
import org.broad.tribble.util.LittleEndianOutputStream;
import org.broad.tribble.util.PositionalStream;
import org.broadinstitute.sting.utils.GenomeLoc;


public class BedWriter {


    protected BufferedWriter mWriter;
    protected PositionalStream positionalStream = null;
    protected DynamicIndexCreator indexer = null;
    protected File indexFile = null;
    LittleEndianOutputStream idxStream = null;
    File location = null;
	
	public BedWriter(File location) {
        this(location, openOutputStream(location));
    }
	
	public BedWriter(OutputStream output) {
        mWriter = new BufferedWriter(new OutputStreamWriter(output));
    }
	
	public BedWriter(File location, OutputStream output) {
        this.location = location;

       
                positionalStream = new PositionalStream(output);
                output = positionalStream;

        //mWriter = new BufferedWriter(new OutputStreamWriter(new PositionalStream(output)));
        mWriter = new BufferedWriter(new OutputStreamWriter(output));

    }
	
	private String locationString() {
        return location == null ? mWriter.toString() : location.getAbsolutePath();
    }

    /**
     * attempt to close the VCF file
     */
    public void close() {
        // try to close the vcf stream
        try {
            mWriter.flush();
            mWriter.close();
        } catch (IOException e) {
            throw new TribbleException("Unable to close " + locationString() + " because of " + e.getMessage());
        }

    }

    protected static OutputStream openOutputStream(File location) {
        try {
            return new FileOutputStream(location);
        } catch (FileNotFoundException e) {
            throw new TribbleException("Unable to create bed file at location: " + location);
        }
    }

    /**
     * add a record to the file
     *
     * @param loc      location
     * @param value 	the value for this location
     */
    public void add(GenomeLoc loc, double value) {
    	String bedLine = String.format("%s\t%d\t%d\t%.2f\n",loc.getContig(), loc.getStart(), loc.getStop(), value);
    	try {
			mWriter.write(bedLine);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    }
    
    public void add(String contig, long start, long end, double value) {
    	String bedLine = String.format("%s\t%d\t%d\t%.2f\n",contig, start, end, value);
    	try {
			mWriter.write(bedLine);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    }

}
