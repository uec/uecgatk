package edu.usc.epigenome.uecgatk.YapingWriter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;

import org.broad.tribble.TribbleException;
import org.broad.tribble.index.DynamicIndexCreator;
import org.broad.tribble.util.LittleEndianOutputStream;
import org.broad.tribble.util.PositionalStream;
import org.broadinstitute.sting.utils.GenomeLoc;

public abstract class FormatWriterBase {
	protected BufferedWriter mWriter;
    protected PositionalStream positionalStream = null;
    protected DynamicIndexCreator indexer = null;
    protected File indexFile = null;
    LittleEndianOutputStream idxStream = null;
    File location = null;
	
	public FormatWriterBase(File location) {
        this(location, openOutputStream(location));
    }
	
	public FormatWriterBase(OutputStream output) {
        mWriter = new BufferedWriter(new OutputStreamWriter(output));
    }

	public FormatWriterBase(File location, OutputStream output) {
        this.location = location;

       
                positionalStream = new PositionalStream(output);
                output = positionalStream;

        //mWriter = new BufferedWriter(new OutputStreamWriter(new PositionalStream(output)));
        mWriter = new BufferedWriter(new OutputStreamWriter(output));

    }
	
	private String locationString() {
        return location == null ? mWriter.toString() : location.getAbsolutePath();
    }


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


    public abstract void add(genomeObject obj);
    
    
    public abstract void addHeader(Object o);
    
}
