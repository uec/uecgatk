/*
 * Copyright (c) 2011 USC Epigenome Center
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package edu.usc.epigenome.uecgatk.readmetrics;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import net.sf.samtools.SAMRecord;
import java.util.*;
import java.io.*;

/**
 * Zack Ramjan
 * USC Epigenome Center 
 * 03/23/2011
 */

/**
 * Walks over the input data set, counts the number of <cycle,base,quality> tuples
 */
@Requires({DataSource.READS})
public class BaseQualityWalker extends ReadWalker<HashMap<String,Long>,HashMap<String,Long>> {
    @Output
    protected PrintStream out;

    @Argument(fullName="mappedOnly", shortName="mo", doc="when this flag is set (default), statistics will be collected "+
                "on ALIGNED reads only, while unmapped reads will be discarded", required=false)
    protected boolean MAPPED_ONLY = false;
    
    @Argument(fullName="label",shortName="p",doc="label for the csv report",required=true)
    protected String LABEL = null;

    public void initialize() 
    {
    	 if ( LABEL == null ) throw new ReviewedStingException("Prefix for output file(s) must be specified");
    }


    public HashMap<String,Long> map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) 
    {
        if ( AlignmentUtils.isReadUnmapped(read) && MAPPED_ONLY) return null;
        
        byte [] quals = read.getBaseQualities();
        byte [] bases =  read.getReadBases();

        HashMap<String,Long> countRead = new HashMap<String,Long>();
        //postAlignment,a,+,8,27,156335
        for (short i = 0; i < bases.length; i++)
        	countRead.put((char)bases[i] + ",+," + (i+1) + "," + quals[i] , 1L);
        
        return countRead;
    }

    /**
     * Provide an initial value for reduce computations.
     *
     * @return Initial value of reduce.
     */
    public HashMap<String,Long> reduceInit() 
    {
        return  new HashMap<String,Long>();  //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * Reduces a single map with the accumulator provided as the ReduceType.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return accumulator with result of the map taken into account.
     */
    public HashMap<String,Long> reduce(HashMap<String,Long> value, HashMap<String,Long> sum) 
    {
    	for(String s : value.keySet())
    	{
    		if(sum.containsKey(s))
    			sum.put(s, sum.get(s) + 1L);
    		else
    			sum.put(s, 1L);
    	}
    	return sum;
    }

    public void onTraversalDone(HashMap<String,Long> result) 
    {
    	Vector<String> sortedKeys = new Vector<String>(result.keySet());
    	Collections.sort(sortedKeys);
    	for(String s : sortedKeys)
    		out.println(LABEL + "," + s + "," + result.get(s));    	
    }
 }
