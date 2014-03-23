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

package edu.usc.epigenome.uecgatk.qcmetrics.read;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;

import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;

import java.util.*;
import java.io.*;

/**
 * Zack Ramjan
 * USC Epigenome Center 
 * 03/23/2011
 */

/**
 * Walks over the input data set, counts the number of nmers of a specified length (size given as input on cmd line)
 */
@Requires({DataSource.READS})
public class CountNmerWalker extends ReadWalker<HashMap<String,Long>,HashMap<String,Long>> {
    @Output
    protected PrintStream out;

    @Argument(fullName="mappedOnly", shortName="mo", doc="when this flag is set (default), statistics will be collected "+
                "on ALIGNED reads only, while unmapped reads will be discarded", required=false)
    protected boolean MAPPED_ONLY = false;
    
    @Argument(fullName="nmer",shortName="p",doc="the \"n\" in nmer, ie length nmers you want",required=true)
    protected int NMER = 0;

    public void initialize() 
    {
    	 if ( NMER == 0 ) throw new ReviewedStingException("must specify an NMER greater than zero");
    }


    public HashMap<String,Long> map(ReferenceContext ref, GATKSAMRecord read, RefMetaDataTracker metaDataTracker) 
    {
        if ( AlignmentUtils.isReadUnmapped(read) && MAPPED_ONLY) return null;
        
        String bases =  new String(read.getReadBases());
        HashMap<String,Long> countRead = new HashMap<String,Long>();

        for (short i = 0; i < bases.length() - NMER; i++)
        {
        	String nmer = bases.substring(i, i+NMER);
        	countRead.put(nmer, 1L + (countRead.containsKey(nmer) ? countRead.get(nmer)  : 0));
        }
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
    			sum.put(s, sum.get(s) + value.get(s));
    		else
    			sum.put(s, value.get(s));
    	}
    	return sum;
    }

    public void onTraversalDone(HashMap<String,Long> result) 
    {
    	Vector<String> sortedKeys = new Vector<String>(result.keySet());
    	Collections.sort(sortedKeys);
    	for(String s : sortedKeys)
    		out.println(s + "\t" + result.get(s));    	
    }
 }
