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
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;

import java.io.*;

/**
 * Zack Ramjan
 * USC Epigenome Center 
 * 05/04/2011
 */

/**
 * Walks over the input data set, counts the number of <cycle,base,quality> tuples
 */
@Requires({DataSource.READS})
public class ReadCounterWalker extends ReadWalker <Long,Long> implements TreeReducible<Long> 
{
    @Output
    protected PrintStream out;

    @Argument(fullName="mappedOnly", shortName="mo", doc="when this flag is set (default), statistics will be collected "+
                "on ALIGNED reads only, while unmapped reads will be discarded", required=false)
    protected boolean MAPPED_ONLY = true;
    
   

    public void initialize() 
    {
    	
    }


    public Long map(ReferenceContext ref, GATKSAMRecord read, RefMetaDataTracker metaDataTracker) 
    {
    	return ( AlignmentUtils.isReadUnmapped(read) && MAPPED_ONLY) ? 0L : 1L;
    }

    /**
     * Provide an initial value for reduce computations.
     *
     * @return Initial value of reduce.
     */
    @Override
    public Long reduceInit() 
    { 
       	return 0L;
    }

    /**
     * Reduces a single map with the accumulator provided as the ReduceType.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return accumulator with result of the map taken into account.
     */
    @Override
    public Long reduce(Long value, Long sum) 
    {
    	return value + sum;
    	
    }
    
    @Override
	public Long treeReduce(Long lval, Long rval)
	{
    	return reduce(lval,rval);
	}

    public void onTraversalDone(Long result) 
    {
    	 out.println(result);
    		
    }
 }