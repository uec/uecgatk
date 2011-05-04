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
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import net.sf.samtools.SAMRecord;
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
public class InvertedReadPairDupsWalker extends ReadWalker <Long[],Long[]> implements TreeReducible<Long[]> 
{
    @Output
    protected PrintStream out;

    @Argument(fullName="mappedOnly", shortName="mo", doc="when this flag is set (default), statistics will be collected "+
                "on ALIGNED reads only, while unmapped reads will be discarded", required=false)
    protected boolean MAPPED_ONLY = true;
    static private int arrayLen = 2;
   

    public void initialize() 
    {
    	
    }


    public Long[] map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) 
    {
    	Long[] result = new Long[arrayLen];
    	result[0] = result[1] = 0L;
        if ( AlignmentUtils.isReadUnmapped(read) && MAPPED_ONLY) return result;
        
    	result[1] = 1L;
    	if(read.getFirstOfPairFlag() && read.getAlignmentStart() == read.getMateAlignmentStart() && read.getReadNegativeStrandFlag() == read.getMateNegativeStrandFlag())
    		result[0] = 1L;
        
    	return result;
    }

    /**
     * Provide an initial value for reduce computations.
     *
     * @return Initial value of reduce.
     */
    @Override
    public Long[] reduceInit() 
    { 
       	Long[] result = new Long[arrayLen];
    	for(int i = 0; i < result.length; i++)
    		result[i] = 0L;
    	return result;
    }

    /**
     * Reduces a single map with the accumulator provided as the ReduceType.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return accumulator with result of the map taken into account.
     */
    @Override
    public Long[] reduce(Long[] value, Long sum[]) {
        
    	Long[] result = new Long[arrayLen];
    	for(int i = 0; i < result.length; i++)
    	{
    		result[i] = value[i] + sum[i];
    	}
    	return result;
    }
    
    @Override
	public Long[] treeReduce(Long[] lval, Long[] rval)
	{
    	return reduce(lval,rval);
	}

    public void onTraversalDone(Long[] result) 
    {
    	this.logger.info("total mapped reads: " + result[1]);
    	this.logger.info("Inverted Pairs: " + result[0]);
    	this.logger.info("inverted Pair %: " + (2.0 * result[0] / (1.0 * result[1]) * 100.0) + "%");
    	out.printf("mapped reads=%d%n" + "Inverted Read Pairs=%d%n" + "inverted Pair Percentage=%f%n", result[1], result[0], (2.0 * result[0] / (1.0 * result[1]) * 100.0) );
    		
    }
 }