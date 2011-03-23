/*
 * Copyright (c) 2010 The Broad Institute
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
import org.broadinstitute.sting.utils.collections.PrimitivePair;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMReadGroupRecord;

import java.util.*;
import java.io.*;

/**
 * Zack Ramjan
 * USC Epigenome Center
 
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
    protected boolean MAPPED_ONLY = true;

    private HashMap<String,Long> bqmap = null;
    
    public void initialize() 
    {
    	bqmap = new HashMap<String,Long>();
    }


    public HashMap<String,Long> map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) 
    {

        if ( AlignmentUtils.isReadUnmapped(read) && MAPPED_ONLY) return null;
        
        byte [] quals = read.getBaseQualities();
        byte [] bases =  read.getReadBases();
        
        
        HashMap<String,Long> countRead = new HashMap<String,Long>();
        //a,+,8,27,156335
        for (short i = 0; i < bases.length; i++)
        {
        	countRead.put(bases[i] + ",+," + (i+1) + "," + quals[i] , 1L);
        }
        
        
        return countRead;
    }

    /**
     * Provide an initial value for reduce computations.
     *
     * @return Initial value of reduce.
     */
    public HashMap<String,Long> reduceInit() {
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
    	for(String s : result.keySet())
    	{
    		out.println(s + "," + result.get(s));
    	}
    }

 

   


   


    static class CycleStats {
        private long readCount = 0;
        private double[] cycleQualsAv = null;
        private double[] cycleQualsSd = null;
        private int minL = 1000000000; // read min. length
        private int maxL = 0; // read max. length

        public CycleStats(int N) {
            readCount = 0;
            cycleQualsAv = new double[N];
            cycleQualsSd = new double[N];
        }

        public void add(byte[] quals) {
            if ( quals.length > cycleQualsAv.length )
                throw new UserException("A read of length "+quals.length+" encountered, which exceeds specified maximum read length");
            if ( quals.length > maxL ) maxL = quals.length;
            if ( quals.length < minL ) minL = quals.length;
            readCount++;
            for ( int i = 0 ; i < quals.length ; i++ ) {
                // NOTE: in the update equaltions below, there is no need to check if readCount == 1 (i.e.
                // we are initializing with the very first record) or not. Indeed, the arrays are initialized with
                // 0; when the very first value arrives, readCount is 1 and cycleQuals[i] gets set to quals[i] (correct!);
                // this will also make the second term in the update equation for Sd (quals[i]-cycleQualsAv[i]) equal
                // to 0, so Sd will be initially set to 0.
                double oldAvg = cycleQualsAv[i]; // save old mean, will need it for calculation of the variance
                cycleQualsAv[i] += ( quals[i] - cycleQualsAv[i] ) / readCount; // update mean
                cycleQualsSd[i] += ( quals[i] - oldAvg ) * ( quals[i] - cycleQualsAv[i] );
            }
        }

        public long getReadCount() { return readCount; }
        public int getMaxReadLength() { return maxL; }
        public int getMinReadLength() { return minL; }
//        long [] getCycleQualSums() { return cycleQuals; }
//        long getCycleQualSum(int i) { return cycleQuals[i]; }
        double getCycleQualAverage(int i) { return cycleQualsAv[i]; }
        double getCycleQualStdDev(int i) { return Math.sqrt( cycleQualsSd[i]/(readCount-1) ); }
    }
}
