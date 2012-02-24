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

import org.apache.commons.math.stat.descriptive.SummaryStatistics;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.Argument;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import java.util.*;
import java.io.*;

/**
 * Zack Ramjan
 * USC Epigenome Center 
 * 03/23/2011
 */

/**
 * Walks over the input data set, counts various stats regarding alignment length of the reads
 */
@Requires({DataSource.READS})
public class ReadLengthWalker extends ReadWalker<ReadLengthWalker.ReadLenInfo,Integer> {
    @Output
    protected PrintStream out;

    @Argument(fullName="mappedOnly", shortName="mo", doc="when this flag is set (default), statistics will be collected on mapped reads only, while unmapped reads will be discarded", required=false)
    protected boolean MAPPED_ONLY = true;

    @Argument(fullName="skip", shortName="skip", doc="When provided, only every skip reads are analyzed", required=false)
    protected int SKIP = 1;

    HashMap<String,SummaryStatistics> stats;
    
    
    public void initialize() 
    {
    	stats=new HashMap<String,SummaryStatistics>();
    	for (String s : Arrays.asList("ReadLength", "NumClippingEvents", "NumClippedBases", "PercentClipped"))
    		stats.put(s, new SummaryStatistics());
    	
    }

    public class ReadLenInfo 
    {
      
        int readLength, nClippingEvents, nClippedBases;
    }

    public ReadLenInfo map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) 
    {
        if ( AlignmentUtils.isReadUnmapped(read) && MAPPED_ONLY)
            return null;

        ReadLenInfo info = new ReadLenInfo();
        for ( CigarElement elt : read.getCigar().getCigarElements() ) {
            if ( elt.getOperator() != CigarOperator.N )

            switch ( elt.getOperator()) {
                case H : // ignore hard clips
                case S : // soft clip
                    info.nClippingEvents++;
                    info.nClippedBases += elt.getLength();
                    // note the fall through here
                case M :
                case D : // deletion w.r.t. the reference
                case P : // ignore pads
                case I : // insertion w.r.t. the reference
                    info.readLength += elt.getLength(); // Unless we have a reference skip, the read gets longer
                    break;
                case N : // reference skip (looks and gets processed just like a "deletion", just different logical meaning)
                    break;
                default : throw new IllegalStateException("Case statement didn't deal with cigar op: " + elt.getOperator());
            }
        }

        return info;  //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * Provide an initial value for reduce computations.
     *
     * @return Initial value of reduce.
     */
    public Integer reduceInit() {
        return 0;
    }

    /**
     * Reduces a single map with the accumulator provided as the ReduceType.
     *
     * @param info  result of the map.
     * @param sum   accumulator for the reduce.
     * @return accumulator with result of the map taken into account.
     */
    public Integer reduce(ReadLenInfo info, Integer sum) {
        if ( info != null ) {
            if ( sum % SKIP == 0 ) 
            {
           
            	stats.get("ReadLength").addValue(info.readLength);
            	stats.get("NumClippingEvents").addValue(info.nClippingEvents);
            	stats.get("NumClippedBases").addValue(info.nClippedBases);
            	stats.get("PercentClipped").addValue( 100.0 * MathUtils.ratio(info.nClippedBases, info.readLength));           
            }
            return sum + 1;  //To change body of implemented methods use File | Settings | File Templates.
        } else {
            return sum;
        }
    }

    public void onTraversalDone(Integer result) 
    {
    	out.print("NumReads=" + stats.get("ReadLength").getN() + "	"); 
    	for (String s : Arrays.asList("ReadLength", "NumClippingEvents", "NumClippedBases", "PercentClipped"))
    		out.print(s + ".avg=" + stats.get(s).getMean() + "	" + s + ".stdv=" + stats.get(s).getStandardDeviation() + "	");
    	out.println();
    }
    
}