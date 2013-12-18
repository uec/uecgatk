package edu.usc.epigenome.uecgatk.qcmetrics.loci;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.fitting.PolynomialFitter;
import org.apache.commons.math3.ml.clustering.CentroidCluster;
import org.apache.commons.math3.ml.clustering.Cluster;
import org.apache.commons.math3.ml.clustering.Clusterable;
import org.apache.commons.math3.ml.clustering.KMeansPlusPlusClusterer;
import org.apache.commons.math3.ml.clustering.MultiKMeansPlusPlusClusterer;
import org.apache.commons.math3.optim.nonlinear.vector.jacobian.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.stat.descriptive.*;

import edu.usc.epigenome.uecgatk.filters.NonUniqueFilter;



/**
 * Zack Ramjan
 * USC Epigenome Center 
 * 06/10/2011
 * 08/08/2012
 */

/**
 * Bin Depths walker. calculate coverage in windows across genome.
 * report stats upon these windows
 */
@By(DataSource.REFERENCE)
@ReadFilters( {NonUniqueFilter.class} ) // Filter out all reads with zero mapping quality
public class ZPeakzWalker extends LocusWalker<Boolean,Boolean>  
{
    @Output
    PrintStream out;
    
    @Argument(fullName="winSize", shortName="winsize", doc="window width", required=false)
    protected int WINSIZE = 50;
    @Argument(fullName="dumpvals", shortName="wig", doc="dumps coverage for each window", required=false)
    protected Boolean dump = false;
    

    protected int current_contig;
    int current_index;
    protected String current_contigName;
    class Stats
    {
    	SummaryStatistics widthStats = new SummaryStatistics();
    	SummaryStatistics heightStats = new SummaryStatistics();
    };
    Stats stats = new Stats();
    protected double PROBABILITY = 1.0;
    short[] coverage = new short[100_000_000];
    float[] coverageDeriv = new float[100_000_000];
    int curveWindow = 3;
    int threads = 4;

    public void initialize() 
    {
    	  	 current_contig = -1;
    	 current_contigName = "";
    	 current_index = 0;
    	 out.println("track type=wiggle_0 name=\"peaks\" description=\"winsize=" + WINSIZE + "\"");
    }
    
    /**
     * We pretty much bypass the whole map/reduce stuff since we use a global to keep track of counts
     * we scan the genome linearly, saving counts and resetting when ever we cross a threshold.
     * since we skipped the M/R, this will not be parallel/treereducable.
     * 
     * basically, we have a single for loop
     * 
     * @param tracker The accessor for reference metadata.
     * @param ref The reference base that lines up with this locus.
     * @param context Information about reads aligning to this locus.
     * @return return average dup count of all trials.
     */
    @Override
    public Boolean map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) 
    {
    	int contig = ref.getLocus().getContigIndex();
    	String contigName = ref.getLocus().getContig();
    	long pos = context.getPosition();
    	int depth = context.getBasePileup().getBases().length;
    	
    	if(pos % WINSIZE == 0)
    		current_index++;
    	
    	
    	if(current_contig != contig)
    	{
    		if(current_contig >= 0)
    			processContig(pos,current_index);
    		
    		//reset counts for new contig
    		current_contig=contig;
    		current_contigName = contigName;
    		current_index = 0;	    		
    	}
    	
    		
    	else
    		if(depth > coverage[current_index]) 
    			coverage[current_index] = (short) depth;
    	
    	return null;
    }

    
    
    /**
     * Provides an initial value for the reduce function. and array of 0's
     * @return always return true.
     */
    @Override
    public Boolean reduceInit() 
    { 
    	return true;
    }

    /**
     * Combines the result of the latest map with the accumulator.  In inductive terms,
     * this represents the step loci[x + 1] = loci[x] + 1
     * @param just a bogus param for the override. 
     * @param just a bogus param for the override. 
     * @return always truen true.
     */
    @Override
    public Boolean reduce(Boolean value, Boolean val2) {
        
    	return true;
    }

    /**
     * Retrieves the final result of the traversal.
     * @param just a bogus param for the override. 
     */
    @Override
    public void onTraversalDone(Boolean t) 
    {
    	processContig(0,current_index);
    }
    
    protected void processContig(long pos,int contigLength)
    {
    	if(dump)
    	{
			out.printf("fixedStep  chrom=%s start=%d step=%d span=%d%n", current_contigName, pos, WINSIZE,WINSIZE);
			for(int i=0; i < (curveWindow-1)/2;i++)
				out.printf("%f%n", 0.0f);
    	}
    	
    	//Phase One: polynomial fit, then derive
    	calculateDeriv(coverage,coverageDeriv,contigLength);
    	
    	//PHASE two, detect transitions
    	System.err.println("Phase Two: Thread on " + current_contigName);
    	float sum = 0; 
    	float prevSum = 0;
    	ArrayList<Peak> peaks = new ArrayList<Peak>();
    	Peak p = null;
    	for(int i= (curveWindow-1)/2;i<=contigLength;i++)
    	{	
    		sum += coverageDeriv[i];
    		
    		//starting peak
    		if(p==null && prevSum <= 0 && sum > 0)
    		{
    			p=new Peak();
    			p.setStart(i);
    		}
    		
    		//max height;
    		if(p!=null && coverage[i] > p.height)
    		{
    			p.setHeight(coverage[i]);
    			p.setSummit(i);
    		}
    		
    		//ending peak
    		if(p != null && sum <= 0 && prevSum > sum )
    		{
    			p.setEnd(i);
    			peaks.add(p);
    			stats.heightStats.addValue(p.getHeight());
    			stats.widthStats.addValue(p.getEnd() - p.getStart());
    			p=null;
    			sum=0;
    		}
    		prevSum = sum;    		
    	}
    	
    	System.err.println("Width: " + stats.widthStats.getMean() + " +- "+ stats.widthStats.getStandardDeviation());
    	System.err.println("Height" + stats.heightStats.getMean() + " +- "+ stats.heightStats.getStandardDeviation());
    	
    	//Phase 3
    	//remove noise
    	
    	//using clustering
    	System.err.println("Clustering results");
    	KMeansPlusPlusClusterer<Peak> clusterer = new KMeansPlusPlusClusterer<>(3);
    	MultiKMeansPlusPlusClusterer<Peak> clusterOptimizer = new MultiKMeansPlusPlusClusterer<Peak>(clusterer,20);
    	
    	
    	List<CentroidCluster<Peak>> clusters = clusterOptimizer.cluster(peaks);
    	Collections.sort(clusters,new Comparator<CentroidCluster<Peak>>(){
			@Override
			public int compare(CentroidCluster<Peak> o1, CentroidCluster<Peak> o2)
			{
				return new Double(o1.getCenter().getPoint()[0]).compareTo(o2.getCenter().getPoint()[0]);
			}});
    	
    	for(CentroidCluster<Peak> cp : clusters)
    		System.err.println("Cluster Center:" + cp.getCenter().getPoint()[0]);
    	
    	ArrayList<Peak> passFilter = new ArrayList<>();
    	//passFilter.addAll(clusters.get(0).getPoints());
    	passFilter.addAll(clusters.get(1).getPoints());
    	passFilter.addAll(clusters.get(2).getPoints());
    	
    	//sort the filtered collection
    	Collections.sort(passFilter,new Comparator<Peak>(){
			@Override
			public int compare(Peak o1, Peak o2)
			{
				return new Integer(o1.getStart()).compareTo(o2.getStart());
			}});    	
//    	int removed = 0;
//    	for (int i=peaks.size()-1; i> -1; i--) 
//    	{
//    		
//    		if(
//        	peaks.get(i).getWidth() < stats.widthStats.getMean() + stats.widthStats.getStandardDeviation() &&
//        	peaks.get(i).getHeight() < stats.heightStats.getMean() + stats.heightStats.getStandardDeviation())
//        	{
//        			peaks.remove(i);
//        			removed++;
//        	}    	    
//    	}
//    	System.err.println(removed + "/" + stats.heightStats.getN() + " removed = " + peaks.size());
    	System.err.println((peaks.size() - passFilter.size()) + "/" + peaks.size() + " removed = " + passFilter.size());
    	
    	
    	peaks=passFilter;
    	
    	
    	//phase 4 optimize peak centers
    	for(Peak peak : peaks)
    	{
    		//using polynomial fitting
//    		double[] guess = {10.0,10.0,10.0};
//        	PolynomialFitter fitter = new PolynomialFitter(new LevenbergMarquardtOptimizer());
//        	for(int i=peak.getStart();i<=peak.getEnd();i++)
//        		fitter.addObservedPoint(i,coverage[i]);
//        	double[] results = fitter.fit(guess);
//        	final PolynomialFunction fitted = new PolynomialFunction(results);
    		float integral = 0;
    		for(int i=peak.getStart();i<=peak.getEnd();i++)
    			integral+=coverage[i];
    		float half_integral = 0;
    		for(int i=peak.getStart();half_integral < 0.5 * integral && i<=peak.getEnd();i++)
    		{
    			half_integral+=coverage[i];
    			peak.setSummit(i);
    		}
    		
    		
    	}
    	
    	
    	System.err.println("writing contig vals to wig for " + current_contigName );
    	int j=0;
    	for(int i= (curveWindow-1)/2;i<=contigLength;i++)
    	{
    		if(i > peaks.get(j).getEnd() && j < (peaks.size()-1))
    			j++;
    		if(i >= peaks.get(j).getStart() && i <= peaks.get(j).getEnd())
    		{
    			if(i == peaks.get(j).getSummit())
    				out.printf("%f%n", peaks.get(j).getHeight());
    			else
    				out.printf("%f%n", peaks.get(j).getHeight() / 2);
    		}
    		else
    			out.printf("%f%n", 0f);
    	}
   
    }
   
    float[] calculateDeriv(final short[] coverageWindows, final float[] coverageWindowsDerivitive, int lastIndex)    
    {
    	class PhaseOne implements Runnable 
    	{
    		int start;
    		int end;
    		PhaseOne(int s, int e) {start=s;end=e;}

			@Override
			public void run() { computePhaseOne(start,end,coverageWindows,coverageWindowsDerivitive);}    	
    	}

    	
    	//phase one
    	//create polynomial and derivative array
    	List<Runnable> jobs = new ArrayList<Runnable>();
    	for (int i = 0; i < threads; i++)
    		jobs.add(new PhaseOne(i*(lastIndex/threads),(i+1) == threads ? lastIndex : ((i+1)*(lastIndex/threads))));    		
    	runParallel(jobs);	
    	return coverageWindowsDerivitive;
    }
    
    
    void computePhaseOne(int start,int end,short[] coverageWindows, float[] coverageWindowsDerivitive)
    {
    	System.err.println("Phase One: Thread on " + current_contigName + ": " + start + ":" + end);
    	double[] guess = {10.0,10.0,10.0};
    	PolynomialFitter fitter = new PolynomialFitter(new LevenbergMarquardtOptimizer());
    	for(int i= start + (curveWindow-1)/2;i<=end;i++)
    	{
    		fitter.clearObservations();
    		for(int j = 0;j<curveWindow;j++)
    			fitter.addObservedPoint(j, coverage[ (i - ((curveWindow-1)/2)) + j ]);
    		double[] results = fitter.fit(guess);
    		final PolynomialFunction fitted = new PolynomialFunction(results);
    		coverageDeriv[i] = (float) fitted.derivative().value((curveWindow-1)/2);
    		//System.err.println(i +  "\t" + coverageDeriv[i]);
    	}
    	System.err.println("Phase One: Thread on " + current_contigName + ": " + start + ":" + end + " complete");
    }
    
    private void runParallel(List<Runnable> jobs)
    {
    	ExecutorService  pool = Executors.newFixedThreadPool(threads);
    	for(Runnable r : jobs)
    		pool.execute(r);
    	pool.shutdown();
    	try
		{	
			pool.awaitTermination(5,TimeUnit.DAYS);
		} catch (InterruptedException e)
		{
			System.err.println("RAN OUT OF TIME");
			e.printStackTrace();
		}
    }
 }