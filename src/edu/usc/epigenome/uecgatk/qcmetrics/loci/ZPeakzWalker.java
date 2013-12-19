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
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.fitting.PolynomialFitter;
import org.apache.commons.math3.ml.clustering.CentroidCluster;
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
    protected int WINSIZE = 30;
    @Argument(fullName="dumpvals", shortName="wig", doc="dumps coverage for each window", required=false)
    protected Boolean dump = false;
    
    @Argument(fullName="numClusters", shortName="clust", doc="number of clusters to use for filtering", required=false)
    protected int NUMBER_OF_CLUSTERS = 3;

    @Argument(fullName="curveWindow", shortName="curve", doc="number of windows to use for polynomial fitting", required=false)
    int CURVE_WINDOW = 3;
    
    protected int current_contig;
    int current_index;
    protected String current_contigName;
    class Stats
    {
    	SummaryStatistics widthStats = new SummaryStatistics();
    	SummaryStatistics heightStats = new SummaryStatistics();
    };
    Stats stats = new Stats();
    HashMap<String,short[]> cov = new HashMap<>();
    HashMap<String,float[]> covDv = new HashMap<>();
    

    
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
    		cov.clear();
    		covDv.clear();
    		for(String s : context.getBasePileup().getSamples())
    		{
    			cov.put(s,new short[100_000_000]);
    			covDv.put(s,new float[100_000_000]);
    		}
    	}
    	
    		
    	else
    		for(String s : context.getBasePileup().getSamples())
    		{
    			if(!cov.containsKey(s))
    			{
    				cov.put(s,new short[100_000_000]);
        			covDv.put(s,new float[100_000_000]);
    			}
    			if(context.getBasePileup().getPileupForSample(s).getBases().length > cov.get(s)[current_index]) 
    				 cov.get(s)[current_index] = (short) context.getBasePileup().getPileupForSample(s).getBases().length;
    		}
    		
    	
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
			for(int i=0; i < (CURVE_WINDOW-1)/2;i++)
				out.printf("%f%n", 0.0f);
    	}
    	
    	HashMap<String,ArrayList<Peak>> samplePeaks = new HashMap<>();
    	//for each sample or bam
    	for(String s : cov.keySet())
    	{
	    	//Phase 1: polynomial fit, then derive
	    	calculateDeriv(cov.get(s),covDv.get(s),contigLength);
	    	
	    	//PHASE 2, detect transitions
	    	ArrayList<Peak> peaks = detectTransitions(cov.get(s),covDv.get(s),contigLength);
	    	
	    	//Phase 3: remove noise using clustering
	    	ArrayList<Peak> passFilter = optimizePeaks(peaks);
	    	
	    	//Phase 4: optimize (inline) peak summit locations
	    	optimizePeakCenters(passFilter,cov.get(s));
	    	samplePeaks.put(s, passFilter);
	    	
	    	//cleanup, set sample names
	    	for(Peak p : passFilter)
	    	{
	    		p.setSample(s);
	    		p.setContig(current_contigName);
	    	}
	    	
	    	
    	}
    	
    	
    	
    	//Find overlaps
    	ArrayList<Peak> mergedPeaks = new ArrayList<>();
    	for(String s : samplePeaks.keySet())
    	{
    		mergedPeaks.addAll(samplePeaks.get(s));
    	}
    	
    	ArrayList<ArrayList<Peak>> mergedBucketedPeaks = new ArrayList<>();
    	Collections.sort(mergedPeaks);
    	for(Peak p : mergedPeaks)
    	{
    		if(mergedBucketedPeaks.size() == 0)
    		{
    			ArrayList<Peak> toAdd = new ArrayList<>();
    			toAdd.add(p);
    			mergedBucketedPeaks.add(toAdd);
    		}
    		else
    		{
    			Boolean matched = false;
    			for(Peak x : mergedBucketedPeaks.get(mergedBucketedPeaks.size() -1))
    				if (p.getSummit() - x.getSummit() < 3)
    						matched = true;
    			if(matched)
    				 mergedBucketedPeaks.get(mergedBucketedPeaks.size() -1).add(p);
    			else
    			{
        			ArrayList<Peak> toAdd = new ArrayList<>();
        			toAdd.add(p);
        			mergedBucketedPeaks.add(toAdd);
        		}
    		}
    	}
    	
    	System.err.println(mergedPeaks.size() + " individual peaks, grouped into " + mergedBucketedPeaks.size());
    	
    	
    	
    	
    	
    	
    	
    	
    	if(samplePeaks.keySet().size() == 1)
    	{
    		String onlyKey = samplePeaks.keySet().toArray(new String[1])[0];
    		System.err.println("writing contig vals to wig for " + current_contigName );
	    	int j=0;
	    	for(int i= (CURVE_WINDOW-1)/2;i<=contigLength;i++)
	    	{
	    		if(i > samplePeaks.get(onlyKey).get(j).getEnd() && j < (samplePeaks.get(onlyKey).size()-1))
	    			j++;
	    		if(i >= samplePeaks.get(onlyKey).get(j).getStart() && i <= samplePeaks.get(onlyKey).get(j).getEnd())
	    		{
	    			if(i == samplePeaks.get(onlyKey).get(j).getSummit())
	    				out.printf("%f%n", samplePeaks.get(onlyKey).get(j).getHeight());
	    			else
	    				out.printf("%f%n", samplePeaks.get(onlyKey).get(j).getHeight() / 2);
	    		}
	    		else
	    			out.printf("%f%n", 0f);
	    	}
    	}
    }
    
    
    
    ArrayList<Peak> optimizePeaks( ArrayList<Peak> peaks)
    {
    	System.err.println("Clustering results");
    	MultiKMeansPlusPlusClusterer<Peak> clusterOptimizer = new MultiKMeansPlusPlusClusterer<Peak>(new KMeansPlusPlusClusterer<Peak>(NUMBER_OF_CLUSTERS),100);
    	
    	
    	List<CentroidCluster<Peak>> clusters = clusterOptimizer.cluster(peaks);
    	Collections.sort(clusters,new Comparator<CentroidCluster<Peak>>(){
			@Override
			public int compare(CentroidCluster<Peak> o1, CentroidCluster<Peak> o2)
			{
				return new Double(o1.getCenter().getPoint()[0]).compareTo(o2.getCenter().getPoint()[0]);
			}});
    	
    	for(CentroidCluster<Peak> cp : clusters)
    		System.err.println("\tCluster Center:" + cp.getCenter().getPoint()[0] + ", size=" + cp.getPoints().size());
    	
    	//the smallest group contains the "noise", trash it.
    	ArrayList<Peak> passFilter = new ArrayList<>();
    	for(int i = 1; i < clusters.size();i++)
    		passFilter.addAll(clusters.get(i).getPoints());
    	
    	//sort the filtered peaks by position
    	Collections.sort(passFilter);

    	System.err.println((peaks.size() - passFilter.size()) + "/" + peaks.size() + " removed = " + passFilter.size() + " total peaks");
    	
    	// merge peaks that are too close
    	int removed = 0;
    	for (int i=passFilter.size()-1; i > 0; i--) 
    	{
    		
    		if(passFilter.get(i).getStart() - passFilter.get(i-1).getEnd() < 3)
        	{
        		//merge	
    			passFilter.get(i-1).setEnd(passFilter.get(i).getEnd());
    			if(passFilter.get(i-1).getHeight() < passFilter.get(i).getHeight())
    				passFilter.get(i-1).setHeight(passFilter.get(i).getHeight());
    			passFilter.remove(i);
        		removed++;
        	}    	    
    	}
    	
    	System.err.println(removed + "/" + (removed + passFilter.size()) + " merged = " + passFilter.size() + " total peaks");
    	return passFilter;
    }
    	
    void optimizePeakCenters(ArrayList<Peak> passFilter, final short[] coverageWindows)
    {
    	System.err.println("Optimizing Peak Summits");
    	//phase 5: optimize peak centers
    	for(Peak peak : passFilter)
    	{
    		//get total area, ignore the smallest 20% of peak in this calc
     		float integral = 0;
    		for(int i=peak.getStart();i<=peak.getEnd();i++)
    			if(1.0f * coverageWindows[i] > peak.getHeight() * 0.2)
    				integral+=coverageWindows[i];
    		float half_integral = 0;
    		
    		//get the midpoint
    		for(int i=peak.getStart();half_integral < 0.5 * integral && i<=peak.getEnd();i++)
    		{
    			if(1.0f * coverageWindows[i] > peak.getHeight() * 0.2)
    				half_integral+=coverageWindows[i];
    			peak.setSummit(i);
    		}    		
    	}
    }
    
    
    ArrayList<Peak> detectTransitions(final short[] coverageWindows, final float[] coverageWindowsDerivitive,int lastIndex)
    {
       	System.err.println("Phase Two: Thread on " + current_contigName);
    	float sum = 0; 
    	float prevSum = 0;
    	ArrayList<Peak> peaks = new ArrayList<Peak>();
    	Peak p = null;
    	for(int i= (CURVE_WINDOW-1)/2;i<=lastIndex;i++)
    	{	
    		sum += coverageWindowsDerivitive[i];
    		
    		//starting peak
    		if(p==null && prevSum <= 0 && sum > 0)
    		{
    			p=new Peak();
    			p.setStart(i);
    		}
    		
    		//max height;
    		if(p!=null && coverageWindows[i] > p.getHeight())
    		{
    			p.setHeight(coverageWindows[i]);
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
    	
    	System.err.println("Width Avg: " + stats.widthStats.getMean() + " +- "+ stats.widthStats.getStandardDeviation());
    	System.err.println("Height Avg: " + stats.heightStats.getMean() + " +- "+ stats.heightStats.getStandardDeviation());
    	
    	
    	return peaks;
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
    	double[] guess = {10.0,10.0,10.0};
    	PolynomialFitter fitter = new PolynomialFitter(new LevenbergMarquardtOptimizer());
    	for(int i= start + (CURVE_WINDOW-1)/2;i<=end;i++)
    	{
    		fitter.clearObservations();
    		for(int j = 0;j<CURVE_WINDOW;j++)
    			fitter.addObservedPoint(j, coverageWindows[ (i - ((CURVE_WINDOW-1)/2)) + j ]);
    		double[] results = fitter.fit(guess);
    		final PolynomialFunction fitted = new PolynomialFunction(results);
    		coverageWindowsDerivitive[i] = (float) fitted.derivative().value((CURVE_WINDOW-1)/2);
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