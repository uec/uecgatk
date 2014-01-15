package edu.usc.epigenome.uecgatk.qcmetrics.loci;

import java.util.HashMap;
import java.util.List;

public class PeakDiffCalcRPKM implements PeakDiffCalc
{
	String[] sample_names;
	HashMap<String,Long> totalReads;
	HashMap<String,short[]> cov;
	
	PeakDiffCalcRPKM(HashMap<String,Long> totalReads,HashMap<String,short[]> cov)
	{
		this.totalReads = totalReads;
		this.cov = cov;
		sample_names = totalReads.keySet().toArray(new String[ totalReads.keySet().size()]);
	}
	
	@Override
	public float calcDiff(List<SinglePeak> peaks)
	{
		HashMap<String,Float> sampleTotal = new HashMap<>();
		int minCoord = 1_000_000_000;
		int maxCoord = -1;
		
		for(String s : sample_names)
             sampleTotal.put(s, 0f);
		 
		
		for(Peak p : peaks)
		{
			minCoord = p.getStart() <= minCoord ? p.getStart() : minCoord;
			maxCoord = p.getEnd() >= maxCoord ? p.getEnd() : maxCoord;
		}
		
		for(String s : sample_names)
		{
			float area = 0f;
			int width = maxCoord - minCoord; 
			for(int i = minCoord; i <= maxCoord; i++)
				area += cov.get(s)[i];
			float rpkm = 1_000_000_000 * area / (width * totalReads.get(s));
			sampleTotal.put(s, sampleTotal.get(s) + rpkm);
		}
		
		float max = -1;
		float min = 10000000000f;
		
		for(String s : sample_names)
		{
			if(sampleTotal.get(s) >  max)
				max = sampleTotal.get(s);
			if(sampleTotal.get(s) <  min)
				min  = sampleTotal.get(s);
		}		
			
		//if there is only 1 sample, return the rpkm * -1
		if(sample_names.length == 1)
			return max;
		
		for(String s : sample_names)
			if(sampleTotal.get(s) == 0f)
				return max * -1;
		
		return (max/min) - 1;
	}
}
