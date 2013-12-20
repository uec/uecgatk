package edu.usc.epigenome.uecgatk.qcmetrics.loci;

import java.util.HashMap;
import java.util.List;

public class PeakDiffCalcRPKM implements PeakDiffCalc
{
	String[] sample_names;
	HashMap<String,Long> totalReads;
	PeakDiffCalcRPKM(HashMap<String,Long> totalReads)
	{
		this.totalReads = totalReads;
		sample_names = totalReads.keySet().toArray(new String[ totalReads.keySet().size()]);
	}
	
	@Override
	public float calcDiff(List<SinglePeak> peaks)
	{
		HashMap<String,Float> sampleTotal = new HashMap<>();
		float max = 0;
		float min = 10000000000f;
		for(String s : sample_names)
			sampleTotal.put(s, 0f);
		
		for(Peak p : peaks)
		{
			float rpkm = 1_000_000_000 * p.getArea() / (p.getWidth() * totalReads.get(p.getSample()));
			sampleTotal.put(p.getSample(), sampleTotal.get(p.getSample()) + rpkm);
		}
		
		for(String s : sample_names)
		{
			if(sampleTotal.get(s) >  max)
				max = sampleTotal.get(s);
			if(sampleTotal.get(s) <  min)
				min  = sampleTotal.get(s);
		}		
			
		//if there is only 1 sample, return the rpkm * -1
		for(String s : sample_names)
			if(sampleTotal.get(s) == 0f)
				return max * -1;
		
		return max/min;
	}
}
