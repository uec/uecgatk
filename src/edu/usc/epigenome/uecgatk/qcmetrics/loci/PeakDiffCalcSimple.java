package edu.usc.epigenome.uecgatk.qcmetrics.loci;

import java.util.HashMap;
import java.util.List;

public class PeakDiffCalcSimple implements PeakDiffCalc
{
	String[] sample_names;
	PeakDiffCalcSimple(String[] samples)
	{
		sample_names = samples;
	}
	
	@Override
	public float calcDiff(List<SinglePeak> peaks)
	{
		HashMap<String,Float> sampleTotal = new HashMap<>();
		for(String s : sample_names)
			sampleTotal.put(s, 0f);
		
		for(Peak p : peaks)
			sampleTotal.put(p.getSample(), sampleTotal.get(p.getSample()) + (p.getHeight() * this.getSampleIndex(p.getSample())));
		
		float total = 0f;
		for(float f : sampleTotal.values())
			total += f;
		return total;
	
	}
	
	private int getSampleIndex(String mySampleName)
	{
		for (int i=0;i< sample_names.length; i++)
		{
			if(sample_names[i].equals(mySampleName))
				return i % 2 == 0 ? 1 : -1;
		}
		return 0;
	}
}
