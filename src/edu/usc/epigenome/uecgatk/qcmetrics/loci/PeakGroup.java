package edu.usc.epigenome.uecgatk.qcmetrics.loci;

import java.util.ArrayList;
import java.util.HashMap;

public class PeakGroup extends ArrayList<SinglePeak> implements Peak,Comparable<Peak>
{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private String[] sample_names;
	
	
	public PeakGroup()
	{
		super();
	}

	public PeakGroup(String[] samples)
	{
		super();
		this.sample_names = samples;
	}
	
	@Override
	public String getContig()
	{
		return this.get(0).getContig();
	}

	@Override
	public void setContig(String contig){ }

	@Override
	public int getStart()
	{
		int start = this.get(0).getStart();
		for(Peak p : this)
		{
			if(p.getStart() < start)
				start = p.getStart();
		}
		return start;
	}

	@Override
	public void setStart(int start){ }

	@Override
	public int getEnd()
	{
		int start = this.get(0).getEnd();
		for(Peak p : this)
		{
			if(p.getEnd() > start)
				start = p.getEnd();
		}
		return start;
	}

	@Override
	public void setEnd(int end){ }

	@Override
	public float getHeight()
	{
		HashMap<String,Float> sampleTotal = new HashMap<>();
		for(String s : sample_names)
			sampleTotal.put(s, 0f);
		
		for(Peak p : this)
			sampleTotal.put(p.getSample(), sampleTotal.get(p.getSample()) + (p.getHeight() * this.getSampleIndex(p.getSample())));
		
		float total = 0f;
		for(float f : sampleTotal.values())
			total += f;
		return total;
	}

	@Override
	public void setHeight(float height){ }

	@Override
	public int getWidth()
	{
		return this.getEnd() - this.getStart();
	}

	@Override
	public int getSummit()
	{
		float start = this.get(0).getHeight();
		int pos = this.get(0).getSummit();
		for(Peak p : this)
		{
			if(p.getHeight() > start)
			{
				start = p.getHeight();
				pos = p.getSummit();
			}
		}
		return pos;
	}

	@Override
	public void setSummit(int summit){ }

	@Override
	public String getSample()
	{
		String sample = "";
		for(String s : this.sample_names)	
			sample += s;
		return sample;
	}

	@Override
	public void setSample(String sample) 	{ }
	
	private int getSampleIndex(String mySampleName)
	{
		for (int i=0;i< sample_names.length; i++)
		{
			if(sample_names[i].equals(mySampleName))
				return i % 2 == 0 ? 1 : -1;
		}
		return 0;
	}

	@Override
	public int compareTo(Peak o)
	{
		return new Integer(this.getStart()).compareTo(o.getStart());
	}
}