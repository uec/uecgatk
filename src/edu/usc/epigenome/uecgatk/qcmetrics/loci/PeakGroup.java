package edu.usc.epigenome.uecgatk.qcmetrics.loci;

import java.util.ArrayList;

public class PeakGroup extends ArrayList<SinglePeak> implements Peak,Comparable<Peak>
{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private String[] sample_names;
	PeakDiffCalc diffcalc;
	
	
	public PeakGroup()
	{
		super();
	}

	public PeakGroup(String[] samples)
	{
		super();
		this.sample_names = samples;
	}
	
	
	public PeakGroup(PeakDiffCalc diffcalc)
	{
		super();
		this.diffcalc = diffcalc;
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
		return diffcalc.calcDiff(this);
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
	
	@Override
	public int compareTo(Peak o)
	{
		return new Integer(this.getStart()).compareTo(o.getStart());
	}

	@Override
	public void setArea(float area) {}

	@Override
	public float getArea()
	{
		// TODO Auto-generated method stub
		return 0;
	}
}