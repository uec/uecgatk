package edu.usc.epigenome.uecgatk.qcmetrics.loci;

import org.apache.commons.math3.ml.clustering.Clusterable;

public class Peak implements Clusterable,Comparable<Peak>
{

	String contig;
	String sample;
	int start;
	int end;
	float height;
	int summit;
	
	public Peak(){};
	public Peak(String contig,int start,int end, int peak, float height, int summit)
	{
		this.contig = contig;
		this.start = start;
		this.end = end;
		this.height = height;
		this.summit = summit;
	};
	public String getContig()
	{
		return contig;
	}
	public void setContig(String contig)
	{
		this.contig = contig;
	}
	public int getStart()
	{
		return start;
	}
	public void setStart(int start)
	{
		this.start = start;
	}
	public int getEnd()
	{
		return end;
	}
	public void setEnd(int end)
	{
		this.end = end;
	}
	public float getHeight()
	{
		return height;
	}
	public void setHeight(float height)
	{
		this.height = height;
	}
	public int getWidth()
	{
		return end-start;
	}
	public int getSummit()
	{
		return summit;
	}
	public void setSummit(int summit)
	{
		this.summit = summit;
	}
	
	
	@Override
	public double[] getPoint()
	{
		return new double[] {this.getHeight()};
	}
	
	@Override
	public int compareTo(Peak o)
	{
		// TODO Auto-generated method stub
		return new Integer(this.getStart()).compareTo(o.getStart());
		
	}
	public String getSample()
	{
		return sample;
	}
	public void setSample(String sample)
	{
		this.sample = sample;
	}
}
