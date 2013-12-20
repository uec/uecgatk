package edu.usc.epigenome.uecgatk.qcmetrics.loci;

import org.apache.commons.math3.ml.clustering.Clusterable;

public class SinglePeak implements Peak,Clusterable,Comparable<Peak>
{

	String contig;
	String sample;
	int start;
	int end;
	float height;
	int summit;
	float area;
	
	public SinglePeak(){};
	public SinglePeak(String contig,int start,int end, int peak, float height, int summit)
	{
		this.contig = contig;
		this.start = start;
		this.end = end;
		this.height = height;
		this.summit = summit;
	};
	@Override
	public String getContig()
	{
		return contig;
	}
	@Override
	public void setContig(String contig)
	{
		this.contig = contig;
	}
	@Override
	public int getStart()
	{
		return start;
	}
	@Override
	public void setStart(int start)
	{
		this.start = start;
	}
	@Override
	public int getEnd()
	{
		return end;
	}
	@Override
	public void setEnd(int end)
	{
		this.end = end;
	}
	@Override
	public float getHeight()
	{
		return height;
	}
	@Override
	public void setHeight(float height)
	{
		this.height = height;
	}
	@Override
	public int getWidth()
	{
		return end-start;
	}
	@Override
	public int getSummit()
	{
		return summit;
	}
	@Override
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
	@Override
	public void setArea(float area)
	{
		this.area = area;
		
	}
	@Override
	public float getArea()
	{
		return area;
	}
}
