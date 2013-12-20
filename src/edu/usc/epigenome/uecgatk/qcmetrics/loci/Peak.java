package edu.usc.epigenome.uecgatk.qcmetrics.loci;

public interface Peak 
{
	public String getContig();
	public void setContig(String contig);
	public int getStart();
	public void setStart(int start);
	public int getEnd();
	public void setEnd(int end);
	public float getHeight();
	public void setHeight(float height);
	public int getWidth();
	public int getSummit();
	public void setSummit(int summit);
	public String getSample();
	public void setSample(String sample);
	public void setArea(float area);
	public float getArea();
}
