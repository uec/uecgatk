package edu.usc.epigenome.uecgatk.YapingWriter;

import java.util.ArrayList;

public class bedObject implements genomeObject {

	private String chr;
	private int start;
	private int end;
	private ArrayList<Object> values;
	
	public bedObject(String chr, int start, int end, ArrayList<Object> values) {
		// TODO Auto-generated constructor stub
		this.chr = chr;
		this.start = start;
		this.end = end;
		this.values = values;
		
	}

	@Override
	public int getStart() {
		// TODO Auto-generated method stub
		return start;
	}

	@Override
	public String getChr() {
		// TODO Auto-generated method stub
		return chr;
	}
	
	public int getEnd() {
		// TODO Auto-generated method stub
		return end;
	}
	
	public ArrayList<Object> getValueObject() {
		// TODO Auto-generated method stub
		return values;
	}
	
	

}
