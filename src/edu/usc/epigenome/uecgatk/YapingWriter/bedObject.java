package edu.usc.epigenome.uecgatk.YapingWriter;

import java.util.ArrayList;
import java.util.List;

import org.broad.tribble.annotation.Strand;

public class bedObject implements genomeObject {

	private String chr;
	private int start;
	private int end;
	private Strand strand = Strand.NONE;
	private List<Object> values;
	
	public bedObject(String chr, int start, int end, List<Object> values) {
		// TODO Auto-generated constructor stub
		this.chr = chr;
		this.start = start;
		this.end = end;
		this.values = values;
		
	}
	
	public bedObject(String chr, int start, int end, Strand strand, List<Object> values) {
		// TODO Auto-generated constructor stub
		this.chr = chr;
		this.start = start;
		this.end = end;
		this.values = values;
		this.strand = strand;
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
	
	public Strand getStrand() {
		// TODO Auto-generated method stub
		return strand;
	}
	
	public String getStrandAsString() {
		// TODO Auto-generated method stub
		return strand.toString();
	}
	
	public List<Object> getValueObject() {
		// TODO Auto-generated method stub
		return values;
	}
	
	

}
