package edu.usc.epigenome.uecgatk.YapingWriter;

import java.util.zip.CRC32;

public class cpgReads implements genomeObject {

	private String chr;
	private int genomeLoc;
	private char methyStatus;
	private byte baseQ;
	private char strand;
	private String readID;
	private CRC32 encrypt;
	
	
	public cpgReads(String chr, int genomeLoc, char methyStatus, byte baseQ, char strand, String readID){
		this.chr = chr;
		this.genomeLoc = genomeLoc;
		this.methyStatus = methyStatus;
		this.baseQ = baseQ;
		this.strand = strand;
		this.readID = readID;
		this.encrypt = new CRC32();
		encrypt.update(readID.getBytes());
	}
	
	
	public char getMethyStatus(){
		return this.methyStatus;
	}
	
	public byte getbaseQ(){
		return this.baseQ;
	}
	
	public char getstrand(){
		return this.strand;
	}
	
	public String getReadID(){
		return this.readID;
	}
	
	public long getEncryptID(){	
		return this.encrypt.getValue();
	}
	
	@Override
	public int getStart() {
		// TODO Auto-generated method stub
		return this.genomeLoc;
	}

	@Override
	public String getChr() {
		// TODO Auto-generated method stub
		return this.chr;
	}

}
