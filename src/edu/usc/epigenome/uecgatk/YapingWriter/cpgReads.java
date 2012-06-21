package edu.usc.epigenome.uecgatk.YapingWriter;

import java.io.UnsupportedEncodingException;
import java.util.zip.CRC32;

import jonelo.jacksum.algorithm.Crc64;
import jonelo.jacksum.util.Service;

public class cpgReads implements genomeObject {

	private String chr;
	private int genomeLoc;
	private char methyStatus;
	private byte baseQ;
	private char strand;
	private String readID;
	//private CRC32 encrypt;
	private Crc64 encrypt64;
	
	
	public cpgReads(String chr, int genomeLoc, char methyStatus, byte baseQ, char strand, String readID){
		this.chr = chr;
		this.genomeLoc = genomeLoc;
		this.methyStatus = methyStatus;
		this.baseQ = baseQ;
		this.strand = strand;
		this.readID = readID;
		//this.encrypt = new CRC32();
		//encrypt.update(readID.getBytes());
		this.encrypt64 = new Crc64();
		encrypt64.update(readID.getBytes());
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
	
	//public long getEncryptID(){	
	//	return this.encrypt.getValue();
	//}
	
	//public String getEncryptIDAscii(){	
	//	return Long.toHexString(this.encrypt.getValue());
	//}
	
	public long getEncryptID64(){	
		return this.encrypt64.getValue();
	}
	
	public String getEncryptID64Ascii(){
		
		return Service.format(this.encrypt64.getByteArray());
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
