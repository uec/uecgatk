package edu.usc.epigenome.uecgatk.BisSNP;

/*
 * Bis-SNP/BisSNP: It is a genotyping and methylation calling in bisulfite treated 
 * massively parallel sequencing (Bisulfite-seq and NOMe-seq) on Illumina platform
 * Copyright (C) <2011>  <Yaping Liu: lyping1986@gmail.com>

 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

import edu.usc.epigenome.uecgatk.YapingWriter.SortingVCFWriterOwn;

public class SortingTcgaVCFWriter extends SortingVCFWriterOwn {

	protected TcgaVCFWriter tcgaInnerWriter = null;
	public SortingTcgaVCFWriter(TcgaVCFWriter innerWriter,
			int maxCachingStartDistance, boolean takeOwnershipOfInner) {
		super(innerWriter, maxCachingStartDistance, takeOwnershipOfInner);
		tcgaInnerWriter = innerWriter;
		// TODO Auto-generated constructor stub
	}

	public SortingTcgaVCFWriter(TcgaVCFWriter innerWriter,
			int maxCachingStartDistance) {
		super(innerWriter, maxCachingStartDistance);
		tcgaInnerWriter = innerWriter;
		// TODO Auto-generated constructor stub
	}

	public TcgaVCFWriter getInnerWriter(){
		return this.tcgaInnerWriter;
	}
	
	public void writerFlush(){
		this.tcgaInnerWriter.writeFlush();
	}
	
}
