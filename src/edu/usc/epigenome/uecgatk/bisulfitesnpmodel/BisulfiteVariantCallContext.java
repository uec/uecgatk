package edu.usc.epigenome.uecgatk.bisulfitesnpmodel;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;

import edu.usc.epigenome.uecgatk.bisulfitesnpmodel.CytosineTypeStatus;

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

public class BisulfiteVariantCallContext{

	public VariantContext vc = null;
    public ReferenceContext refBase;
    public CytosineTypeStatus cts = null;
    public AlignmentContext rawContext = null;

    // Was the site called confidently, either reference or variant?
    public boolean confidentlyCalled = false;
    public boolean emited = false;

    public BisulfiteVariantCallContext(VariantContext vc, AlignmentContext rawContext, boolean confidentlyCalledP, CytosineTypeStatus cts, boolean emited) {
        this.vc = vc;
        this.rawContext = rawContext;
        this.confidentlyCalled = confidentlyCalledP;
        this.cts = cts;
        this.emited = emited;
    }

    public BisulfiteVariantCallContext(VariantContext vc, AlignmentContext rawContext, ReferenceContext ref, boolean confidentlyCalledP, CytosineTypeStatus cts, boolean emited) {
        this.vc = vc;
        this.rawContext = rawContext;
        this.refBase = ref;
        this.confidentlyCalled = confidentlyCalledP;
        this.cts = cts;
        this.emited = emited;
    }

    // blank variant context => we're a ref site
    public BisulfiteVariantCallContext(boolean confidentlyCalledP, CytosineTypeStatus cts, boolean emited) {
        this.confidentlyCalled = confidentlyCalledP;
        this.cts = cts;
        this.emited = emited;
    }

    public void setRefBase(ReferenceContext ref) {
        this.refBase = ref;
    }

    public boolean isVariant() {
        if(this.vc.hasGenotypes()){
        	//System.err.println(this.vc.getGenotype(0).toString());
        	return (!this.vc.getGenotype(0).isHomRef()) && this.emited;
        }
        return false;
    }
    
    public boolean isHetSnp() {
    	 if(this.vc.hasGenotypes()){
         	return (this.vc.getGenotype(0).isHet()) && this.emited;
         }
         return false;
    }
    
}
