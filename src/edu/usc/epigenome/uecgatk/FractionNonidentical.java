package edu.usc.epigenome.uecgatk;

/**
 * @author benb
 *
 * Representation of a fraction where 1/2 does not equal 2/4
 */
public class FractionNonidentical extends Number {
    protected int numerator;
    protected int denominator;

    public FractionNonidentical(int numerator, int denominator) {
        if(denominator == 0) {
            throw new IllegalArgumentException("denominator is zero");
        }
        if(denominator < 0) {
            numerator *= -1;
            denominator *= -1;
        }
        this.numerator = numerator;
        this.denominator = denominator;
    }

    public FractionNonidentical(int numerator) {
        this.numerator = numerator;
        this.denominator = 1;
    }

    public FractionNonidentical() {
        this.numerator = 0;
        this.denominator = 0;
    }
    
    public int getNumerator() {
        return this.numerator;
    }

    public int getDenominator() {
        return this.denominator;
    }
    
    public void incNumerator()
    {
    	this.numerator++;
    }
    
    public void incDenominator()
    {
    	this.denominator++;
    }
    

    public byte byteValue() {
        return (byte) this.doubleValue();
    }

    public double doubleValue() {
        return ((double) numerator)/((double) denominator);
    }

    public float floatValue() {
        return (float) this.doubleValue();
    }

    public int intValue() {
        return (int) this.doubleValue();
    }

    public long longValue() {
        return (long) this.doubleValue();
    }

    public short shortValue() {
        return (short) this.doubleValue();
    }

    public boolean equals(FractionNonidentical frac) {
        return this.compareTo(frac) == 0;
    }

    public int compareTo(FractionNonidentical frac) {
        long t = this.getNumerator() * frac.getDenominator();
        long f = frac.getNumerator() * this.getDenominator();
        int result = 0;
        if(t>f) {
            result = 1;
        }
        else if(f>t) {
            result = -1;
        }
        else
        {
        	// Equal values, sort by denominator
        	result = (new Integer(this.getDenominator())).compareTo(new Integer(frac.getDenominator()));
        }
        return result;
    }
}