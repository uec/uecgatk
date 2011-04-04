package edu.usc.epigenome.uecgatk;

/**
 * @author benb
 *
 * Representation of a fraction where 1/2 does not equal 2/4
 */
public class FractionNonidentical extends Number implements Comparable<Number> {
    protected int numerator;
    protected int denominator;
    
    protected boolean useIdentical = false;

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

    public boolean isUseIdentical() {
		return useIdentical;
	}

	public void setUseIdentical(boolean useIdentical) {
		this.useIdentical = useIdentical;
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
    
    
    @Override
	public String toString() {
		String out = String.format("%d/%d", this.getNumerator(), this.getDenominator());
		return out;
	}

    @Override
	public byte byteValue() {
        return (byte) this.doubleValue();
    }

    @Override
    public double doubleValue() {
        return ((double) numerator)/((double) denominator);
    }

    @Override
    public float floatValue() {
        return (float) this.doubleValue();
    }

    @Override
    public int intValue() {
        return (int) this.doubleValue();
    }

    @Override
    public long longValue() {
        return (long) this.doubleValue();
    }

    @Override
    public short shortValue() {
        return (short) this.doubleValue();
    }

	@Override
	public boolean equals(Object obj) {
		if (this.getClass().isAssignableFrom(obj.getClass()))
		{
			//System.err.println("Equals 1");
			return (this.compareTo((FractionNonidentical)obj) == 0);
		}
		else
		{
			System.err.println("Equals 2");
	       return super.equals(obj);
		}
	       
	}

	
	
	@Override
	public int hashCode() {
		return (new Float(this.floatValue())).hashCode();
	}

	@Override
	public int compareTo(Number num) 
	{
		int result = 0;
	
		if (this.getClass().isAssignableFrom(num.getClass()))
		{
			FractionNonidentical frac = (FractionNonidentical)num;
			
			long t = this.getNumerator() * frac.getDenominator();
			long f = frac.getNumerator() * this.getDenominator();
			
			if(t>f) {
				result = 1;
			}
			else if(f>t) {
				result = -1;
			}
			else
			{
				// Equal values, sort by denominator
				// Fraction value is identical .  what do we want to do
				if (!this.useIdentical)
				{
					result = (new Integer(this.getDenominator())).compareTo(new Integer(frac.getDenominator()));
				}
			}
			//System.err.printf("compareTo(%s,%s)=%d\n",this, frac, result);
		}
		else
		{
				result = (new Float(this.floatValue())).compareTo(new Float(num.floatValue()));
		}

		return result;	
	}


    
}