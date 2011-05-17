package edu.usc.epigenome.uecgatk.bisulfitesnpmodel;


public class BaseUtilsMore {

	public final static byte A = (byte)'A';
    public final static byte C = (byte)'C';
    public final static byte G = (byte)'G';
    public final static byte T = (byte)'T';
    public final static byte R = (byte)'R'; //A,G
    public final static byte Y = (byte)'Y'; //C,T
    public final static byte S = (byte)'S'; //G,C
    public final static byte W = (byte)'W'; //A,T
    public final static byte K = (byte)'K'; //G,T
    public final static byte M = (byte)'M'; //A,C
    public final static byte B = (byte)'B'; //C,G,T
    public final static byte H = (byte)'H'; //A,C,T
    public final static byte V = (byte)'V'; //A,C,G

    public final static byte N = (byte)'N';
    public final static byte D = (byte)'D'; //A,G,T
	
	
	static public byte[] toUpperCase(byte[] in)
	{
		byte[] out = new byte[in.length];
		for (int i = 0; i < in.length; i++)
		{
			out[i] = toUpperCase(in[i]);
		}
		return out;
	}

	static public byte toUpperCase(byte b)
	{
		return (byte)Character.toUpperCase((char)b);
	}
	
	static public boolean iupacCodeEqual(byte pattern, byte observe, boolean negStrand){
		pattern = toUpperCase(pattern);
		observe = toUpperCase(observe);
		switch(observe){
			case 'A':
				switch(pattern){
					case 'A':
					case 'R':
					case 'W':
					case 'M':
					case 'H':
					case 'V':
					case 'D':
					case 'N':
						return true;
					case 'G':
					case 'S':
					case 'K':
					case 'B':
						if(negStrand){
							return true;
						}	
						else{
							return false;
						}				
					default:
						return false;
				}
				
			case 'C':
				switch(pattern){
					case 'C':
					case 'Y':
					case 'S':
					case 'M':
					case 'H':
					case 'V':
					case 'B':
					case 'N':
						return true;	
					default:
						return false;
				}
				
			case 'G':
				switch(pattern){
					case 'G':
					case 'R':
					case 'S':
					case 'K':
					case 'D':
					case 'V':
					case 'B':
					case 'N':
						return true;	
					default:
						return false;
				}
				
			case 'T':
				switch(pattern){
					case 'T':
					case 'W':
					case 'K':
					case 'D':
					case 'Y':
					case 'H':
					case 'B':
					case 'N':
						return true;
					case 'C':
					case 'S':
					case 'M':
					case 'V':
						if(negStrand){
							return false;
						}	
						else{
							return true;
						}		
				default:
					return false;
				}
				
			case 'N':
				return true;
			default:
				System.err.println("error! wrong observed base!");
				return false;
		}
		
	}
	
	static public byte iupacCodeComplement(byte base) {
		base = toUpperCase(base);
		switch (base) {
        	case 'A':
        		return 'T';
        	case 'C':
        		return 'G';
        	case 'G':
        		return 'C';
        	case 'T':
        		return 'A';
        	case 'R':
        		return 'Y';
        	case 'Y':
        		return 'R';
        	case 'S':
        		return 'S';
        	case 'W':
        		return 'W';
        	case 'K':
        		return 'M';
        	case 'M':
        		return 'K';
        	case 'B':
        		return 'V';
        	case 'H':
        		return 'D';
        	case 'D':
        		return 'H';
        	case 'V':
        		return 'B';
            default: return base;
        }
    }
	
}
