//    <TMSEG: Prediction of Transmembrane Helices in Proteins.>
//    Copyright (C) 2014  Michael Bernhofer
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.

package data;

import java.util.Arrays;

import util.ErrorUtils;

public class Protein {
	
	
	private String 		name 		= null;
	private String 		header 		= null;
	private char[] 		sequence 	= null;
	private char[] 		structure 	= null;
	private Pssm 		pssm 		= null;
	
	private char[] 		prediction 	= null;
	private int[] 		confidence 	= null;
	private int[] 		solRaw 		= null;
	private int[] 		tmhRaw 		= null;
	private int[] 		sigRaw 		= null;
	private int[] 		segmentRaw	= null;
	private int 		topologyRaw = -1;
	
	private boolean 	isRealTmp 	= false;
	private boolean 	isPredTmp 	= false;
	private boolean 	hasRealSigP = false;
	private boolean 	hasPredSigP = false;
	
	private Protein(String name, String header, String sequence, String structure)
	{
		this.name 			= name;
		this.header 		= header;
		this.sequence 		= sequence.toCharArray();
		this.structure 		= structure.toCharArray();
		this.prediction 	= new char[this.structure.length];
		this.isRealTmp 		= structure.contains("H") || structure.contains("h");
		this.hasRealSigP 	= structure.contains("S") || structure.contains("s");
	}
	
	
	public static Protein newProtein(String name, String header, String sequence, String structure)
	{
		//check mandatory data
		if (name == null || header == null || sequence == null)
		{
			ErrorUtils.printError(Protein.class, "Insufficient information for protein " + name + ".", null);
			
			return null;
		}
		
		//check annotated structure
		if (structure != null)
		{
			//check sequence and structure length
			if (sequence.length() != structure.length())
			{
				ErrorUtils.printError(Protein.class, "Sequence and structure for protein " + name + " do not match.", null);
				
				return null;
			}
			
			//replace signalP symbols
			structure = structure.replaceAll("\\.", "N");
			structure = structure.replaceAll("t", "h");
			structure = structure.replaceAll("T", "H");
		}
		else
		{
			char[] str 	= new char[sequence.length()];
			
			Arrays.fill(str, 'N');
			
			structure = new String(str);
		}
		
		return new Protein(name, header, sequence, structure);
	}
	
	
	public void setPssm(Pssm pssm)
	{
		this.pssm = pssm;
	}
	
	
	public void setPredTmp(boolean isTmp)
	{
		this.isPredTmp = isTmp;
	}
	
	
	public void setPredSigP(boolean hasSigP)
	{
		this.hasPredSigP = hasSigP;
	}
	
	
	public void setPrediction(char[] prediction)
	{
		if (prediction.length == this.sequence.length)
		{
			this.prediction = prediction;
		}
	}
	
	
	public void setConfidence(int[] confidence)
	{
		if (confidence.length == this.sequence.length)
		{
			this.confidence = confidence;
		}
	}
	
	
	public void setSolRaw(int[] solRaw)
	{
		if (solRaw.length == this.sequence.length)
		{
			this.solRaw = solRaw;
		}
	}
	
	
	public void setTmhRaw(int[] tmhRaw)
	{
		if (tmhRaw.length == this.sequence.length)
		{
			this.tmhRaw = tmhRaw;
		}
	}
	
	
	public void setSigRaw(int[] sigRaw)
	{
		if (sigRaw.length == this.sequence.length)
		{
			this.sigRaw = sigRaw;
		}
	}
	
	
	public void setSegmentRaw(int[] segmentRaw)
	{
		if (segmentRaw.length == this.sequence.length)
		{
			this.segmentRaw = segmentRaw;
		}
	}
	
	
	public void setTopologyRaw(int topologyRaw)
	{
		this.topologyRaw = topologyRaw;
	}
	
	
	public String getName()
	{
		return this.name;
	}
	
	
	public String getHeader()
	{
		return this.header;
	}
	
	
	public char[] getSequence()
	{
		return this.sequence;
	}
	
	
	public char[] getStructure()
	{
		return this.structure;
	}
	
	
	public Pssm getPssm()
	{
		return this.pssm;
	}
	
	
	public char[] getPrediction()
	{
		return this.prediction;
	}
	
	
	public int[] getConfidence()
	{
		return this.confidence;
	}
	
	
	public int[] getSolRaw()
	{
		return this.solRaw;
	}
	
	
	public int[] getTmhRaw()
	{
		return this.tmhRaw;
	}
	
	
	public int[] getSigRaw()
	{
		return this.sigRaw;
	}
	
	
	public int[] getSegmentRaw()
	{
		return this.segmentRaw;
	}
	
	
	public int getTopologyRaw()
	{
		return this.topologyRaw;
	}
	
	
	public boolean isRealTmp()
	{
		return this.isRealTmp;
	}
	
	
	public boolean isPredTmp()
	{
		return this.isPredTmp;
	}
	
	
	public boolean hasRealSigP()
	{
		return this.hasRealSigP;
	}
	
	
	public boolean hasPredSigP()
	{
		return this.hasPredSigP;
	}

}
