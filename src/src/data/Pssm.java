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

import io.FReader;

import java.io.File;

import util.ErrorUtils;
import util.Mappings;

public class Pssm {
	
	
	private int 	length 	= 0;
	private int[][] pssm 	= null; //matrix with the PSSM scores
	
	
	private Pssm(String pssmFile, int seqLength)
	{
		this.length = seqLength;
		this.pssm 	= new int[seqLength][20];
		
		if (!this.readPssmFile(pssmFile)) {this.length = -1;}
	}
	
	
	public static Pssm newPssm(String pssmFile, int seqLength)
	{
		if ((new File(pssmFile)).exists())
		{
			Pssm pssm = new Pssm(pssmFile, seqLength);
			
			if (pssm.getLength() > 0)
			{
				return pssm;
			}
			else
			{
				return null;
			}
		}
		else
		{
			ErrorUtils.printError(Pssm.class, "Could not find pssm file: " + pssmFile, null);
			
			return null;
		}
	}
	
	
	private boolean readPssmFile(String pssmFile)
	{
		int lines = 0;
		
		FReader.openFile(pssmFile);
		
		String line = FReader.readLine();
		
		while (line != null)
		{
			line = line.trim();
			
			String[] content = line.split("\\s+");
			
			if (content.length == 40) //Matrix header
			{
				int[] newIndexPssm = new int[20];
				
				for (int i = 0; i < 20; ++i)
				{
					char 	aa 		= content[i].charAt(0);
					int 	index 	= Mappings.aaToInt(aa);
					
					if (index != Mappings.defaultValue)
					{
						newIndexPssm[i] = index;
					}
					else
					{
						ErrorUtils.printError(this.getClass(), "Malformed pssm matrix header in " + pssmFile, null);
						
						FReader.closeFile();
						
						return false;
					}
				}
				
				line 	= FReader.readLine();
				line 	= line.trim();
				content = line.split("\\s+");
				
				while (content.length == 44)
				{
					int pos = Integer.valueOf(content[0]) - 1;
					
					for (int i = 0; i < 20; ++i)
					{
						int value = Integer.valueOf(content[i+2]);
						
						this.pssm[pos][newIndexPssm[i]] = value;
					}
					
					line 	= FReader.readLine();
					line 	= line.trim();
					content = line.split("\\s+");
					
					++lines;
				}
				
				break;
			}
			
			line = FReader.readLine();
		}
		
		FReader.closeFile();
		
		if (lines != this.length)
		{
			ErrorUtils.printError(Pssm.class, "Sequence length and pssm file size do not match for " + pssmFile, null);
			
			return false;
		}
		
		return true;
	}
	
	
	public int getLength()
	{
		return this.length;
	}
	
	
	public int getScore(int pos, int aa)
	{
		return this.pssm[pos][aa];
	}

}
