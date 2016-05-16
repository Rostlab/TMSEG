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

package util;


/**
 * Maps amino acids and secondary structures to indices and back.
 * Includes physicochemical properties of amino acids.
 */
public class Mappings {
	
	public static final int defaultValue 	= 255;
	
	public static final int indexNotTmh 	= 0;
	public static final int indexTmh 		= 1;
	public static final int indexLoop 		= 2;
	public static final int indexSignal 	= 3;
	public static final int indexUnknown 	= 4;
	public static final int indexInside 	= 5;
	public static final int indexOutside 	= 6;
	
	
	public static int aaToInt(char aminoAcid)
	{
		aminoAcid = Character.toUpperCase(aminoAcid);
		
		switch (aminoAcid)
		{
			case 'A': 	return 0;
			case 'C': 	return 1;
			case 'D': 	return 2;
			case 'E': 	return 3;
			case 'F': 	return 4;
			case 'G': 	return 5;
			case 'H': 	return 6;
			case 'I': 	return 7;
			case 'K': 	return 8;
			case 'L': 	return 9;
			case 'M': 	return 10;
			case 'N': 	return 11;
			case 'P': 	return 12;
			case 'Q': 	return 13;
			case 'R': 	return 14;
			case 'S': 	return 15;
			case 'T': 	return 16;
			case 'V': 	return 17;
			case 'W': 	return 18;
			case 'Y': 	return 19;
			case '-': 	return 20;
			default: 	return defaultValue;
		}
	}
	
	
	public static char intToAa(int aminoAcid)
	{
		switch (aminoAcid)
		{
			case 0: 	return 'A';
			case 1: 	return 'C';
			case 2: 	return 'D';
			case 3: 	return 'E';
			case 4: 	return 'F';
			case 5: 	return 'G';
			case 6: 	return 'H';
			case 7: 	return 'I';
			case 8: 	return 'K';
			case 9: 	return 'L';
			case 10: 	return 'M';
			case 11: 	return 'N';
			case 12: 	return 'P';
			case 13: 	return 'Q';
			case 14: 	return 'R';
			case 15: 	return 'S';
			case 16: 	return 'T';
			case 17: 	return 'V';
			case 18: 	return 'W';
			case 19: 	return 'Y';
			case 20: 	return '-';
			default: 	return 'X';
		}
	}
	
	
	public static int ssToInt(char secondaryStructure)
	{
		secondaryStructure = Character.toUpperCase(secondaryStructure);
		
		switch (secondaryStructure)
		{
			case '1': 	return indexNotTmh;
			case '2': 	return indexNotTmh;
			case 'N': 	return indexNotTmh;
			case 'H': 	return indexTmh;
			case 'L': 	return indexLoop;
			case 'S': 	return indexSignal;
			case 'U': 	return indexUnknown;
			case ' ': 	return indexUnknown;
			default: 	return indexNotTmh;
		}
	}
	
	
	public static char intToSs(int secondaryStructure)
	{
		switch (secondaryStructure)
		{
			case indexNotTmh: 	return 'N'; //check Protein class if changed!
			case indexTmh: 		return 'H'; //check Protein class if changed!
			case indexSignal: 	return 'S'; //check Protein class if changed!
			default: 			return 'X';
		}
	}
	
	
	public static int topToInt(char secondaryStructure)
	{
		secondaryStructure = Character.toUpperCase(secondaryStructure);
		
		switch (secondaryStructure)
		{
			case '1': 	return indexInside;
			case '2': 	return indexOutside;
			default: 	return defaultValue;
		}
	}
	
	
	public static char intToTop(int secondaryStructure)
	{
		switch (secondaryStructure)
		{
			case indexInside: 	return '1';
			case indexOutside: 	return '2';
			default: 			return 'X';
		}
	}
	
	
	public static double hydrophobicity(char aminoAcid)
	{
		aminoAcid = Character.toUpperCase(aminoAcid);
		
		switch (aminoAcid)
		{
			case 'A': 	return 1.8;
			case 'C': 	return 2.5;
			case 'D': 	return -3.5;
			case 'E': 	return -3.5;
			case 'F': 	return 2.8;
			case 'G': 	return -0.4;
			case 'H': 	return -3.2;
			case 'I': 	return 4.5;
			case 'K': 	return -3.9;
			case 'L': 	return 3.8;
			case 'M': 	return 1.9;
			case 'N': 	return -3.5;
			case 'P': 	return -1.6;
			case 'Q': 	return -3.5;
			case 'R': 	return -4.5;
			case 'S': 	return -0.8;
			case 'T': 	return -0.7;
			case 'V': 	return 4.2;
			case 'W': 	return -0.9;
			case 'Y': 	return -1.3;
			default: 	return defaultValue;
		}
	}
	
	
	public static int charge(char aminoAcid)
	{
		aminoAcid = Character.toUpperCase(aminoAcid);
		
		switch (aminoAcid)
		{
			case 'A': 	return 0;
			case 'C': 	return 0;
			case 'D': 	return -1;
			case 'E': 	return -1;
			case 'F': 	return 0;
			case 'G': 	return 0;
			case 'H': 	return 0;
			case 'I': 	return 0;
			case 'K': 	return 1;
			case 'L': 	return 0;
			case 'M': 	return 0;
			case 'N': 	return 0;
			case 'P': 	return 0;
			case 'Q': 	return 0;
			case 'R': 	return 1;
			case 'S': 	return 0;
			case 'T': 	return 0;
			case 'V': 	return 0;
			case 'W': 	return 0;
			case 'Y': 	return 0;
			default: 	return defaultValue;
		}
	}
	
	
	public static int polarity(char aminoAcid)
	{
		aminoAcid = Character.toUpperCase(aminoAcid);
		
		switch (aminoAcid)
		{
			case 'A': 	return 0;
			case 'C': 	return 1;
			case 'D': 	return 1;
			case 'E': 	return 1;
			case 'F': 	return 0;
			case 'G': 	return 0;
			case 'H': 	return 1;
			case 'I': 	return 0;
			case 'K': 	return 1;
			case 'L': 	return 0;
			case 'M': 	return 0;
			case 'N': 	return 1;
			case 'P': 	return 0;
			case 'Q': 	return 1;
			case 'R': 	return 1;
			case 'S': 	return 1;
			case 'T': 	return 1;
			case 'V': 	return 0;
			case 'W': 	return 0;
			case 'Y': 	return 1;
			default: 	return defaultValue;
		}
	}

}
