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

import java.util.ArrayList;


public class FastaReader {
	
	
	/**
	 * Read (multi-line) FASTA files
	 * 
	 * @param filename
	 * @return
	 */
	public static ArrayList<Protein> readFastaFile(String filename)
	{
		ArrayList<Protein> 	proteins 	= new ArrayList<Protein>();
		StringBuilder 		sequence 	= new StringBuilder();
		String[] 			content 	= null;
		String 				line 		= null;
		String 				name 		= null;
		String 				header 		= null;
		
		FReader.openFile(filename);
		
		line = FReader.readLine();
		
		while (line != null)
		{
			if (line.startsWith(">"))
			{
				line 	= line.trim();
				header 	= line;
				content = line.split("\\s");
				name 	= content[0].replaceFirst(">", "").trim().split("\\|")[0].trim();
				line 	= FReader.readLine();
				
				while (line != null && !line.startsWith(">"))
				{
					sequence.append(line.trim());
					
					line = FReader.readLine();
				}
				
				Protein protein = Protein.newProtein(name, header, sequence.toString().replaceAll("\\s", ""), null);
				
				if (protein != null) {proteins.add(protein);}
				
				sequence = new StringBuilder();
			}
			else
			{
				line = FReader.readLine();
			}
		}
		
		FReader.closeFile();
		
		return proteins;
	}
	
	
	/**
	 * Read (multi-line) structure files
	 * 
	 * @param filename
	 * @return
	 */
	public static ArrayList<Protein> readStructureFile(String filename)
	{
		ArrayList<Protein> 	proteins 	= new ArrayList<Protein>();
		String[] 			content 	= null;
		String 				line 		= null;
		String 				name 		= null;
		String 				header 		= null;
		String 				sequence 	= null;
		String 				structure 	= null;
		
		FReader.openFile(filename);
		
		line = FReader.readLine();
		
		while (line != null)
		{
			if (line.startsWith(">"))
			{
				line 		= line.trim();
				header 		= line;
				content 	= line.split("\\s");
				name 		= content[0].replaceFirst(">", "").trim().split("\\|")[0].trim();
				sequence 	= FReader.readLine();
				structure 	= FReader.readLine();
				
				if (sequence != null) 	{sequence = sequence.trim().replaceAll("\\s", "");}
				if (structure != null) 	{structure = structure.trim().replaceAll("\\s", "");}
				
				Protein protein = Protein.newProtein(name, header, sequence, structure);
				
				if (protein != null) {proteins.add(protein);}
			}
			
			line = FReader.readLine();
		}
		
		FReader.closeFile();
		
		return proteins;
	}

}
