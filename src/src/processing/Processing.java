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

package processing;

import java.util.Arrays;
import java.util.List;

import util.Mappings;
import data.Protein;


/**
 * Class for the post-processing of the predicted protein structure.
 * Builds the consensus over all prediction methods and removes 
 * helices which are too short.
 */
public class Processing {
	
	
	/**
	 * Processes a list of proteins.
	 * 
	 * @param proteins
	 * @param onlyTmp
	 * @param minLength
	 * @param wSol
	 * @param wTmh
	 * @param wSig
	 */
	public static void process(List<Protein> proteins, boolean onlyTmp, int minLength, int wSol, int wTmh, int wSig)
	{
		for (Protein protein : proteins)
		{
			process(protein, onlyTmp, minLength, wSol, wTmh, wSig);
		}
	}
	
	
	/**
	 * Processes a protein.
	 * 
	 * @param protein
	 * @param onlyTmp
	 * @param minLength
	 * @param wSol
	 * @param wTmh
	 * @param wSig
	 */
	public static void process(Protein protein, boolean onlyTmp, int minLength, int wSol, int wTmh, int wSig)
	{
		if (protein == null || protein.getPrediction() == null) {return;}
		
		if (onlyTmp && !protein.isPredTmp()) {return;}
		
		processPrediction(protein, minLength, wSol, wTmh, wSig);
	}
	
	
	/**
	 * Generates the structure prediction based on the raw output of the HelixIndexer
	 * and the given probability penalties for the three structural states:
	 * soluble, transmembrane, and signal peptide.
	 * 
	 * @param protein
	 * @param minLength
	 * @param wSol
	 * @param wTmh
	 * @param wSig
	 */
	private static void processPrediction(Protein protein, int minLength, int wSol, int wTmh, int wSig)
	{
		int[] 	scoresSol 	= protein.getSolRaw();
		int[] 	scoresTmh 	= protein.getTmhRaw();
		int[] 	scoresSig 	= protein.getSigRaw();
		char[] 	prediction 	= protein.getPrediction();
		
		if (scoresSol == null || scoresTmh == null || scoresSig == null) {return;}
		
		scoresSol = medianFilter(scoresSol, 5);
		scoresTmh = medianFilter(scoresTmh, 5);
		scoresSig = medianFilter(scoresSig, 5);
		
		int 	lastSPIndex 		= 0;
		boolean hasSignalPeptide 	= false;
		
		for (int i = 0; i < prediction.length; ++i)
		{
			int sSol = scoresSol[i] - wSol;
			int sTmh = scoresTmh[i] - wTmh;
			int sSig = scoresSig[i] - wSig;
			
			if (sSol >= sTmh && sSol >= sSig)
			{
				prediction[i] = Mappings.intToSs(Mappings.indexNotTmh);
			}
			else if (sTmh >= sSig)
			{
				prediction[i] = Mappings.intToSs(Mappings.indexTmh);
			}
			else
			{
				lastSPIndex 		= i;
				hasSignalPeptide 	= true;
				
				prediction[i] = Mappings.intToSs(Mappings.indexSignal);
			}
		}
		
		protein.setPredTmp(cutShortHelices(protein, minLength));
		
		if (hasSignalPeptide)
		{
			protein.setPredSigP(processSignalPeptide(prediction, lastSPIndex));
		}
		
		protein.setPrediction(prediction);
	}
	
	
	/**
	 * Refines the signal peptide prediction by removing short 
	 * signal peptides (less than 4 consecutive residues).
	 * 
	 * @param prediction
	 * @param lastSPIndex
	 * @return
	 */
	private static boolean processSignalPeptide(char[] prediction, int lastSPIndex)
	{
		boolean hasSignalPeptide = false;
		
		for (int i = 0; i <= lastSPIndex; ++i)
		{
			char type = prediction[i];
			
			if (Mappings.ssToInt(type) == Mappings.indexSignal)
			{
				int start = i;
				
				while (i <= lastSPIndex && prediction[i] == type) {++i;}
				
				--i;
				
				if ((i - start + 1) >= 4)
				{
					hasSignalPeptide = true;
					
					break;
				}
			}
		}
		
		if (hasSignalPeptide)
		{
			for (int i = 0; i <= lastSPIndex; ++i)
			{
				prediction[i] = Mappings.intToSs(Mappings.indexSignal);
			}
		}
		else
		{
			for (int i = 0; i <= lastSPIndex; ++i)
			{
				if (Mappings.ssToInt(prediction[i]) == Mappings.indexSignal)
				{
					prediction[i] = Mappings.intToSs(Mappings.indexNotTmh);
				}
			}
		}
		
		return hasSignalPeptide;
	}
	
	
	/**
	 * Removes all predicted TMHs which are shorter than the minimal size.
	 * 
	 * @param protein
	 * @param minLength
	 */
	private static boolean cutShortHelices(Protein protein, int minLength)
	{
		boolean isTmp = false;
		
		char[] prediction = protein.getPrediction();
		
		for (int i = 0; i < prediction.length; ++i)
		{
			if (Mappings.ssToInt(prediction[i]) == Mappings.indexTmh)
			{
				int 	start 	= i;
				char 	type 	= prediction[i];
				
				while (i < prediction.length && prediction[i] == type) {++i;}
				
				int length = i - start;
				
				if (length < minLength)
				{
					for (int j = start; j < i; ++j)
					{
						prediction[j] = Mappings.intToSs(Mappings.indexNotTmh);
					}
				}
				else
				{
					isTmp = true;
				}
				
				--i;
			}
		}
		
		return isTmp;
	}
	
	
	/**
	 * Checks for a list of proteins if these proteins are to be
	 * classified as transmembrane or not.
	 * 
	 * @param proteins
	 */
	public static void tmpCheck(List<Protein> proteins)
	{
		for (Protein protein : proteins)
		{
			tmpCheck(protein);
		}
	}
	
	
	/**
	 * Checks for a protein if this protein is to be
	 * classified as transmembrane or not.
	 * 
	 * @param protein
	 */
	public static void tmpCheck(Protein protein)
	{
		if (protein == null || protein.getPrediction() == null) {return;}
		
		boolean isTmp = isTmp(protein);
		
		protein.setPredTmp(isTmp);
	}
	
	
	/**
	 * Searches a protein for TMH-residues and returns if a TMH was found.
	 * 
	 * @param protein
	 * @return
	 */
	private static boolean isTmp(Protein protein)
	{
		char[] prediction = protein.getPrediction();
		
		for (int i = 0; i < prediction.length; ++i)
		{
			if (Mappings.ssToInt(prediction[i]) == Mappings.indexTmh) {return true;}
		}
		
		return false;
	}
	
	
	/**
	 * Tests if the topology is consistent.
	 * 
	 * @param structure
	 * @return
	 */
	public static boolean testTopology(char[] structure)
	{
		if (structure == null) {return false;}
		
		int 	firstTop = Mappings.defaultValue;
		boolean switched = false;
		
		for (int i = 0; i < structure.length; ++i)
		{
			char 	type 		= structure[i];
			int 	ssIndex 	= Mappings.ssToInt(type);
			int 	topIndex 	= Mappings.topToInt(type);
			
			if (firstTop == Mappings.defaultValue && topIndex != Mappings.defaultValue)
			{
				firstTop = topIndex;
			}
			else if (firstTop != Mappings.defaultValue && ssIndex == Mappings.indexTmh)
			{
				switched = !switched;
				
				while (i < structure.length && structure[i] == type) {++i;}
				
				--i;
			}
			else if (topIndex != Mappings.defaultValue)
			{
				if ((switched && topIndex == firstTop) || (!switched && topIndex != firstTop)) {return false;}
			}
		}
		
		return true;
	}
	
	
	/**
	 * Extrapolates the topology for unknown regions.
	 * 
	 * @param structure
	 * @return
	 */
	public static char[] extrapolateTopology(char[] structure)
	{
		if (structure == null) {return structure;}
		
		char[] 	newStructure 	= structure.clone();
		int 	firstTop 		= Mappings.defaultValue;
		int 	firstPos 		= 0;
		
		for (int i = 0; i < newStructure.length; ++i)
		{
			int topIndex = Mappings.topToInt(newStructure[i]);
			
			if (firstTop == Mappings.defaultValue && topIndex != Mappings.defaultValue)
			{
				firstTop = topIndex;
				firstPos = i;
			}
		}
		
		boolean switched = false;
		
		for (int i = firstPos; i >= 0; --i)
		{
			char 	type 		= newStructure[i];
			int 	ssIndex 	= Mappings.ssToInt(type);
			int 	topIndex 	= Mappings.topToInt(type);
			
			if (ssIndex != Mappings.indexUnknown)
			{
				if (ssIndex == Mappings.indexTmh)
				{
					switched = !switched;
					
					while (i >= 0 && newStructure[i] == type) {--i;}
					
					++i;
				}
				else if (ssIndex == Mappings.indexNotTmh && topIndex == Mappings.defaultValue)
				{
					if (!switched)
					{
						newStructure[i] = Mappings.intToTop(firstTop);
					}
					else if (firstTop == Mappings.indexInside)
					{
						newStructure[i] = Mappings.intToTop(Mappings.indexOutside);
					}
					else if (firstTop == Mappings.indexOutside)
					{
						newStructure[i] = Mappings.intToTop(Mappings.indexInside);
					}
				}
			}
		}
		
		switched = false;
		
		for (int i = firstPos; i < newStructure.length; ++i)
		{
			char 	type 		= newStructure[i];
			int 	ssIndex 	= Mappings.ssToInt(type);
			int 	topIndex 	= Mappings.topToInt(type);
			
			if (ssIndex != Mappings.indexUnknown)
			{
				if (ssIndex == Mappings.indexTmh)
				{
					switched = !switched;
					
					while (i < newStructure.length && newStructure[i] == type) {++i;}
					
					--i;
				}
				else if (ssIndex == Mappings.indexNotTmh && topIndex == Mappings.defaultValue)
				{
					if (!switched)
					{
						newStructure[i] = Mappings.intToTop(firstTop);
					}
					else if (firstTop == Mappings.indexInside)
					{
						newStructure[i] = Mappings.intToTop(Mappings.indexOutside);
					}
					else if (firstTop == Mappings.indexOutside)
					{
						newStructure[i] = Mappings.intToTop(Mappings.indexInside);
					}
				}
			}
		}
		
		return newStructure;
	}
	
	
	/**
	 * Assigns confidence scores to all predicted TMH within a list of proteins.
	 * 
	 * @param proteins
	 */
	public static void assignConfidence(List<Protein> proteins)
	{
		for (Protein protein : proteins)
		{
			assignConfidence(protein);
		}
	}
	
	
	/**
	 * Assigns confidence scores to all predicted TMH within a protein.
	 * 
	 * @param protein
	 */
	public static void assignConfidence(Protein protein)
	{
		if (protein == null) {return;}
		
		int[] 	tmhScores 	= protein.getTmhRaw();
		int[] 	solScores 	= protein.getSolRaw();
		int[] 	segScores 	= protein.getSegmentRaw();
		char[] 	prediction 	= protein.getPrediction();
		
		if (tmhScores == null || segScores == null || prediction == null) {return;}
		
		int[] confidence = new int[prediction.length];
		
		for (int i = 0; i < prediction.length; ++i)
		{
			char type = prediction[i];
			
			if (Mappings.ssToInt(type) == Mappings.indexTmh)
			{
				int conf 	= 0;
				int start 	= i;
				
				while (i < prediction.length && prediction[i] == type)
				{
					conf += Math.max((tmhScores[i] - solScores[i] + 125) , segScores[i]);
					
					++i;
				}
				
				int end = i;
				
				conf = conf / (end - start);
				conf = conf / 100;
				
				for (int j = start; j < end; ++j)
				{
					confidence[j] = Math.min(9, conf);
				}
				
				--i;
			}
		}
		
		protein.setConfidence(confidence);
	}
	
	
	/**
	 * Runs a median filter over an array with a given window size.
	 * 
	 * @param scores
	 * @param windowSize
	 * @return
	 */
	public static int[] medianFilter(int[] scores, int windowSize)
	{
		if (scores == null || scores.length == 0) {return scores;}
		
		int 	offset 		= windowSize / 2;
		int[] 	newScores 	= new int[scores.length];
		
		newScores[0] = scores[0];
		
		for (int i = 1; i < offset; ++i)
		{
			int[] window = Arrays.copyOfRange(scores, 0, i+i+1);
			
			Arrays.sort(window);
			
			newScores[i] = window[i];
		}
		
		for (int i = offset; i < scores.length-offset; ++i)
		{
			int[] window = Arrays.copyOfRange(scores, i-offset, i+offset+1);
			
			Arrays.sort(window);
			
			newScores[i] = window[offset];
		}
		
		for (int i = 1; i < offset; ++i)
		{
			int 	j 		= scores.length-(i+1);
			int[] 	window 	= Arrays.copyOfRange(scores, j-i, scores.length);
			
			Arrays.sort(window);
			
			newScores[j] = window[i];
		}
		
		newScores[newScores.length-1] = scores[scores.length-1];
		
		return newScores;
	}

}
