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

package predictors;

import io.ModelHandler;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import util.ErrorUtils;
import util.Globals;
import util.Mappings;
import weka.classifiers.Classifier;
import weka.classifiers.functions.MultilayerPerceptron;
import weka.core.Attribute;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.SparseInstance;
import weka.core.converters.ArffSaver;
import data.Protein;
import data.Pssm;


/**
 * Class to predict transmembrane helices within a protein.
 */
public class HelixPredictor {
	
	
	private ArrayList<Attribute> 	attributes 		= null;
	private Instances 				dataset 		= null;
	private Classifier 				classifier 		= null;
	private boolean 				isTrained 		= false;
	
	
	public HelixPredictor()
	{
		this.buildAttributes();
		this.initialize();
	}
	
	
	/**
	 * Initializes the attributes list for the Weka arff-format.
	 */
	private void buildAttributes()
	{
		this.attributes = new ArrayList<Attribute>();
		
		
		//helix features
		for (int i = 0; i < 20; ++i)
		{
			this.attributes.add(new Attribute("conserved_composition_" + Mappings.intToAa(i)));
		}
		
		for (int i = 0; i < 20; ++i)
		{
			this.attributes.add(new Attribute("non-conserved_composition_" + Mappings.intToAa(i)));
		}
		
		this.attributes.add(new Attribute("length"));
		
		this.attributes.add(new Attribute("conserved_hydrophobicity"));
		this.attributes.add(new Attribute("non-conserved_hydrophobicity"));
		
		this.attributes.add(new Attribute("conserved_hydrophobic"));
		this.attributes.add(new Attribute("non-conserved_hydrophobic"));
		
		this.attributes.add(new Attribute("conserved_charged"));
		this.attributes.add(new Attribute("non-conserved_charged"));
		
		ArrayList<String> classes = new ArrayList<String>();
		
		classes.add(String.valueOf(Mappings.indexNotTmh));
		classes.add(String.valueOf(Mappings.indexTmh));
		
		this.attributes.add(new Attribute("class", classes));
	}
	
	
	/**
	 * Initializes the classifier and dataset.
	 */
	public void initialize()
	{
		this.isTrained 	= false;
		this.classifier = null;
		this.dataset 	= new Instances("HelixPredictor Model", this.attributes, 0);
		
		this.dataset.setClassIndex(this.attributes.size()-1);
	}
	
	
	/**
	 * Inputs a given list of proteins for the training data.
	 * 
	 * @param proteins
	 */
	public void input(List<Protein> proteins)
	{
		for (Protein protein : proteins)
		{
			this.input(protein);
		}
	}
	
	
	/**
	 * Inputs a given protein for the training data.
	 * 
	 * @param protein
	 */
	public void input(Protein protein)
	{
		if (protein == null) 				{return;}
		if (protein.getStructure() == null) {return;}
		if (protein.getPssm() == null) 		{return;}
		
		Pssm 	pssm 		= protein.getPssm();
		char[] 	structure 	= protein.getStructure();
		
		char 	type 	= Mappings.intToSs(Mappings.defaultValue);
		int 	start 	= 0;
		int 	end 	= 0;
		
		for (int i = 0; i < structure.length; ++i)
		{
			//find consecutive transmembrane helices and non-transmembrane regions
			if (Mappings.ssToInt(structure[i]) == Mappings.indexTmh)
			{
				start 	= i;
				type 	= structure[i];
				
				while (i < structure.length && structure[i] == type)
				{
					++i;
				}
				
				--i;
				
				end = i;
				
				this.addSegmentToDatabse(pssm, start, end, Mappings.indexTmh);
				
				//oversampling to compensate for the negative samples (see below)
				this.addSegmentToDatabse(pssm, start, end, Mappings.indexTmh);
				this.addSegmentToDatabse(pssm, start, end, Mappings.indexTmh);
				this.addSegmentToDatabse(pssm, start, end, Mappings.indexTmh);
				this.addSegmentToDatabse(pssm, start, end, Mappings.indexTmh);
				this.addSegmentToDatabse(pssm, start, end, Mappings.indexTmh);
				this.addSegmentToDatabse(pssm, start, end, Mappings.indexTmh);
				
				//generate incomplete/shifted helices
				if (start < end - 6)
				{
					this.addSegmentToDatabse(pssm, start+6, end, Mappings.indexNotTmh);
					this.addSegmentToDatabse(pssm, start, end-6, Mappings.indexNotTmh);
				}
				
				if (end < structure.length - 6)
				{
					this.addSegmentToDatabse(pssm, start, end+6, Mappings.indexNotTmh);
					this.addSegmentToDatabse(pssm, start+6, end+6, Mappings.indexNotTmh);
				}
				
				if (start >= 6)
				{
					this.addSegmentToDatabse(pssm, start-6, end, Mappings.indexNotTmh);
					this.addSegmentToDatabse(pssm, start-6, end-6, Mappings.indexNotTmh);
				}
			}
			else if (Mappings.ssToInt(structure[i]) != Mappings.indexUnknown)
			{
				start 	= i;
				type 	= structure[i];
				
				while (i < structure.length && structure[i] == type)
				{
					++i;
				}
				
				--i;
				
				end = i;
				
				this.addSegmentToDatabse(pssm, start, end, Mappings.indexNotTmh);
			}
		}
	}
	
	
	/**
	 * Analyzes and optimizes predicted transmembrane helices
	 * within a given list of proteins.
	 * 
	 * @param proteins
	 * @param cutoff
	 */
	public void predict(List<Protein> proteins, double cutoff)
	{
		for (Protein protein : proteins)
		{
			this.predict(protein, cutoff);
		}
	}
	
	
	/**
	 * Analyzes and optimizes predicted transmembrane helices
	 * within a given protein.
	 * 
	 * @param protein
	 * @param cutoff
	 */
	public void predict(Protein protein, double cutoff)
	{
		if (protein == null) 					{return;}
		if (protein.getPssm() == null) 			{return;}
		if (protein.getPrediction() == null) 	{return;}
		
		//analyse transmembrane proteins only (predicted)
		if (!protein.isPredTmp()) {return;}
		
		if (protein.getSegmentRaw() == null) {protein.setSegmentRaw(new int[protein.getPrediction().length]);}
		
		this.splitTMHs(protein, cutoff);
		this.adjustTMHs(protein, cutoff);
		
		for (int i = 0; i < 5; ++i)
		{
			//continue for up to five additional rounds (only adjust after splitting)
			if (!this.splitTMHs(protein, cutoff) || !this.adjustTMHs(protein, cutoff)) {return;}
		}
		
		//final split
		this.splitTMHs(protein, cutoff);
	}
	
	
	/**
	 * Adjusts predicted transmembrane helices within a given protein.
	 * 
	 * @param protein
	 * @param cutoff
	 * @return
	 */
	private boolean adjustTMHs(Protein protein, double cutoff)
	{
		boolean adjust 		= false;
		Pssm 	pssm 		= protein.getPssm();
		char[] 	structure 	= protein.getPrediction();
		int[] 	segmentRaw 	= protein.getSegmentRaw();
		
		for (int i = 0; i < structure.length; ++i)
		{
			try
			{
				if (Mappings.ssToInt(structure[i]) == Mappings.indexTmh)
				{
					int 	start 	= i;
					char 	type 	= structure[i];
					
					//go to end of transmembrane helix
					while (i < structure.length && structure[i] == type) {++i;}
					
					--i;
					
					int end = i;
					
					Instance window = this.buildInstance(pssm, start, end);
					
					window.isMissing((Attribute)this.attributes.get(this.attributes.size()-1));
					window.setDataset(this.dataset);
					
					double[] probabilities = this.classifier.distributionForInstance(window);
					
					double 	bestProb 	= probabilities[Mappings.indexTmh];
					int 	bestStart 	= -1;
					int 	bestEnd 	= -1;
					
					//shift TMH start/end around and find best position
					for (int newStart = start-Globals.PREDICTOR_MAX_SHIFT; newStart <= start+Globals.PREDICTOR_MAX_SHIFT; ++newStart)
					{
						if (newStart < 0) {continue;}
						
						for (int newEnd = end-Globals.PREDICTOR_MAX_SHIFT; newEnd <= end+Globals.PREDICTOR_MAX_SHIFT; ++newEnd)
						{
							if (newEnd >= structure.length) {break;}
							
							window = this.buildInstance(pssm, newStart, newEnd);
							
							window.isMissing((Attribute)this.attributes.get(this.attributes.size()-1));
							window.setDataset(this.dataset);
							
							probabilities = this.classifier.distributionForInstance(window);
							
							if (probabilities[Mappings.indexTmh] > bestProb)
							{
								bestProb 	= probabilities[Mappings.indexTmh];
								bestStart 	= newStart;
								bestEnd 	= newEnd;
							}
						}
					}
					
					//adjust the TMH!
					if (bestProb < cutoff)
					{
						for (int j = start; j <= end; ++j)
						{
							structure[j] 	= Mappings.intToSs(Mappings.indexNotTmh);
							segmentRaw[j] 	= 0;
						}
					}
					else if (bestStart != -1 && bestEnd != -1)
					{
						start 	= Math.min(start, bestStart);
						end 	= Math.max(end, bestEnd);
						
						for (int j = start; j <= end; ++j)
						{
							if (j >= bestStart && j <= bestEnd)
							{
								structure[j] 	= Mappings.intToSs(Mappings.indexTmh);
								segmentRaw[j] 	= (int)(1000 * bestProb);
							}
							else
							{
								structure[j] 	= Mappings.intToSs(Mappings.indexNotTmh);
								segmentRaw[j] 	= 0;
							}
						}
						
						adjust 	= true;
						i 		= end;
					}
					else
					{
						for (int j = start; j <= end; ++j)
						{
							segmentRaw[j] = (int)(1000 * bestProb);
						}
					}
				}
				else
				{
					segmentRaw[i] = 0;
				}
			}
			catch (Exception e)
			{
				ErrorUtils.printError(HelixPredictor.class, "Prediction failed for " + protein.getHeader(), e);
				
				return false;
			}
		}
		
		return adjust;
	}
	
	
	/**
	 * Splits predicted transmembrane helices within a given protein.
	 * 
	 * @param protein
	 * @param cutoff
	 * @return
	 */
	private boolean splitTMHs(Protein protein, double cutoff)
	{
		boolean split 		= false;
		Pssm 	pssm 		= protein.getPssm();
		char[] 	structure 	= protein.getPrediction();
		int[] 	segmentRaw 	= protein.getSegmentRaw();
		int 	minLength 	= 2*Globals.PREDICTOR_HELIX_MIN_SIZE+Globals.PREDICTOR_GAP_MIN_SIZE;
		
		for (int i = 0; i < structure.length; ++i)
		{
			try
			{
				if (Mappings.ssToInt(structure[i]) == Mappings.indexTmh)
				{
					int 	start 	= i;
					char 	type 	= structure[i];
					
					//go to end of transmembrane helix
					while (i < structure.length && structure[i] == type) {++i;}
					
					--i;
					
					int end = i;
					
					//if TMH is too short jump to the next one
					if (end-start+1 < minLength) {continue;}
					
					Instance window = this.buildInstance(pssm, start, end);
					
					window.isMissing((Attribute)this.attributes.get(this.attributes.size()-1));
					window.setDataset(this.dataset);
					
					double[] probabilities = this.classifier.distributionForInstance(window);
					
					double 	bestProb 	= probabilities[Mappings.indexTmh];
					double 	bestProb1 	= 0;
					double 	bestProb2 	= 0;
					int 	bestBreak1 	= -1;
					int 	bestBreak2 	= -1;
					
					//insert a variable gap into the TMH and find best constellation
					for (int break1 = start+(Globals.PREDICTOR_HELIX_MIN_SIZE-1); break1 < end; ++break1)
					{
						for (int break2 = break1+Globals.PREDICTOR_GAP_MIN_SIZE; break2 < end-(Globals.PREDICTOR_HELIX_MIN_SIZE-1); ++break2)
						{
							if (break2 == break1) {continue;}
							
							Instance window1 = this.buildInstance(pssm, start, break1);
							Instance window2 = this.buildInstance(pssm, break2+1, end);
							
							window1.isMissing((Attribute)this.attributes.get(this.attributes.size()-1));
							window2.isMissing((Attribute)this.attributes.get(this.attributes.size()-1));
							window1.setDataset(this.dataset);
							window2.setDataset(this.dataset);
							
							double prob1 = this.classifier.distributionForInstance(window1)[Mappings.indexTmh];
							double prob2 = this.classifier.distributionForInstance(window2)[Mappings.indexTmh];
							
							if (prob1 >= cutoff && prob2 >= cutoff)
							{
								double avgProb = (prob1 + prob2) / 2.0;
								
								if (avgProb > bestProb)
								{
									bestProb 	= avgProb;
									bestProb1 	= prob1;
									bestProb2 	= prob2;
									bestBreak1 	= break1;
									bestBreak2 	= break2;
								}
							}
						}
					}
					
					//split the TMH!
					if (bestBreak1 != -1 && bestBreak2 != -1)
					{
						for (int j = start; j <= bestBreak1; ++j)
						{
							structure[j] 	= Mappings.intToSs(Mappings.indexTmh);
							segmentRaw[j] 	= (int)(1000 * bestProb1);
						}
						
						for (int j = bestBreak1+1; j <= bestBreak2; ++j)
						{
							structure[j] 	= Mappings.intToSs(Mappings.indexNotTmh);
							segmentRaw[j] 	= 0;
						}
						
						for (int j = bestBreak2+1; j <= end; ++j)
						{
							structure[j] 	= Mappings.intToSs(Mappings.indexTmh);
							segmentRaw[j] 	= (int)(1000 * bestProb2);
						}
						
						split = true;
					}
				}
			}
			catch (Exception e)
			{
				ErrorUtils.printError(HelixPredictor.class, "Prediction failed for " + protein.getHeader(), e);
				
				return false;
			}
		}
		
		return split;
	}
	
	
	/**
	 * Analyzes a given segment (TMH or not) and saves it in the database.
	 * 
	 * @param pssm
	 * @param start
	 * @param end
	 * @param structureIndex
	 */
	private void addSegmentToDatabse(Pssm pssm, int start, int end, int structureIndex)
	{
		Instance segment = this.buildInstance(pssm, start, end);
		
		segment.setValue((Attribute)this.attributes.get(this.attributes.size()-1), structureIndex);
		
		segment.setDataset(this.dataset);
		
		this.dataset.add(segment);
	}
	
	
	/**
	 * Analyzes a given segment and returns the TMH probability.
	 * 
	 * @param pssm
	 * @param start
	 * @param end
	 * @return
	 */
	public double getSegmentProbability(Pssm pssm, int start, int end)
	{
		double tmhProbability = -1;
		
		try
		{
			Instance window = this.buildInstance(pssm, start, end);
			
			window.isMissing((Attribute)this.attributes.get(this.attributes.size()-1));
			window.setDataset(this.dataset);
			
			tmhProbability = this.classifier.distributionForInstance(window)[Mappings.indexTmh];;
		}
		catch (Exception e)
		{
			ErrorUtils.printError(HelixPredictor.class, "Prediction failed for segment (" + start + "-" + end + ")", e);
			
			return -1.0;
		}
		
		return tmhProbability;
	}
	
	
	/**
	 * Converts a given segment (TMH or not) into a Weka Instance.
	 * 
	 * @param pssm
	 * @param start
	 * @param end
	 * @return
	 */
	private Instance buildInstance(Pssm pssm, int start, int end)
	{
		SparseInstance 	window 			= new SparseInstance(this.attributes.size());
		int 			length 			= end - start + 1;
		int 			attIndex 		= 0;
		int 			conserved 		= 0;
		int 			nonConserved 	= 0;
		double 			consAvgHydro 	= 0;
		double 			nonConsAvgHydro = 0;
		double 			consHydro 		= 0;
		double 			nonConsHydro 	= 0;
		double 			consCharged 	= 0;
		double 			nonConsCharged 	= 0;
		double[] 		consAaComp 		= new double[20];
		double[] 		nonConsAaComp 	= new double[20];
		
		//amino acid composition, hydrophobicity, and charge
		for (int i = start; i <= end; ++i)
		{
			for (int j = 0; j < 20; ++j)
			{
				int 	score 	= pssm.getScore(i, j);
				char 	aa 		= Mappings.intToAa(j);
				
				if (score > 0)
				{
					consAvgHydro += Mappings.hydrophobicity(aa);
					
					if (Mappings.hydrophobicity(aa) > 0) 	{++consHydro;}
					if (Mappings.charge(aa) != 0) 			{++consCharged;}
					
					++consAaComp[j];
					++conserved;
				}
				else if (score < 0)
				{
					nonConsAvgHydro += Mappings.hydrophobicity(aa);
					
					if (Mappings.hydrophobicity(aa) > 0) 	{++nonConsHydro;}
					if (Mappings.charge(aa) != 0) 			{++nonConsCharged;}
					
					++nonConsAaComp[j];
					++nonConserved;
				}
			}
		}
		
		conserved 		= Math.max(conserved, 1);
		nonConserved 	= Math.max(nonConserved, 1);
		
		for (int i = 0; i < consAaComp.length; ++i)
		{
			consAaComp[i] = consAaComp[i] / conserved;
			
			window.setValue((Attribute)this.attributes.get(attIndex++), consAaComp[i]);
		}
		
		for (int i = 0; i < nonConsAaComp.length; ++i)
		{
			nonConsAaComp[i] = nonConsAaComp[i] / nonConserved;
			
			window.setValue((Attribute)this.attributes.get(attIndex++), nonConsAaComp[i]);
		}
		
		window.setValue((Attribute)this.attributes.get(attIndex++), length);
		
		window.setValue((Attribute)this.attributes.get(attIndex++), consAvgHydro / conserved);
		window.setValue((Attribute)this.attributes.get(attIndex++), nonConsAvgHydro / nonConserved);
		
		window.setValue((Attribute)this.attributes.get(attIndex++), consHydro / conserved);
		window.setValue((Attribute)this.attributes.get(attIndex++), nonConsHydro / nonConserved);
		
		window.setValue((Attribute)this.attributes.get(attIndex++), consCharged / conserved);
		window.setValue((Attribute)this.attributes.get(attIndex++), nonConsCharged / nonConserved);
		
		return window;
	}
	
	
	/**
	 * Trains the Weka Classifer.
	 */
	public void trainClassifier()
	{
		try
		{
			MultilayerPerceptron 	classifier 	= new weka.classifiers.functions.MultilayerPerceptron();
			Instances 				data 		= this.dataset;
			
			if (data.classIndex() == -1)
			{
				data.setClassIndex(data.numAttributes() - 1);
			}
			
			data.randomize(new Random(data.size()));
			
			String[] optClassifier 	= weka.core.Utils.splitOptions("-L 0.01 -M 0.8 -N 256 -V 20 -S 0 -E 5 -H 25 -B -I -D -C");
			
			classifier.setOptions(optClassifier);
			classifier.setSeed(data.size());
			
			classifier.buildClassifier(data);
			
			this.classifier = classifier;
			this.isTrained 	= true;
		}
		catch (Exception e)
		{
			ErrorUtils.printError(HelixPredictor.class, "Training failed", e);
		}
	}
	
	
	/**
	 * Saves the Weka arff data and model into files.
	 * Extensions: Arff -> .arff; Model -> .[model|model.gz]
	 * 
	 * @param filename
	 */
	public void dumpModelToFile(String filename)
	{
		//if classifier is still untrained, train it!
		if (!this.isTrained)
		{
			this.trainClassifier();
		}
		
		try
		{
			ArffSaver saver = new ArffSaver();
			
			saver.setInstances(this.dataset);
			saver.setFile(new File(filename + ".arff"));
			saver.writeBatch();
			
			weka.core.SerializationHelper.write(filename + ".model", this.classifier);
			
			ModelHandler.saveGZip(this.classifier, filename + ".model.gz");
		}
		catch (Exception e)
		{
			ErrorUtils.printError(HelixPredictor.class, "Failed to write model file(s) " + filename + ".[arrf|model|model.gz]", e);
		}
	}
	
	
	/**
	 * Loads a previously saved model.
	 * 
	 * @param filename
	 */
	public void loadModelFromFile(String filename)
	{
		try
		{
			File model 		= new File(filename + ".model");
			File modelGzip 	= new File(filename + ".model.gz");
			
			if (model.exists())
			{
				this.classifier = (Classifier)weka.core.SerializationHelper.read(filename + ".model");
			}
			else if (modelGzip.exists())
			{
				this.classifier = ModelHandler.loadGZip(filename + ".model.gz");
			}
			else
			{
				throw new FileNotFoundException(filename + ".[model|model.gz]");
			}
		}
		catch (Exception e)
		{
			ErrorUtils.printError(HelixPredictor.class, "Failed to load model file(s) " + filename + ".[model|model.gz]", e);
		}
	}
}
