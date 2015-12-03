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
import weka.classifiers.trees.RandomForest;
import weka.core.Attribute;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.SparseInstance;
import weka.core.converters.ArffSaver;
import data.Protein;
import data.Pssm;


/**
 * Class to predict transmembrane residues within a protein.
 */
public class HelixIndexer {
	
	
	public static final int 		indexNotTmh 	= 0;
	public static final int 		indexTmh 		= 1;
	public static final int 		indexSignal 	= 2;
	
	private ArrayList<Attribute> 	attributes 		= null;
	private Instances 				dataset 		= null;
	private Classifier 				classifier 		= null;
	private boolean 				isTrained 		= false;
	
	private double[] 				globalConsAa 	= null;
	private double[] 				globalNonConsAa = null;
	
	
	public HelixIndexer()
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
		
		for (int i = 1; i <= ((2*Globals.INDEXER_WINDOW_SIZE)+1); ++i)
		{
			for (int j = 0; j < 21; ++j)
			{
				this.attributes.add(new Attribute(Mappings.intToAa(j)+"_"+i));
			}
		}
		
		this.attributes.add(new Attribute("conserved_hydrophobicity"));
		this.attributes.add(new Attribute("non-conserved_hydrophobicity"));
		
		this.attributes.add(new Attribute("conserved_hydrophobic"));
		this.attributes.add(new Attribute("non-conserved_hydrophobic"));
		
		this.attributes.add(new Attribute("conserved_pos_charged"));
		this.attributes.add(new Attribute("non-conserved_pos_charged"));
		
		this.attributes.add(new Attribute("conserved_neg_charged"));
		this.attributes.add(new Attribute("non-conserved_neg_charged"));
		
		this.attributes.add(new Attribute("conserved_polar"));
		this.attributes.add(new Attribute("non-conserved_polar"));
		
		ArrayList<String> lengths = new ArrayList<String>();
		
		lengths.add(String.valueOf("0"));
		lengths.add(String.valueOf("1"));
		lengths.add(String.valueOf("2"));
		lengths.add(String.valueOf("3"));
		lengths.add(String.valueOf("4"));
		
		this.attributes.add(new Attribute("n-term_distance", lengths));
		this.attributes.add(new Attribute("c-term_distance", lengths));
		this.attributes.add(new Attribute("global_length", lengths));
		
		for (int j = 0; j < 20; ++j)
		{
			this.attributes.add(new Attribute("global_conserved_"+Mappings.intToAa(j)));
			this.attributes.add(new Attribute("global_non-conserved_"+Mappings.intToAa(j)));
		}
		
		ArrayList<String> classes = new ArrayList<String>();
		
		classes.add(String.valueOf(HelixIndexer.indexNotTmh));
		classes.add(String.valueOf(HelixIndexer.indexTmh));
		classes.add(String.valueOf(HelixIndexer.indexSignal));
		
		this.attributes.add(new Attribute("class", classes));
	}
	
	
	/**
	 * Initializes the classifier and dataset.
	 */
	public void initialize()
	{
		this.isTrained 	= false;
		this.classifier = null;
		this.dataset 	= new Instances("HelixIndexer Model", this.attributes, 0);
		
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
		int 	length 		= pssm.getLength();
		char[] 	structure 	= protein.getStructure();
		
		if (pssm.getLength() != structure.length)
		{
			ErrorUtils.printError(HelixIndexer.class, "PSSM and structure annotation length do not match for " + protein.getName(), null);
			
			return;
		}
		
		this.globalComposition(pssm);
		
		for (int i = 0; i < length; ++i)
		{
			if (Mappings.ssToInt(structure[i]) != Mappings.indexUnknown)
			{
				this.addWindowToDatabase(pssm, i, structure);
			}
		}
	}
	
	
	/**
	 * Predicts transmembrane residues for a given list of proteins.
	 * 
	 * @param proteins
	 */
	public void predict(List<Protein> proteins)
	{
		for (Protein protein : proteins)
		{
			this.predict(protein);
		}
	}
	
	
	/**
	 * Predicts transmembrane residues for a given protein.
	 * 
	 * @param protein
	 */
	public void predict(Protein protein)
	{
		if (protein == null || protein.getPssm() == null) {return;}
		
		Pssm 		pssm 		= protein.getPssm();
		int 		length 		= pssm.getLength();
		int[] 		scoresSol 	= new int[length];
		int[] 		scoresTmh 	= new int[length];
		int[] 		scoresSig 	= new int[length];
		
		this.globalComposition(pssm);
		
		//slide window along the sequence
		for (int i = 0; i < length; ++i)
		{
			try
			{
				Instance window = this.buildInstance(pssm, i);
				
				window.isMissing((Attribute)this.attributes.get(this.attributes.size()-1));
				window.setDataset(this.dataset);
				
				double[] probabilities = this.classifier.distributionForInstance(window);
				
				if (i < 40)
				{
					scoresSol[i] = (int)(1000 * probabilities[HelixIndexer.indexNotTmh]);
					scoresTmh[i] = (int)(1000 * probabilities[HelixIndexer.indexTmh]);
					scoresSig[i] = (int)(1000 * probabilities[HelixIndexer.indexSignal]);
				}
				else
				{
					scoresSol[i] = (int)(1000 * probabilities[HelixIndexer.indexNotTmh]);
					scoresTmh[i] = (int)(1000 * probabilities[HelixIndexer.indexTmh]);
					scoresSig[i] = 0;
				}
			}
			catch (Exception e)
			{
				ErrorUtils.printError(HelixIndexer.class, "Prediction failed for " + protein.getHeader(), e);
				
				return;
			}
		}
		
		//save scores into the protein
		protein.setSolRaw(scoresSol);
		protein.setTmhRaw(scoresTmh);
		protein.setSigRaw(scoresSig);
	}
	
	
	/**
	 * Analyzes a given window and saves it in the database.
	 * 
	 * @param pssm
	 * @param windowCenter
	 * @param structure
	 */
	private void addWindowToDatabase(Pssm pssm, int windowCenter, char[] structure)
	{
		int 		index 	= Mappings.ssToInt(structure[windowCenter]);
		Instance 	window 	= this.buildInstance(pssm, windowCenter);
		
		if 		(index == Mappings.indexTmh) 	{index = HelixIndexer.indexTmh;}
		else if (index == Mappings.indexSignal) {index = HelixIndexer.indexSignal;}
		else 									{index = HelixIndexer.indexNotTmh;}
		
		window.setValue((Attribute)this.attributes.get(this.attributes.size()-1), index);
		window.setDataset(this.dataset);
		
		this.dataset.add(window);
	}
	
	
	/**
	 * Calculates the amino acid composition of the whole protein.
	 * 
	 * @param pssm
	 */
	private void globalComposition(Pssm pssm)
	{
		int 		conserved 			= 0;
		int 		nonConserved 		= 0;
		double[] 	globalConserved 	= new double[20];
		double[] 	globalNonConserved 	= new double[20];
		
		for (int i = 0; i < pssm.getLength(); ++i)
		{
			for (int j = 0; j < 20; ++j)
			{
				int score = pssm.getScore(i, j);
				
				if (score > 0)
				{
					++globalConserved[j];
					++conserved;
				}
				else if (score < 0)
				{
					++globalNonConserved[j];
					++nonConserved;
				}
			}
		}
		
		for (int i = 0; i < 20; ++i)
		{
			globalConserved[i] 		= globalConserved[i] / conserved;
			globalNonConserved[i] 	= globalNonConserved[i] / nonConserved;
		}
		
		this.globalConsAa 		= globalConserved;
		this.globalNonConsAa 	= globalNonConserved;
	}
	
	
	/**
	 * Converts a given window into a Weka Instance.
	 * 
	 * @param pssm
	 * @param windowCenter
	 * @return
	 */
	private Instance buildInstance(Pssm pssm, int windowCenter)
	{
		SparseInstance 	window 			= new SparseInstance(this.attributes.size());
		int 			windowStart 	= windowCenter - Globals.INDEXER_WINDOW_SIZE;
		int 			windowStop 		= windowCenter + Globals.INDEXER_WINDOW_SIZE;
		int 			globalLenght 	= pssm.getLength();
		int 			nTermDist 		= windowCenter + 1;
		int 			cTermDist 		= globalLenght - windowCenter;
		int 			attIndex 		= 0;
		int 			conserved 		= 0;
		int 			nonConserved 	= 0;
		double 			consAvgHydro 	= 0;
		double 			nonConsAvgHydro = 0;
		double 			consHydro 		= 0;
		double 			nonConsHydro 	= 0;
		double 			consPCharged 	= 0;
		double 			nonConsPCharged = 0;
		double 			consNCharged 	= 0;
		double 			nonConsNCharged = 0;
		double 			consPolar 		= 0;
		double 			nonConsPolar 	= 0;
		
		//amino acid at position i in window
		for (int i = windowStart; i <= windowStop; ++i)
		{
			if (i >= 0 && i < globalLenght)
			{
				for (int j = 0; j < 20; ++j)
				{
					int score = pssm.getScore(i, j);
					
					if (Math.abs(i - windowCenter) <= Globals.INDEXER_INNER_WINDOW_SIZE)
					{
						char aa = Mappings.intToAa(j);
						
						if (score > 0)
						{
							consAvgHydro += Mappings.hydrophobicity(aa);

							if (Mappings.hydrophobicity(aa) > 0) 	{++consHydro;}
							if (Mappings.charge(aa) > 0) 			{++consPCharged;}
							if (Mappings.charge(aa) < 0) 			{++consNCharged;}
							if (Mappings.polarity(aa) > 0) 			{++consPolar;}
							
							++conserved;
						}
						else if (score < 0)
						{
							nonConsAvgHydro += Mappings.hydrophobicity(aa);
							
							if (Mappings.hydrophobicity(aa) > 0) 	{++nonConsHydro;}
							if (Mappings.charge(aa) > 0) 			{++nonConsPCharged;}
							if (Mappings.charge(aa) < 0) 			{++nonConsNCharged;}
							if (Mappings.polarity(aa) > 0) 			{++nonConsPolar;}
							
							++nonConserved;
						}
					}
					
					window.setValue((Attribute)this.attributes.get(attIndex++), score);
				}
				
				window.setValue((Attribute)this.attributes.get(attIndex++), -10);
			}
			else
			{
				for (int j = 0; j < 20; ++j)
				{
					window.setValue((Attribute)this.attributes.get(attIndex++), 0);
				}
				
				window.setValue((Attribute)this.attributes.get(attIndex++), 10);
			}
		}
		
		conserved 		= Math.max(conserved, 1);
		nonConserved 	= Math.max(nonConserved, 1);
		
		window.setValue((Attribute)this.attributes.get(attIndex++), consAvgHydro / conserved);
		window.setValue((Attribute)this.attributes.get(attIndex++), nonConsAvgHydro / nonConserved);
		
		window.setValue((Attribute)this.attributes.get(attIndex++), consHydro / conserved);
		window.setValue((Attribute)this.attributes.get(attIndex++), nonConsHydro / nonConserved);
		
		window.setValue((Attribute)this.attributes.get(attIndex++), consPCharged / conserved);
		window.setValue((Attribute)this.attributes.get(attIndex++), nonConsPCharged / nonConserved);
		
		window.setValue((Attribute)this.attributes.get(attIndex++), consNCharged / conserved);
		window.setValue((Attribute)this.attributes.get(attIndex++), nonConsNCharged / nonConserved);
		
		window.setValue((Attribute)this.attributes.get(attIndex++), consPolar / conserved);
		window.setValue((Attribute)this.attributes.get(attIndex++), nonConsPolar / nonConserved);
		
		if 		(nTermDist > 40) 	{nTermDist = 4;}
		else if (nTermDist > 30) 	{nTermDist = 3;}
		else if (nTermDist > 20) 	{nTermDist = 2;}
		else if (nTermDist > 10) 	{nTermDist = 1;}
		else 						{nTermDist = 0;}
		
		if 		(cTermDist > 40) 	{cTermDist = 4;}
		else if (cTermDist > 30) 	{cTermDist = 3;}
		else if (cTermDist > 20) 	{cTermDist = 2;}
		else if (cTermDist > 10) 	{cTermDist = 1;}
		else 						{cTermDist = 0;}
		
		if 		(globalLenght > 240) 	{globalLenght = 4;}
		else if (globalLenght > 180) 	{globalLenght = 3;}
		else if (globalLenght > 120) 	{globalLenght = 2;}
		else if (globalLenght > 60) 	{globalLenght = 1;}
		else 							{globalLenght = 0;}
		
		window.setValue((Attribute)this.attributes.get(attIndex++), nTermDist);
		window.setValue((Attribute)this.attributes.get(attIndex++), cTermDist);
		window.setValue((Attribute)this.attributes.get(attIndex++), globalLenght);
		
		for (int i = 0; i < 20; ++i)
		{
			window.setValue((Attribute)this.attributes.get(attIndex++), this.globalConsAa[i]);
			window.setValue((Attribute)this.attributes.get(attIndex++), this.globalNonConsAa[i]);
		}
		
		return window;
	}
	
	
	/**
	 * Trains the Weka Classifer.
	 */
	public void trainClassifier()
	{
		try
		{
			RandomForest 	classifier 	= new weka.classifiers.trees.RandomForest();
			Instances 		data 		= this.dataset;
			
			if (data.classIndex() == -1)
			{
				data.setClassIndex(data.numAttributes() - 1);
			}
			
			data.randomize(new Random(data.size()));
			
			String[] optClassifier 	= weka.core.Utils.splitOptions("-I 100 -K 9 -S 1 -num-slots 3");
			
			classifier.setOptions(optClassifier);
			classifier.setSeed(data.size());
			
			classifier.buildClassifier(data);
			
			this.classifier = classifier;
			this.isTrained 	= true;
		}
		catch (Exception e)
		{
			ErrorUtils.printError(HelixIndexer.class, "Training failed", e);
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
			ErrorUtils.printError(HelixIndexer.class, "Failed to write model file(s) " + filename + ".[arrf|model|model.gz]", e);
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
			ErrorUtils.printError(HelixIndexer.class, "Failed to load model file(s) " + filename + ".[model|model.gz]", e);
		}
	}
}
