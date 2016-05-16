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

import processing.Processing;
import util.ErrorUtils;
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
import data.Segment;


/**
 * Class to predict the N-terminal topology of transmembrane proteins.
 */
public class TopologyPredictor {
	
	
	private static final int 		indexInside 	= 0;
	private static final int 		indexOutside 	= 1;
	
	private ArrayList<Attribute> 	attributes 		= null;
	private Instances 				dataset 		= null;
	private Classifier 				classifier 		= null;
	private boolean 				isTrained 		= false;
	
	
	public TopologyPredictor()
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
		
		for (int i = 0; i < 20; ++i)
		{
			this.attributes.add(new Attribute("conserved_composition_side1_" + Mappings.intToAa(i)));
		}
		
		for (int i = 0; i < 20; ++i)
		{
			this.attributes.add(new Attribute("non-conserved_composition_side1_" + Mappings.intToAa(i)));
		}
		
		for (int i = 0; i < 20; ++i)
		{
			this.attributes.add(new Attribute("conserved_composition_side2_" + Mappings.intToAa(i)));
		}
		
		for (int i = 0; i < 20; ++i)
		{
			this.attributes.add(new Attribute("non-conserved_composition_side2_" + Mappings.intToAa(i)));
		}
		
		this.attributes.add(new Attribute("conserved_pos_charged_side1"));
		this.attributes.add(new Attribute("non-conserved_pos_charged_side1"));
		
		this.attributes.add(new Attribute("conserved_pos_charged_side2"));
		this.attributes.add(new Attribute("non-conserved_pos_charged_side2"));
		
		this.attributes.add(new Attribute("conserved_pos_charged_diff"));
		this.attributes.add(new Attribute("non-conserved_pos_charged_diff"));
		
		ArrayList<String> classes = new ArrayList<String>();
		
		classes.add(String.valueOf(TopologyPredictor.indexInside));
		classes.add(String.valueOf(TopologyPredictor.indexOutside));
		
		this.attributes.add(new Attribute("class", classes));
	}
	
	
	/**
	 * Initializes the classifier and dataset.
	 */
	public void initialize()
	{
		this.isTrained 	= false;
		this.classifier = null;
		this.dataset 	= new Instances("TopologyPredictor Model", this.attributes, 0);
		
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
		if (protein.getPssm() == null) 		{return;}
		if (protein.getStructure() == null) {return;}
		
		Pssm 	pssm 		= protein.getPssm();
		char[] 	structure 	= protein.getStructure();
		
		if (Processing.testTopology(structure))
		{
			structure = Processing.extrapolateTopology(structure);
			
			int startPos 		= 0;
			int structureIndex 	= -1;
			
			for (int i = 0; i < structure.length; ++i)
			{
				int topIndex = Mappings.topToInt(structure[i]);
				
				if (topIndex != Mappings.defaultValue)
				{
					if (topIndex == Mappings.indexInside)
					{
						structureIndex = TopologyPredictor.indexInside;
					}
					else
					{
						structureIndex = TopologyPredictor.indexOutside;
					}
					
					startPos = i;
					
					break;
				}
			}
			
			if (structureIndex != -1)
			{
				this.addProteinToDatabse(pssm, structure, structureIndex, startPos);
			}
		}
	}
	
	
	/**
	 * Predicts the N-terminal topology for a given list of proteins.
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
	 * Predicts the N-terminal topology for a given protein.
	 * 
	 * @param protein
	 * @param cutoff
	 */
	public void predict(Protein protein, double cutoff)
	{
		if (protein == null) 					{return;}
		if (protein.getPssm() == null) 			{return;}
		if (protein.getPrediction() == null) 	{return;}
		
		if (!protein.isPredTmp()) {return;}
		
		Pssm 	pssm 		= protein.getPssm();
		char[] 	prediction 	= protein.getPrediction();
		
		try
		{
			ArrayList<Segment> 	solSegments = findSegments(prediction);
			Instance 			instance 	= this.buildInstance(pssm, prediction, solSegments, 0);
			
			instance.isMissing((Attribute)this.attributes.get(this.attributes.size()-1));
			instance.setDataset(this.dataset);
			
			double[] probabilities = this.classifier.distributionForInstance(instance);
			
			char top = Character.UNASSIGNED;
			
			protein.setTopologyRaw((int)(1000 * probabilities[TopologyPredictor.indexInside]));
			
			if (!protein.hasPredSigP() && probabilities[TopologyPredictor.indexInside] >= cutoff)
			{
				top = Mappings.intToTop(Mappings.indexInside);
			}
			else
			{
				top = Mappings.intToTop(Mappings.indexOutside);
			}
			
			for (int i = 0; i < prediction.length; ++i)
			{
				char type = prediction[i];
				
				if (Mappings.ssToInt(type) == Mappings.indexNotTmh)
				{
					prediction[i] = top;
				}
				else if (Mappings.ssToInt(type) == Mappings.indexTmh)
				{
					if (top == Mappings.intToTop(Mappings.indexInside))
					{
						top = Mappings.intToTop(Mappings.indexOutside);
					}
					else
					{
						top = Mappings.intToTop(Mappings.indexInside);
					}
					
					while (i < prediction.length && type == prediction[i]) {++i;}
					
					--i;
				}
			}
		}
		catch (Exception e)
		{
			ErrorUtils.printError(TopologyPredictor.class, "Prediction failed for " + protein.getHeader(), e);
			
			return;
		}
	}
	
	
	/**
	 * Analyzes a given window and saves it in the database.
	 * 
	 * @param pssm
	 * @param structure
	 * @param structureIndex
	 * @param startPos
	 */
	private void addProteinToDatabse(Pssm pssm, char[] structure, int structureIndex, int startPos)
	{
		ArrayList<Segment> 	solSegments = findSegments(structure);
		Instance 			segment 	= this.buildInstance(pssm, structure, solSegments, startPos);
		
		segment.setValue((Attribute)this.attributes.get(this.attributes.size()-1), structureIndex);
		
		segment.setDataset(this.dataset);
		
		this.dataset.add(segment);
	}
	
	
	/**
	 * Finds all transmembrane helices within a protein and generates
	 * the list of segments that will be used for the features of each
	 * side of the membrane.
	 * 
	 * @param structure
	 * @return
	 */
	private ArrayList<Segment> findSegments(char[] structure)
	{
		ArrayList<Segment> tmhSegments = new ArrayList<Segment>();
		ArrayList<Segment> solSegments = new ArrayList<Segment>();
		
		for (int i = 0; i < structure.length; ++i)
		{
			char type = structure[i];
			
			if (Mappings.ssToInt(type) == Mappings.indexTmh)
			{
				int start = i;
				
				while (i < structure.length && structure[i] == type) {++i;}
				
				--i;
				
				int end = i;
				
				tmhSegments.add(new Segment(start, end));
			}
		}
		
		int side = 0;
		int off1 = 8;
		int off2 = 15;
		
		for (int i = 0; i < tmhSegments.size(); ++i)
		{
			Segment segment = tmhSegments.get(i);
			
			if (i == 0)
			{
				int start 	= Math.max(segment.start - off2, 0);
				int end 	= Math.min(segment.start + off1, structure.length - 1);
				
				end 		= Math.min(segment.end, end); //do not extend beyond helix end
				
				solSegments.add(new Segment(start, end, side, 0));
			}
			else
			{
				Segment prev = tmhSegments.get(i - 1);
				
				int start 	= Math.max(prev.end - off1, 0);
				int end 	= Math.min(prev.end + off2, structure.length - 1);
				
				start 		= Math.max(prev.start, start); // do not extend beyond helix start
				
				end 		= Math.min(segment.start + off1, end); // do not extend beyond next helix cutoff
				end 		= Math.min(segment.end, end); // do not extend beyond next helix end
				
				solSegments.add(new Segment(start, end, side, 0));
				
				start 		= Math.max(segment.start - off2, 0);
				start 		= Math.max(prev.end - off1, start); // do not extend beyond prev helix cutoff
				start 		= Math.max(prev.start, start); // do not extend beyond prev helix start
				
				end 		= Math.min(segment.start + off1, structure.length - 1);
				end 		= Math.min(segment.end, end); // do not extend beyond helix end
				
				solSegments.add(new Segment(start, end, side, 0));
			}
			
			side = (side + 1) % 2;
			
			if (i == tmhSegments.size() - 1)
			{
				int start 	= Math.max(segment.end - off1, 0);
				int end 	= Math.min(segment.end + off2, structure.length - 1);
				
				start 		= Math.max(segment.start, start); //do not extend beyond helix start
				
				solSegments.add(new Segment(start, end, side, 0));
			}
		}
		
		return solSegments;
	}
	
	
	/**
	 * Converts a list of segments for both sides of the membrane into a WEKA instance.
	 * 
	 * @param pssm
	 * @param structure
	 * @param segments
	 * @param startPos
	 * @return
	 */
	private Instance buildInstance(Pssm pssm, char[] structure, ArrayList<Segment> segments, int startPos)
	{
		SparseInstance 	protein 			= new SparseInstance(this.attributes.size());
		double[] 		consComposition1 	= new double[20];
		double[] 		nonConsComposition1 = new double[20];
		double[] 		consComposition2 	= new double[20];
		double[] 		nonConsComposition2 = new double[20];
		int 			attIndex 			= 0;
		int 			conserved1 			= 0;
		int 			nonConserved1 		= 0;
		int 			conserved2 			= 0;
		int 			nonConserved2 		= 0;
		int 			consPosCharged1 	= 0;
		int 			nonConsPosCharged1 	= 0;
		int 			conPosCharged2 		= 0;
		int 			nonConsPosCharged2 	= 0;
		int 			firstSide 			= -1;
		
		
		for (Segment segment : segments)
		{
			if (segment.end < startPos) {continue;}
			if (firstSide == -1) 		{firstSide = segment.type;}
			
			int start 	= Math.max(segment.start, startPos);
			int end 	= segment.end;
			int side 	= segment.type;
			
			for (int i = start; i <= end; ++i)
			{
				if (Mappings.ssToInt(structure[i]) != Mappings.indexUnknown)
				{
					for (int j = 0; j < 20; ++j)
					{
						int 	score 	= pssm.getScore(i, j);
						char 	aa 		= Mappings.intToAa(j);
						
						if (score > 0)
						{
							if (side == firstSide)
							{
								++consComposition1[j];
								++conserved1;
								
								if (Mappings.charge(aa) > 0) {++consPosCharged1;}
							}
							else
							{
								++consComposition2[j];
								++conserved2;
								
								if (Mappings.charge(aa) > 0) {++conPosCharged2;}
							}
						}
						else if (score < 0)
						{
							if (side == firstSide)
							{
								++nonConsComposition1[j];
								++nonConserved1;
								
								if (Mappings.charge(aa) > 0) {++nonConsPosCharged1;}
							}
							else
							{
								++nonConsComposition2[j];
								++nonConserved2;
								
								if (Mappings.charge(aa) > 0) {++nonConsPosCharged2;}
							}
						}
					}
				}
			}
		}
		
		conserved1 		= Math.max(conserved1, 1);
		nonConserved1 	= Math.max(nonConserved1, 1);
		conserved2 		= Math.max(conserved2, 1);
		nonConserved2 	= Math.max(nonConserved2, 1);
		
		//normalize for length
		for (int i = 0; i < consComposition1.length; ++i)
		{
			consComposition1[i] = consComposition1[i] / conserved1;
			
			protein.setValue((Attribute)this.attributes.get(attIndex++), consComposition1[i]);
		}
		
		for (int i = 0; i < nonConsComposition1.length; ++i)
		{
			nonConsComposition1[i] = nonConsComposition1[i] / nonConserved1;
			
			protein.setValue((Attribute)this.attributes.get(attIndex++), nonConsComposition1[i]);
		}
		
		for (int i = 0; i < consComposition2.length; ++i)
		{
			consComposition2[i] = consComposition2[i] / conserved2;
			
			protein.setValue((Attribute)this.attributes.get(attIndex++), consComposition2[i]);
		}
		
		for (int i = 0; i < nonConsComposition2.length; ++i)
		{
			nonConsComposition2[i] = nonConsComposition2[i] / nonConserved2;
			
			protein.setValue((Attribute)this.attributes.get(attIndex++), nonConsComposition2[i]);
		}
		
		protein.setValue((Attribute)this.attributes.get(attIndex++), (double)consPosCharged1 / (double)conserved1);
		protein.setValue((Attribute)this.attributes.get(attIndex++), (double)nonConsPosCharged1 / (double)nonConserved1);
		
		protein.setValue((Attribute)this.attributes.get(attIndex++), (double)conPosCharged2 / (double)conserved2);
		protein.setValue((Attribute)this.attributes.get(attIndex++), (double)nonConsPosCharged2 / (double)nonConserved2);
		
		protein.setValue((Attribute)this.attributes.get(attIndex++), (consPosCharged1 - conPosCharged2));
		protein.setValue((Attribute)this.attributes.get(attIndex++), (nonConsPosCharged1 - nonConsPosCharged2));
		
		return protein;
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
			
			String[] optClassifier 	= weka.core.Utils.splitOptions("-I 100 -K 7 -S 1 -num-slots 1");
			
			classifier.setOptions(optClassifier);
			classifier.setSeed(data.size());
			
			classifier.buildClassifier(data);
			
			this.classifier = classifier;
			this.isTrained 	= true;
		}
		catch (Exception e)
		{
			ErrorUtils.printError(TopologyPredictor.class, "Training failed", e);
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
			ErrorUtils.printError(TopologyPredictor.class, "Failed to write model file(s) " + filename + ".[arrf|model|model.gz]", e);
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
			ErrorUtils.printError(TopologyPredictor.class, "Failed to load model file(s) " + filename + ".[model|model.gz]", e);
		}
	}
}
