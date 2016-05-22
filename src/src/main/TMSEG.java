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

package main;

import io.FWriter;

import java.io.File;
import java.io.FilenameFilter;
import java.net.MalformedURLException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.ArrayList;

import predictors.HelixIndexer;
import predictors.HelixPredictor;
import predictors.TopologyPredictor;
import processing.Processing;
import util.ErrorUtils;
import util.Globals;
import util.Mappings;
import data.FastaReader;
import data.Protein;
import data.Pssm;

public class TMSEG {
	
	
	private static HelixIndexer 		helixIndexer 		= null;
	private static HelixPredictor 		helixPredictor 		= null;
	private static TopologyPredictor 	topologyPredictor 	= null;
	
	private static boolean 	doMultiJob 	= false;
	private static boolean 	doAdjust 	= false;
	private static boolean 	doTopology 	= false;
	
	private static String 	fastaPath 	= null;
	private static String 	pssmPath 	= null;
	private static String 	outPath 	= null;
	private static String 	outPathRaw 	= null;
	private static String 	rootDir 	= null;
	
	
	public static void main(String[] args)
	{
		//parse parameter
		if (!parseParams(args))
		{
			printHelp();
			
			return;
		}
		
		//find TMSEG root directory
		if (!initRootDir(main.TMSEG.class))
		{
			ErrorUtils.printError(TMSEG.class, "Failed to initialize root directory", null);
			
			return;
		}
		
		//check all files
		if (!checkFiles()) {return;}
		
		//load the predictors
		loadPredictors();
		
		//run prediction(s)
		if (doMultiJob)
		{
			File 			fastaFolder = new File(fastaPath);
			FilenameFilter 	fastaFilter = new FilenameFilter()
			{
				@Override
				public boolean accept(File dir, String name)
				{
					return name.toLowerCase().endsWith(".fasta");
				}
			};
			
			for (File file : fastaFolder.listFiles(fastaFilter))
			{
				String fastaFile 	= file.getAbsolutePath();
				String fileName 	= new File(fastaFile).getName(); fileName = fileName.substring(0, fileName.length()-6);
				String pssmFile 	= new File(pssmPath + "/" + fileName + ".pssm").getAbsolutePath();
				String outFile 		= null;
				String outFileRaw 	= null;
				
				if (outPath != null)
				{
					outFile = new File(outPath + "/" + fileName + ".tmseg").getAbsolutePath();
				}
				
				if (outPathRaw != null)
				{
					outFileRaw = new File(outPathRaw + "/" + fileName + ".tmseg-raw").getAbsolutePath();
				}
				
				doPrediction(fastaFile, pssmFile, outFile, outFileRaw);
			}
		}
		else
		{
			doPrediction(fastaPath, pssmPath, outPath, outPathRaw);
		}
	}
	
	
	private static void doPrediction(String fastaFile, String pssmFile, String outFile, String outFileRaw)
	{
		Protein 			protein = null;
		ArrayList<Protein> 	tmpList = null;
		
		//read input FASTA file
		if (doAdjust)
		{
			tmpList = FastaReader.readStructureFile(fastaFile);
		}
		else
		{
			tmpList = FastaReader.readFastaFile(fastaFile);
		}
		
		if (tmpList != null && tmpList.size() > 0) {protein = tmpList.get(0);}
		
		//read intput PSSM file
		if (protein != null)
		{
			Pssm pssm = Pssm.newPssm(pssmFile, protein.getSequence().length);
			
			if (pssm != null)
			{
				protein.setPssm(pssm);
			}
			else
			{
				ErrorUtils.printError(TMSEG.class, "Failed to parse input PSSM file: " + pssmFile, null);
				
				return;
			}
		}
		else
		{
			ErrorUtils.printError(TMSEG.class, "Failed to parse input FASTA file: " + fastaFile, null);
			
			return;
		}
		
		//do standard prediction
		if (!doAdjust)
		{
			doIndexing(protein);
			doAdjustments(protein);
			doTopology(protein);
			
			Processing.assignConfidence(protein);
		}
		//do post-processing only
		else
		{
			protein.setPrediction(protein.getStructure());
			
			Processing.tmpCheck(protein);
			
			if (!doTopology) {doAdjustments(protein);}
			
			doTopology(protein);
		}
		
		//write prediction files
		if (outFile != null)
		{
			writeRefinedPrediction(protein, outFile);
		}
		
		if (outFileRaw != null)
		{
			writeRawPrediction(protein, outFileRaw);
		}
	}
	
	
	private static void doIndexing(Protein protein)
	{
		helixIndexer.predict(protein);
		
		Processing.process(protein, false, 7, 185, 60, 0);
	}
	
	
	private static void doAdjustments(Protein protein)
	{
		Globals.PREDICTOR_HELIX_MIN_SIZE 	= 17;
		Globals.PREDICTOR_MAX_SHIFT 		= 3;
		Globals.PREDICTOR_GAP_MIN_SIZE 		= 1;
		
		helixPredictor.predict(protein, 0.0);
		
		Processing.tmpCheck(protein);
	}
	
	
	private static void doTopology(Protein protein)
	{
		topologyPredictor.predict(protein, 0.45);
	}
	
	
	private static void writeRefinedPrediction(Protein protein, String outFile)
	{
		int[] 	confidence 	= protein.getConfidence();
		char[] 	sequence 	= protein.getSequence();
		char[] 	prediction 	= protein.getPrediction();
		
		if (sequence == null || prediction == null) {return;}
		
		StringBuilder segments = new StringBuilder();
		
		segments.append("# SEGMENT\tSTART\tEND\tRI\n");
		segments.append("##\n");
		
		for (int i = 0; i < prediction.length; ++i)
		{
			int 	start 	= i+1;
			char 	type 	= prediction[i];
			
			while (i < prediction.length && prediction[i] == type) {++i;}
			
			--i;
			
			int end = i+1;
			int ss 	= Mappings.ssToInt(type);
			int top = Mappings.topToInt(type);
			
			if (ss == Mappings.indexTmh)
			{
				int ri = -1; //RI=-1 if no prediction was made (i.e. adjustment-mode only)
				
				if (confidence != null) {ri = confidence[start];}
				
				segments.append("# TRANSMEM\t" + start + "\t" + end + "\t" + ri + "\n");
			}
			else if (ss == Mappings.indexLoop)
			{
				segments.append("# REENTRANT\t" + start + "\t" + end + "\n");
			}
			else if (ss == Mappings.indexSignal)
			{
				segments.append("# SIGNAL\t" + start + "\t" + end + "\n");
			}
			else if (ss == Mappings.indexNotTmh)
			{
				if (top == Mappings.indexInside)
				{
					segments.append("# INSIDE\t" + start + "\t" + end + "\n");
				}
				else if (top == Mappings.indexOutside)
				{
					segments.append("# OUTSIDE\t" + start + "\t" + end + "\n");
				}
				else
				{
					segments.append("# NON-MEM\t" + start + "\t" + end + "\n");
				}
			}
			else
			{
				segments.append("# UNKNOWN\t" + start + "\t" + end + "\n");
			}
		}
		
		segments.append("##");
		
		FWriter.openFile(outFile);
		FWriter.writeLine(segments.toString());
		FWriter.writeLine(protein.getHeader());
		FWriter.writeLine(new String(sequence));
		FWriter.writeLine(new String(prediction));
		FWriter.closeFile();
	}
	
	
	private static void writeRawPrediction(Protein protein, String outFileRaw)
	{
		int 			topRaw 		= protein.getTopologyRaw();
		int[] 			solRaw 		= protein.getSolRaw();
		int[] 			tmhRaw 		= protein.getTmhRaw();
		int[] 			sigRaw 		= protein.getSigRaw();
		int[] 			segRaw 		= protein.getSegmentRaw();
		char[] 			sequence 	= protein.getSequence();
		char[] 			prediction 	= protein.getPrediction();
		
		if (sequence == null || prediction == null) {return;}
		
		StringBuilder output = new StringBuilder();
		
		output.append("# " + protein.getHeader() + "\n");
		output.append("# SEQ\tSOL\tTMH\tSIG\tSEG\tTOP\tPRED\n");
		
		for (int i = 0; i < sequence.length; ++i)
		{
			output.append(sequence[i]);
			output.append("\t");
			
			if (solRaw != null && solRaw[i] >= 0)
			{
				output.append(solRaw[i]/1000.0);
				output.append("\t");
			}
			else
			{
				output.append('.');
				output.append("\t");
			}
			
			if (tmhRaw != null && tmhRaw[i] >= 0)
			{
				output.append(tmhRaw[i]/1000.0);
				output.append("\t");
			}
			else
			{
				output.append('.');
				output.append("\t");
			}
			
			if (sigRaw != null && sigRaw[i] >= 0)
			{
				output.append(sigRaw[i]/1000.0);
				output.append("\t");
			}
			else
			{
				output.append('.');
				output.append("\t");
			}
			
			if (segRaw != null && segRaw[i] >= 0)
			{
				output.append(segRaw[i]/1000.0);
				output.append("\t");
			}
			else
			{
				output.append('.');
				output.append("\t");
			}
			
			if (i == 0 && topRaw >= 0)
			{
				output.append(topRaw/1000.0);
				output.append("\t");
			}
			else
			{
				output.append('.');
				output.append("\t");
			}
			
			output.append(prediction[i]);
			output.append("\n");
		}
		
		FWriter.openFile(outFileRaw);
		FWriter.write(output.toString());
		FWriter.closeFile();
	}
	
	
	private static void loadPredictors()
	{
		helixIndexer 		= new HelixIndexer();
		helixPredictor 		= new HelixPredictor();
		topologyPredictor 	= new TopologyPredictor();
		
		helixIndexer.initialize();
		helixPredictor.initialize();
		topologyPredictor.initialize();
		
		helixIndexer.loadModelFromFile(rootDir + "/models/hIndexer");
		helixPredictor.loadModelFromFile(rootDir + "/models/hPredictor");
		topologyPredictor.loadModelFromFile(rootDir + "/models/tPredictor");
	}
	
	
	private static boolean parseParams(String[] args)
	{
		if (args == null || args.length < 1) {return false;}
		
		int maxIndex = args.length - 1;
		
		for (int i = 0; i <= maxIndex; ++i)
		{
			String param = args[i].trim();
			String value = null;
			
			if (param.equalsIgnoreCase("-i"))
			{
				if (i == maxIndex) {break;}
				
				value = args[i+1].trim();
				
				if (value.startsWith("-")) {continue;}
				
				fastaPath = new File(value).getAbsolutePath();
				
				++i;
			}
			else if (param.equalsIgnoreCase("-p"))
			{
				if (i == maxIndex) {break;}
				
				value = args[i+1].trim();
				
				if (value.startsWith("-")) {continue;}
				
				pssmPath = new File(value).getAbsolutePath();
				
				++i;
			}
			else if (param.equalsIgnoreCase("-o"))
			{
				if (i == maxIndex) {break;}
				
				value = args[i+1].trim();
				
				if (value.startsWith("-")) {continue;}
				
				outPath = new File(value).getAbsolutePath();
				
				++i;
			}
			else if (param.equalsIgnoreCase("-r"))
			{
				if (i == maxIndex) {break;}
				
				value = args[i+1].trim();
				
				if (value.startsWith("-")) {continue;}
				
				outPathRaw = new File(value).getAbsolutePath();
				
				++i;
			}
			else if (param.equalsIgnoreCase("-m"))
			{
				doMultiJob = true;
			}
			else if (param.equalsIgnoreCase("-x"))
			{
				doAdjust = true;
			}
			else if (param.equalsIgnoreCase("-t"))
			{
				doTopology = true;
			}
			else
			{
				ErrorUtils.printError(TMSEG.class, "Unkown parameter: " + param, null);
				
				return false;
			}
		}
		
		return checkParams();
	}
	
	
	private static boolean checkParams()
	{
		boolean passed = true;
		
		if (fastaPath == null)
		{
			ErrorUtils.printError(TMSEG.class, "Missing input FASTA file/folder", null);
			
			passed = false;
		}
		else if (doMultiJob && !(new File(fastaPath).isDirectory()))
		{
			ErrorUtils.printError(TMSEG.class, "Input FASTA path is not a folder", null);
			
			passed = false;
		}
		else if (!doMultiJob && (new File(fastaPath).isDirectory()))
		{
			ErrorUtils.printError(TMSEG.class, "Input FASTA path is not a file", null);
			
			passed = false;
		}
		
		if (pssmPath == null)
		{
			ErrorUtils.printError(TMSEG.class, "Missing input PSSM file/folder", null);
			
			passed = false;
		}
		else if (doMultiJob && !(new File(pssmPath).isDirectory()))
		{
			ErrorUtils.printError(TMSEG.class, "Input PSSM path is not a folder", null);
			
			passed = false;
		}
		else if (!doMultiJob && (new File(pssmPath).isDirectory()))
		{
			ErrorUtils.printError(TMSEG.class, "Input PSSM path is not a file", null);
			
			passed = false;
		}
		
		if (outPath == null && outPathRaw == null)
		{
			ErrorUtils.printError(TMSEG.class, "Missing output file/folder", null);
			
			passed = false;
		}
		else if (doMultiJob)
		{
			if (outPath != null && !(new File(outPath).isDirectory()))
			{
				ErrorUtils.printError(TMSEG.class, "Output path is not a folder", null);
				
				passed = false;
			}
			
			if (outPathRaw != null && !(new File(outPathRaw).isDirectory()))
			{
				ErrorUtils.printError(TMSEG.class, "Raw output path is not a folder", null);
				
				passed = false;
			}
		}
		else if (!doMultiJob)
		{
			if (outPath != null && (new File(outPath).isDirectory()))
			{
				ErrorUtils.printError(TMSEG.class, "Output path is not a file", null);
				
				passed = false;
			}
			
			if (outPathRaw != null && (new File(outPathRaw).isDirectory()))
			{
				ErrorUtils.printError(TMSEG.class, "Raw output path is not a file", null);
				
				passed = false;
			}
		}
		
		if (doTopology == true && doAdjust != true)
		{
			ErrorUtils.printError(TMSEG.class, "-t set, but -x is not", null);
			
			passed = false;
		}
		
		return passed;
	}
	
	
	private static boolean checkFiles()
	{
		File fasta 	= new File(fastaPath);
		File pssm 	= new File(pssmPath);
		File m1 	= new File(rootDir + "/models/hIndexer.model");
		File m1Gz 	= new File(rootDir + "/models/hIndexer.model.gz");
		File m2 	= new File(rootDir + "/models/hPredictor.model");
		File m2Gz 	= new File(rootDir + "/models/hPredictor.model.gz");
		File m3 	= new File(rootDir + "/models/tPredictor.model");
		File m3Gz 	= new File(rootDir + "/models/tPredictor.model.gz");
		
		if (!fasta.exists())
		{
			ErrorUtils.printError(TMSEG.class, "Could not find FASTA file/folder " + fastaPath, null);
			
			return false;
		}
		
		if (!pssm.exists())
		{
			ErrorUtils.printError(TMSEG.class, "Could not find PSSM file/folder " + pssmPath, null);
			
			return false;
		}
		
		if (!m1.exists() && !m1Gz.exists())
		{
			ErrorUtils.printError(TMSEG.class, "Could not find HelixIndexer model file", null);
			
			return false;
		}
		
		if (!m2.exists() && !m2Gz.exists())
		{
			ErrorUtils.printError(TMSEG.class, "Could not find HelixPredictor model file", null);
			
			return false;
		}
		
		if (!m3.exists() && !m3Gz.exists())
		{
			ErrorUtils.printError(TMSEG.class, "Could not find TopologyPredictor model file", null);
			
			return false;
		}
		
		return true;
	}
	
	
	@SuppressWarnings("rawtypes")
	private static boolean initRootDir(Class target)
	{
		URL 	url 	= null;
		String 	extURL 	= null;
		
		try
		{
			url = target.getProtectionDomain().getCodeSource().getLocation();
		}
		catch (SecurityException e)
		{
			url = target.getResource(target.getSimpleName() + ".class");
		}
		
		extURL = url.toExternalForm();
		
		if (extURL.endsWith(".jar"))
		{
			extURL = extURL.substring(0, extURL.lastIndexOf("/"));
		}
		else
		{
			String suffix = "/"+(target.getName()).replace(".", "/")+".class";
			
			extURL = extURL.replace(suffix, "");
			
			if (extURL.startsWith("jar:") && extURL.endsWith(".jar!"))
			{
				extURL = extURL.substring(4, extURL.lastIndexOf("/"));
			}
		}
		
		try
		{
			url = new URL(extURL);
		}
		catch (MalformedURLException e)
		{
			//do nothing
		}
		
		try
		{
			rootDir = new File(url.toURI()).getAbsolutePath();
		}
		catch (URISyntaxException e)
		{
			rootDir = new File(url.getPath()).getAbsolutePath();
		}
		
		if (rootDir != null && !rootDir.equalsIgnoreCase(""))
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	
	
	private static void printHelp()
	{
		System.out.println(	"TMSEG help. Please use the following parameters:\n" +
							"-i <path>      Input FASTA file/folder\n" +
							"-p <path>      Input PSSM file/folder\n" +
							"-o <path>      Output file/folder (human readable)\n" +
							"-r <path>      Output file/folder (raw prediction scores)\n" +
							"-m FLAG        if set, do multi-job (interpret input/output paths as folders)\n" +
							"-x FLAG        if set, a previous prediction is processed (must be supplied in FASTA file)\n" +
							"-t FLAG        if set, only the topology prediction is performed (-x must be set)");
	}

}
