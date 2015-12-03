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

package io;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import util.ErrorUtils;
import weka.classifiers.Classifier;

public class ModelHandler {
	
	
	public static final int BUFFER_SIZE = 16 * 1024;
	
	
	public static void saveGZip(Classifier classifier, String filename)
	{
		FileOutputStream 		fos = null;
		BufferedOutputStream 	bos = null;
		GZIPOutputStream 		gos = null;
		ObjectOutputStream 		oos = null;
		
		try
		{
			fos = new FileOutputStream(filename);
			bos = new BufferedOutputStream(fos, BUFFER_SIZE);
			gos = new GZIPOutputStream(bos, BUFFER_SIZE);
			oos = new ObjectOutputStream(gos);
			
			oos.writeObject(classifier);
			oos.flush();
		}
		catch (IOException e)
		{
			ErrorUtils.printError(ModelHandler.class, "Failed to read " + filename, e);
		}
		finally
		{
			ErrorUtils.closeQuietly(oos);
			ErrorUtils.closeQuietly(gos);
			ErrorUtils.closeQuietly(bos);
			ErrorUtils.closeQuietly(fos);
		}
	}
	
	
	public static Classifier loadGZip(String filename)
	{
		Classifier 				cls = null;
		FileInputStream 		fis = null;
		BufferedInputStream 	bis = null;
		GZIPInputStream 		gis = null;
		ObjectInputStream 		ois = null;
		
		try
		{
			fis = new FileInputStream(filename);
			bis = new BufferedInputStream(fis, BUFFER_SIZE);
			gis = new GZIPInputStream(bis, BUFFER_SIZE);
			ois = new ObjectInputStream(gis);
			
			cls = (Classifier)ois.readObject();
		}
		catch (IOException e)
		{
			ErrorUtils.printError(ModelHandler.class, "Failed to read " + filename, e);
		}
		catch (ClassNotFoundException e)
		{
			ErrorUtils.printError(ModelHandler.class, "File is not a WEKA model: " + filename, e);
		}
		finally
		{
			ErrorUtils.closeQuietly(ois);
			ErrorUtils.closeQuietly(gis);
			ErrorUtils.closeQuietly(bis);
			ErrorUtils.closeQuietly(fis);
		}
		
		return cls;
	}

}
