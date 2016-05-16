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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

public class FWriter {
	
	
	private static File 			file 	= null;
	private static BufferedWriter 	writer 	= null;
	
	
	public static void openFile(String filename)
	{
		try
		{
			file = new File(filename);
			
			if (file.exists())
			{
				file.delete();
			}
			else
			{
				File parent = file.getParentFile();
				
				if (parent != null && !parent.exists())
				{
					parent.mkdirs();
				}
			}
			
			writer = new BufferedWriter(new FileWriter(file));
		}
		catch (Exception e)
		{
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	
	public static void openFileAppend(String filename)
	{
		try
		{
			file = new File(filename);
			
			if (!file.exists())
			{
				file.createNewFile();
			}
			
			writer = new BufferedWriter(new FileWriter(file, true));
		}
		catch (Exception e)
		{
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	
	public static void closeFile()
	{
		isInitialized();
		
		try
		{
			writer.flush();
			writer.close();
		}
		catch (Exception e)
		{
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	
	public static void writeLine(String line)
	{
		isInitialized();
		
		try
		{
			writer.write(line+"\n");
			writer.flush();
		}
		catch (Exception e)
		{
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	
	public static void write(String line)
	{
		isInitialized();
		
		try
		{
			writer.write(line);
			writer.flush();
		}
		catch (Exception e)
		{
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	
	public static void delete()
	{
		isInitialized();
		
		try
		{
			if (writer != null)
			{
				writer.close();
			}
			
			if (file != null && file.exists())
			{
				file.delete();
			}
		}
		catch (Exception e)
		{
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	
	public static void delete(String filename)
	{
		try
		{
			File tmpFile = new File(filename);
			
			if (tmpFile != null && tmpFile.exists())
			{
				tmpFile.delete();
			}
		}
		catch (Exception e)
		{
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	
	private static boolean isInitialized()
	{
		if (file != null && writer != null)
		{
			return true;
		}
		else
		{
			System.err.println("Error: FileWriter not yet initialized!");
			System.exit(1);
		}
		
		return false;
	}
}
