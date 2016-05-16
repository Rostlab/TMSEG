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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

public class FReader {
	
	
	private static File 			file 	= null;
	private static BufferedReader 	reader 	= null;
	
	
	public static void openFile(String filename)
	{
		try
		{
			file 	= new File(filename);
			reader 	= new BufferedReader(new FileReader(file));
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
			reader.close();
		}
		catch (Exception e)
		{
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	
	public static String readLine()
	{
		isInitialized();
		
		try
		{
			if (reader.ready())
			{
				return reader.readLine();
			}
		}
		catch (Exception e)
		{
			e.printStackTrace();
			System.exit(1);
		}
		
		return null;
	}
	
	
	private static boolean isInitialized()
	{
		if (file != null && reader != null)
		{
			return true;
		}
		else
		{
			System.err.println("Error: FileReader not yet initialized!");
			System.exit(1);
		}
		
		return false;
	}
}
