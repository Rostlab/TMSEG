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

import java.io.Closeable;
import java.io.IOException;

public class ErrorUtils {
	
	
	public static void printError(Class<?> c, String m, Exception e)
	{
		if (c != null)
		{
			System.err.println("<" + c.getCanonicalName() + "> Error: " + m);
			
			if (e != null)
			{
				System.err.print("<" + c.getCanonicalName() + "> ...caused by: ");
				
				e.printStackTrace();
			}
		}
		else
		{
			System.err.println("Error: " + m);
			
			if (e != null)
			{
				System.err.print("...caused by: ");
				
				e.printStackTrace();
			}
		}
	}
	
	
	public static void printWarning(Class<?> c, String m, Exception e)
	{
		if (c != null)
		{
			System.err.println("<" + c.getCanonicalName() + "> Warning: " + m);
			
			if (e != null)
			{
				System.err.print("<" + c.getCanonicalName() + "> ...caused by: ");
				
				e.printStackTrace();
			}
		}
		else
		{
			System.err.println("Warning: " + m);
			
			if (e != null)
			{
				System.err.print("...caused by: ");
				
				e.printStackTrace();
			}
		}
	}
	
	
	public static void closeQuietly(Closeable c)
	{
		if (c != null)
		{
			try
			{
				c.close();
			}
			catch (IOException e) {}
		}
	}

}
