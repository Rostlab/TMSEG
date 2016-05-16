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

public class Segment {
	
	
	public int start 	= -1;
	public int end 		= -1;
	public int type 	= 0;
	public int score 	= 0;
	public int shift 	= Integer.MAX_VALUE;
	
	
	public Segment(int start, int end)
	{
		if (start <= end)
		{
			this.start 	= start;
			this.end 	= end;
		}
	}
	
	
	public Segment(int start, int end, int type, int score)
	{
		if (start <= end)
		{
			this.start 	= start;
			this.end 	= end;
			this.type 	= type;
			this.score 	= score;
		}
	}

}
