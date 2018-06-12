using System;

namespace ConceptStrings.ChaosTools
{
	/// <summary>
	/// hyper point class
	/// </summary>
	internal class HyperPoint
	{
		internal double [] coord;

		internal HyperPoint(int n) 
		{
			coord = new double [n];
		}

		internal HyperPoint(double [] x) 
		{
			coord = new double[x.Length];
			for (int i=0; i<x.Length; ++i) coord[i] = x[i];
		}

		internal Object clone() 
		{
			return new HyperPoint(coord);
		}

		internal bool equals(HyperPoint p) 
		{
			for (int i=0; i<coord.Length; ++i)
				if (coord[i] != p.coord[i])
					return false;
			return true;
		}

		internal static double sqrdist(HyperPoint x, HyperPoint y) 
		{
			double dist = 0;
			for (int i=0; i<x.coord.Length; ++i) 
			{
				double diff = (x.coord[i] - y.coord[i]);
				dist += (diff*diff);
			}
			return dist;
		}

		internal static double eucdist(HyperPoint x, HyperPoint y) 
		{
			return Math.Sqrt(sqrdist(x, y));
		}
	}
}
