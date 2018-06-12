using System;

namespace ConceptStrings.ChaosTools
{
	/// <summary>
	/// hyper rectangle class 
	/// </summary>
	internal class HyperRect
	{
		internal HyperPoint min;
		internal HyperPoint max;

		internal HyperRect(int ndims) 
		{
			min = new HyperPoint(ndims);
			max = new HyperPoint(ndims);
		}

		internal HyperRect(HyperPoint vmin, HyperPoint vmax) 
		{

			min = (HyperPoint)vmin.clone();
			max = (HyperPoint)vmax.clone();
		}

		internal Object clone() 
		{
	
			return new HyperRect(min, max);
		}

		// from Moore's eqn. 6.6
		internal HyperPoint closest(HyperPoint t) 
		{

			HyperPoint p = new HyperPoint(t.coord.Length);

			for (int i=0; i<t.coord.Length; ++i) 
			{
				if (t.coord[i]<=min.coord[i]) 
				{
					p.coord[i] = min.coord[i];
				}
				else if (t.coord[i]>=max.coord[i]) 
				{
					p.coord[i] = max.coord[i];
				}
				else 
				{
					p.coord[i] = t.coord[i];
				}
			}
	
			return p;
		}

		// used in initial conditions of KDTree.nearest()
		internal static HyperRect infiniteHyperRect(int d) 
		{
	
			HyperPoint vmin = new HyperPoint(d);
			HyperPoint vmax = new HyperPoint(d);
	
			for (int i=0; i<d; ++i) 
			{
				vmin.coord[i] = Double.NegativeInfinity;
				vmax.coord[i] = Double.PositiveInfinity;
			}

			return new HyperRect(vmin, vmax);
		}

		// currently unused
		internal HyperRect intersection(HyperRect r) 
		{

			HyperPoint newmin = new HyperPoint(min.coord.Length);
			HyperPoint newmax = new HyperPoint(min.coord.Length);

			for (int i=0; i<min.coord.Length; ++i) 
			{
				newmin.coord[i] = Math.Max(min.coord[i], r.min.coord[i]);
				newmax.coord[i] = Math.Min(max.coord[i], r.max.coord[i]);
				if (newmin.coord[i] >= newmax.coord[i]) return null;
			}

			return new HyperRect(newmin, newmax);
		}

		// currently unused
		internal double area () 
		{

			double a = 1;

			for (int i=0; i<min.coord.Length; ++i) 
			{
				a *= (max.coord[i] - min.coord[i]);
			}

			return a;
		}
	}
}
