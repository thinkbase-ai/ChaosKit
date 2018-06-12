using System;
using System.Collections.Generic;

namespace ConceptStrings.ChaosTools
{
	/// <summary>
	/// Sugihara and May type neighborhood prediction
	/// </summary>
	[Serializable]
	internal class SugMay
	{
		internal SugMay(double [,] trainStim, double [,] trainTarg)
		{
			tree = new KDTree(trainStim);
			simplexVertices = trainStim.GetLength(1) + 1; // number of columns
			this.trainTarg = trainTarg;
		}
		internal double Predict(double [] validStim, int outputIndex, out double confidence)
		{
            List<KDTreePoint> list = new List<KDTreePoint>();
			tree.Search(validStim,simplexVertices, ref list);

			double result = 0.0f;
			double b = 0.0;
			double av = 0.0;
			double scale = 4.0;
			double averageRes = 0.0f;
			double averageDist = 0.0f;
			for(int n = 0; n < simplexVertices; n++)
			{
				KDTreePoint point = (KDTreePoint)list[n];
				averageRes += trainTarg[(int)point.index,outputIndex];
				averageDist += Math.Sqrt(point.distSquared);
			}
			if(averageDist != 0.0)
				av =  scale/averageDist;
			else av = scale;
			averageRes /= (double)simplexVertices;
			averageDist /= (double)simplexVertices;
			double stdDevRes = 0.0f;
			double stdDevDist = 0.0f;
			int nSignCount = 0;
			for(int n = 0; n < simplexVertices; n++)
			{
				KDTreePoint point = (KDTreePoint)list[n];
				double dist = Math.Sqrt(point.distSquared);
				b += Math.Exp(-dist * av);
				double dTemp = trainTarg[(int)point.index,outputIndex];    
				result += Math.Exp(-dist * av) * dTemp;
				double fTemp1 = dTemp - averageRes;
				double fTemp2 = dist - averageDist;
				stdDevRes += fTemp1 * fTemp1;
				stdDevDist += fTemp2 * fTemp2;
				if(dTemp >= 0.0) 
					nSignCount++; 
			}
			stdDevRes = Math.Sqrt(stdDevRes /(double)simplexVertices);
			stdDevDist = Math.Sqrt(stdDevDist /(double)simplexVertices);
			if(nSignCount < simplexVertices / 2) 
				nSignCount = simplexVertices - nSignCount;
			double fSignCount = (double)nSignCount / (double)simplexVertices;
			confidence = fSignCount * (1.0f - stdDevRes) * (1.0f - stdDevDist);
			if(b != 0.0)
				b = 1.0/b;
			else return 0.0;
			result *= (double)b;
			return result;
		}
		int simplexVertices;
		double [,] trainTarg;
		KDTree tree;
	}
}
