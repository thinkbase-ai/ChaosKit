using System;

namespace ConceptStrings.ChaosTools
{
	/// <summary>
	/// 
	/// </summary>
	[Serializable]
	internal class HyperCubeList
	{
		internal HyperCubeList(int wDimension, double dEdge)
		{
			dimension = wDimension;
			edge = dEdge;
		}

		internal bool AddCubes(double [] pFirstVector,double [] pSecondVector)
		{
			// we require that the vector formed by the difference of the two above is bigger than 
			// a hypercube.
			// First find the nearest cube.
			double [] pFirstCubeStart = new double [dimension];
			double [] pLastCubeStart = new double [dimension];
			double fLength = 0.0f;
			for(int n = 0; n < dimension; n++)
			{
				pFirstCubeStart[n] = Math.Floor((pFirstVector[n] / edge)) * edge;
				fLength += (pFirstVector[n] - pSecondVector[n]) * (pFirstVector[n] - pSecondVector[n]);
			}
			fLength = Math.Sqrt(fLength);
			if(fLength < edge)
			{
				return false;
			}
			for(int n = 0; n < dimension; n++)
			{
				pLastCubeStart[n] = Math.Floor((pSecondVector[n] / edge)) * edge;
			}
			// if pLastCubeStart and pFirstCubeStart are adjacent, we still don't have a int enough vector
			bool bAdjacent = true;
			for	(int n = 0; n < dimension; n++)
			{
				double fDist = pFirstCubeStart[n] - pLastCubeStart[n];
				if(fDist > edge || fDist < -edge)
					bAdjacent = false;
			}
			if(bAdjacent)
			{
				return false;
			}
			return false;
		}
		internal double EvaluateVectors()
		{
			return 0.0f;
		}

		int dimension;
		double edge;
	}
}
