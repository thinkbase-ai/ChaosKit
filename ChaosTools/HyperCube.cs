using System;

namespace ConceptStrings.ChaosTools
{
	/// <summary>
	/// 
	/// </summary>
	[Serializable]
	internal class HyperCube
	{
		internal HyperCube(int wDimension, double [] pNearestCorner, double fEdge)
		{
			m_pNearestCorner = new double [wDimension];
			m_pVectorSum = new double [wDimension];
			m_wDimension = wDimension;
			for(int n = 0; n < wDimension; n++)
				m_pNearestCorner[n]  = pNearestCorner[n];
			m_wUseCount = 0;
			m_fEdge = fEdge;
		}

		internal bool AddVector(double [] pFirstPoint, double [] pSecondPoint)
		{
			return true;
		}
		internal double GetMagnitudeVectorSum()
		{
			if(m_wUseCount == 0)
				return 0.0f;
			double fTemp = 0.0f;
			for(int n = 0; n < m_wDimension; n++)								 
				fTemp += m_pVectorSum[n] * m_pVectorSum[n];
			return Math.Sqrt(fTemp);
		}
		internal int  GetUseCount(){return m_wUseCount;}

		int m_wDimension;
		double [] m_pNearestCorner;
//		double [] m_pFarthestCorner;
		double [] m_pVectorSum;
		int m_wUseCount;
		double m_fEdge;
	}
}
