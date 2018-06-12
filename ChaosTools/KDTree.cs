using System;
using System.Collections.Generic;

namespace ConceptStrings.ChaosTools
{
    /// <summary>
    /// Must be implemented by objects that can be sorted by KDTrees.
    /// </summary>
    /// <remarks>KDTrees can be generated for arrays of doubles or objects 
    /// possessing a reference quantity that can be converted into a double. 
    /// Implement this interface in any object you wish to store on KDTrees.</remarks>
    public interface ITreeReferenceable
	{
		/// <summary>
		/// get the reference value
		/// </summary>
		/// <returns>the reference value</returns>
		double GetValue();
	}
	/// <summary>
	/// Data structure used to efficiently find near-neighboring points for numeric vectors from a database of other points
	/// </summary>
	/// <remarks>Finds either the 'n' nearest vectors to a given vector in Euclidean space, or all the neighbors within a given Euclidean distance from a vector.
	/// Vectors can be of arbitrary dimension but must be numeric and consistent in dimensionality. Processing occurs in approximately O(logN) time.</remarks>
	[Serializable]
	public class KDTree
	{
		/// <summary>
		/// number of dimensions
		/// </summary>
		private int dimensions;

		/// <summary>
		/// root of KD-tree
		/// </summary>
		private KDNode root;

		/// <summary>
		/// count of nodes
		/// </summary>
		private int m_count;
		internal int duplicates = 0;

		/// <summary>
		/// Constructor
		/// </summary>
		private KDTree()
		{
			dimensions = 0;
			root = null;
		}
		/// <summary>
		/// Create a tree from parameters
		/// </summary>
		/// <param name="inputMatrix">Matrix of data points</param>
		public KDTree(double [,] inputMatrix)
		{    
			dimensions = inputMatrix.GetLength(1);
			for(int n = 0; n < inputMatrix.GetLength(0);n++)
			{
				double [] val = new double[inputMatrix.GetLength(1)];
				for(int p = 0; p < inputMatrix.GetLength(1);p++)
				{
					val[p] = inputMatrix[n,p];
				}
				this.Insert(val,n);
			}
		}
		/// <summary>
		/// Alternate constructor for creating trees of objects 
		/// implementing the <see cref="ITreeReferenceable"/> interface.
		/// </summary>
		/// <param name="inputMatrix">two dimensional array of objects</param>
		public KDTree(ITreeReferenceable [,] inputMatrix)
		{
			dimensions = inputMatrix.GetLength(1);
			for(int n = 0; n < inputMatrix.GetLength(0);n++)
			{
				double [] val = new double[inputMatrix.GetLength(1)];
				for(int p = 0; p < inputMatrix.GetLength(1);p++)
				{
					val[p] = inputMatrix[n,p].GetValue();
				}
				this.Insert(val,n);
			}
		}

		/// <summary>
		/// Insert a node in a KD-tree.  	
		/// </summary>
		/// <remarks>
		/// Key is a vector of doubles used to locate the value in a cartesian space.
		/// The value parameters should each be unique. Duplicate objects, or those that appear duplicates because the are indistinguishable will create exceptions. 
		/// value can be any object, but for efficiency should be a simple object for which efficient hashing algorithms are defined.
		/// To speed up deletions a reverse index is maintained of values in a hashtable.</remarks>
		/// <exception cref="ArgumentOutOfRangeException">If Key has wrong number of dimensions</exception>
		/// <param name="key">key for KD-tree node</param>
		/// <param name="value">External value associated with that key</param>
		public void Insert(double [] key, Object value)
		{
			if (key.Length != dimensions) 
			{
				throw new ArgumentOutOfRangeException();
			}

			else 
			{
				try
				{
					root = KDNode.ins(new HyperPoint(key), value, root, 0, dimensions);
				}
				catch
				{
					duplicates++;
				}
			}
			m_count++;
		}

		/// <summary>
		/// Look for 'n' neighbour objects to a given point. Neighbours may be duplicates.
		/// </summary>
		/// <param name="inputVector">vector of objects to find neighbours of.</param>
		/// <param name="neighbourCount">number of neighbours to find</param>
		/// <param name="distances">array of distances</param>
		public void Search( ITreeReferenceable [] inputVector, int neighbourCount, ref List<KDTreePoint> distances)  
		{
			double [] newInputVector = new double[inputVector.Length];
			for(int n = 0; n < inputVector.Length; n++)
				newInputVector[n] = inputVector[n].GetValue();
			Search(newInputVector,neighbourCount, ref distances);
		}

		/// <summary>
		/// Look for 'n' neighbour objects to a given point. Neighbours may be duplicates.
		/// </summary>
		/// <param name="inputVector">vector of doubles to find neighbours of.</param>
		/// <param name="neighbourCount">number of neighbours to find</param>
		/// <param name="distances">array of distances</param>
        public void Search(double[] inputVector, int neighbourCount, ref List<KDTreePoint> distances)
		{
			distances.Clear();
			if (neighbourCount < 0 || neighbourCount > m_count) 
			{
				throw new ArgumentOutOfRangeException("Number of neighbors cannot" +
					" be negative or greater than number of nodes");
			}

			if (inputVector.Length != dimensions) 
			{
				throw new ArgumentOutOfRangeException();
			}

			Object [] nbrs = new Object [neighbourCount];
			NearestNeighborList nnl = new NearestNeighborList(neighbourCount);

			// initial call is with infinite hyper-rectangle and max distance
			HyperRect hr = HyperRect.infiniteHyperRect(inputVector.Length);
			double max_dist_sqd = Double.MaxValue;
			HyperPoint keyp = new HyperPoint(inputVector);

			KDNode.nnbr(root, keyp, hr, max_dist_sqd, 0, dimensions, nnl);
			while(distances.Count < neighbourCount)
			{
				KDNode kd = (KDNode)nnl.removeHighest();
				double dist = HyperPoint.eucdist(keyp,kd.k);
				if(kd.duplicates != null)
				{
					kd.duplicates.Add(kd.v);
					foreach(object obj in kd.duplicates)
					{
						distances.Add(new KDTreePoint(obj,dist*dist));
						if(distances.Count >= neighbourCount)
							break;//found enough
					}
				}
				else
				{
					distances.Add(new KDTreePoint(kd.v,dist*dist));
				}
			}
			distances.Sort();
		}


		/// <summary>
		/// Look for 'n' neighbours to a given point. Each neighbour must be unique..
		/// </summary>
		/// <param name="inputVector">vector to find neighbours of</param>
		/// <param name="neighbourCount">number of neighbours to find</param>
		/// <param name="distances">array of distances</param>
        public void UniqueSearch(double[] inputVector, int neighbourCount, ref List<KDTreePoint> distances)  //   
		{
			distances.Clear();
			if (neighbourCount < 0 || neighbourCount > m_count) 
			{
				throw new ArgumentOutOfRangeException("Number of neighbors cannot" +
					" be negative or greater than number of nodes");
			}

			if (inputVector.Length != dimensions) 
			{
				throw new ArgumentOutOfRangeException();
			}

			Object [] nbrs = new Object [neighbourCount];
			NearestNeighborList nnl = new NearestNeighborList(neighbourCount);

			// initial call is with infinite hyper-rectangle and max distance
			HyperRect hr = HyperRect.infiniteHyperRect(inputVector.Length);
			double max_dist_sqd = Double.MaxValue;
			HyperPoint keyp = new HyperPoint(inputVector);

			KDNode.nnbr(root, keyp, hr, max_dist_sqd, 0, dimensions, nnl);
			while(distances.Count < neighbourCount)
			{
				KDNode kd = (KDNode)nnl.removeHighest();
				double dist = HyperPoint.eucdist(keyp,kd.k);
				distances.Add(new KDTreePoint(kd.v,dist*dist));
			}
			distances.Sort();		
		}
		/// <summary>
		/// Look for 'n' neighbour objects to a given point. Each neighbour must be unique..
		/// </summary>
		/// <param name="inputVector">vector of objects to find neighbours of</param>
		/// <param name="neighbourCount">number of neighbours to find</param>
		/// <param name="distances">array of distances</param>
        public void UniqueSearch(ITreeReferenceable[] inputVector, int neighbourCount, ref List<KDTreePoint> distances)  //   
		{
			double [] newInputVector = new double[inputVector.Length];
			for(int n = 0; n < inputVector.Length; n++)
				newInputVector[n] = inputVector[n].GetValue();
			UniqueSearch(newInputVector,neighbourCount, ref distances);
		}


		/// <summary>
		/// find the neighbours within a distance of a point
		/// </summary>
		/// <remarks>Note that the distance array contains the squares of the distances found.</remarks>
		/// <param name="inputVector">Vector to find neighbours of</param>
		/// <param name="maxDistance">distance to search over</param>
		/// <param name="neighbourCount">number of neighbours found within distance</param>
		/// <param name="distances">array of squared distances of neighbours from the input vector</param>
		/// <param name="maxNeighbourCount">maximum permissable size of arrays</param>
        public void SearchByDistance(double[] inputVector, double maxDistance, out int neighbourCount, ref List<KDTreePoint> distances, int maxNeighbourCount) 
		{ 
			double[] lower = new double[inputVector.Length];
			double[] upper = new double[inputVector.Length];
			for(int n = 0; n < inputVector.Length; n++)
			{
				lower[n] = inputVector[n] - maxDistance;
				upper[n] = inputVector[n] + maxDistance;
			}
			KDNode[] list = this.range(lower,upper);
			foreach(KDNode node in list)
			{
				double dist = HyperPoint.eucdist(new HyperPoint(inputVector),node.k);
				if(dist <= maxDistance)
				{
					distances.Add(new KDTreePoint(node.v,dist*dist));
				}
			}
			distances.Sort();
			neighbourCount = distances.Count;
	}
		/// <summary>
		/// find the neighbour objects within a distance of a point
		/// </summary>
		/// <remarks>Note that the distance array contains the squares of the distances found.</remarks>
		/// <param name="inputVector">Vector to find neighbours of</param>
		/// <param name="maxDistance">distance to search over</param>
		/// <param name="neighbourCount">number of neighbours found within distance</param>
		/// <param name="distances">array of squared distances of neighbours from the input vector</param>
		/// <param name="maxNeighbourCount">maximum permissable size of arrays</param>
		public void SearchByDistance(ITreeReferenceable [] inputVector,double maxDistance, out int neighbourCount,ref List<KDTreePoint> distances,int  maxNeighbourCount) 
		{
			double [] newInputVector = new double[inputVector.Length];
			for(int n = 0; n < inputVector.Length; n++)
				newInputVector[n] = inputVector[n].GetValue();
			SearchByDistance(newInputVector, maxDistance, out neighbourCount, ref distances,maxNeighbourCount);
		}

		/// <summary>
		/// 
		/// </summary>
		/// <param name="lowk"></param>
		/// <param name="uppk"></param>
		/// <returns></returns>
		internal KDNode [] range(double [] lowk, double [] uppk)
		{
			if (lowk.Length != uppk.Length) 
			{
				throw new ArgumentOutOfRangeException();
			}
			else if (lowk.Length != dimensions) 
			{
				throw new ArgumentOutOfRangeException();
			}
            List<KDNode> list = new List<KDNode>();
			KDNode.rsearch(new HyperPoint(lowk), new HyperPoint(uppk),root, 0, dimensions, list);
			return list.ToArray();
		}

	}
	/// <summary>
	/// Represents data about a point in the source data.
	/// </summary>
	public class KDTreePoint : IComparable
	{
		/// <summary>
		/// Constructor
		/// </summary>
		/// <param name="Index">index into the source data</param>
		/// <param name="DistSquared">Squared distance of this point from the search pattern</param>
		internal KDTreePoint(object Index, double DistSquared)
		{
			_index = Index;
			_distSquared = DistSquared;
		}
		/// <summary>
		/// Squared distance of this point from the search pattern
		/// </summary>
		public double distSquared
		{
			get
			{
				return _distSquared;
			}
		}
		private double _distSquared;
		/// <summary>
		/// index into the source data
		/// </summary>
		public object index
		{
			get
			{
				return _index;
			}
		}
		private object _index;

		#region IComparable Members
		/// <summary>
		/// Compares one KDTReePoint to another based on distance.
		/// </summary>
		/// <param name="obj"></param>
		/// <returns></returns>
		public int CompareTo(object obj)
		{
			KDTreePoint point = (KDTreePoint)obj;
			return distSquared.CompareTo(point.distSquared);
		}

		#endregion
	}
}
