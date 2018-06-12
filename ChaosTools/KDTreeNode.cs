using System;

namespace ConceptStrings.ChaosTools
{
	/// <summary>
	/// Implements a node of the KDTree structure
	/// </summary>
	[Serializable]
	internal class KDTreeNode
	{
		internal KDTreeNode()
		{
			bucket = staticBucket;       
			left = null;
			right = null;
		}
		internal KDTreeNode(int wNumPoints,int [] pBucket, int inputIndex)
		{                                                     
			bucket = pBucket;
			bucketCount = wNumPoints;
			bucketIndex = inputIndex;
			left = null;
			right = null;
			discriminator = 0;
			partition = 0.0f; 
			nodeCount++;
		}

		internal KDTreeNode(int wNumPoints,int [] pBucket, int inputIndex, int wDisc, double fPart)
		{ 
			bucket = pBucket;
			bucketCount = wNumPoints;
			bucketIndex = inputIndex;
			left = null;
			right = null;
			discriminator = wDisc;
			partition = fPart; 
			nodeCount++;
		}

		internal int GetBucketCount()
		{
			return bucketCount;
		}
		internal int GetDiscrim()
		{
			return discriminator;
		}
		internal double GetPartition()
		{
			return partition;
		}
		internal KDTreeNode GetLeft()
		{
			return left;
		} 
		internal KDTreeNode GetRight()
		{
			return right;
		} 
		internal int GetBucketVal(int index)         
		{
			if(index >= bucketCount)
				throw new Exception("KDTReeNode.GetBucketVal out of range");
			return bucket[bucketIndex + index];
		}  
		internal void SetLeft(KDTreeNode pLeft)
		{
			left = pLeft;
		}
		internal void SetRight(KDTreeNode pRight)
		{
			right = pRight;
		}   
		internal void SetStaticBucket(int [] pStaticBucket)
		{
			staticBucket = pStaticBucket; 
			bucket = staticBucket;
		}
 
		internal bool IsLeafNode()
		{
			return (bucketCount != 0);
		}	

		internal static int nodeCount = 0;
		internal static int [] staticBucket;
		private int [] bucket;
		private int discriminator;  	// dimension discriminated
		private double partition;    // partition value
		private int bucketIndex;        
		private int bucketCount;
		private KDTreeNode left;
		private KDTreeNode right;
	}
}
