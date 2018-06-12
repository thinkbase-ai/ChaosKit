using System;

namespace ConceptStrings.ChaosTools
{
	/// <summary>
	/// 
	/// </summary>
	internal class NearestNeighborList
	{
		internal static int REMOVE_HIGHEST = 1;
		internal static int REMOVE_LOWEST = 2;

		PriorityQueue m_Queue = null;
		int m_Capacity = 0;

		// constructor
		internal NearestNeighborList(int capacity) 
		{
			m_Capacity = capacity;
			m_Queue = new PriorityQueue(m_Capacity,Double.PositiveInfinity);
		}

		internal double getMaxPriority() 
		{
			if (m_Queue.length()==0) 
			{
				return Double.PositiveInfinity;
			}
			return m_Queue.getMaxPriority();
		}

		internal bool insert(Object obj,double priority) 
		{
			if (m_Queue.length()<m_Capacity) 
			{
				// capacity not reached
				m_Queue.add(obj,priority);
				return true;
			}
			if (priority>m_Queue.getMaxPriority()) 
			{
				// do not insert - all elements in queue have lower priority
				return false;
			}
			// remove object with highest priority
			m_Queue.remove();
			// add new object
			m_Queue.add(obj,priority);
			return true;
		}

		internal bool isCapacityReached() 
		{
		return m_Queue.length()>=m_Capacity;
		}

		internal Object getHighest() 
		{
		return m_Queue.front();
		}

		internal bool isEmpty() 
		{
		return m_Queue.length()==0;
		}

		internal int getSize() 
		{
		return m_Queue.length();
		}

		internal Object removeHighest() 
		{
			// remove object with highest priority
		return m_Queue.remove();
		}
	}
}
