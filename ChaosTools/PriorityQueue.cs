using System;

namespace ConceptStrings.ChaosTools
{
	/// <summary>
	/// 
	/// </summary>
	internal class PriorityQueue
	{
		private double maxPriority = Double.MaxValue;
		/// <summary>
		/// This contains the list of objects in the queue.
		/// </summary>
		private Object[] data;
		/// <summary>
		/// This contains the list of prioritys in the queue.
		/// </summary>
		private double[] value;
		/// <summary>
		/// Holds the number of elements currently in the queue.
		/// </summary>
		private int count;
		/// <summary>
		/// This holds the number elements this queue can have.
		/// </summary>
		private int capacity;

		/// <summary>
		/// Creates a new PriorityQueue object. The PriorityQueue object allows 
		/// objects to be entered into the queue and to leave in the order of priority 
		/// i.e the highest priority get's to leave first.
		/// </summary>
		internal PriorityQueue() 
		{
			init(20);
		}

		/// <summary>
		/// Creates a new PriorityQueue object. The PriorityQueue object allows 
		/// objects to be entered into the queue and to leave in the order of priority 
		/// i.e the highest priority get's to leave first.
		/// </summary>
		/// <param name="capacity">the initial capacity of the queue before a resize</param>	
		internal PriorityQueue(int capacity) 
		{
			init(capacity);
		}

		/// <summary>
		/// Creates a new PriorityQueue object. The PriorityQueue object allows 
		/// objects to be entered into the queue and to leave in the order of priority 
		/// i.e the highest priority get's to leave first.
		/// </summary>
		/// <param name="capacity">the initial capacity of the queue before a resize</param>
		/// <param name="maxPriority">param maxPriority is the maximum possible priority for an object</param>
		internal PriorityQueue(int capacity, double maxPriority) 
		{
			this.maxPriority = maxPriority;
			init(capacity);
		}


		/// <summary>
		/// This is an initializer for the object. It basically initializes
		///  an array of long called value to represent the prioritys of the objects, 
		///  it also creates an array of objects to be used in parallel with the array 
		///  of longs, to represent the objects entered, these can be used to sequence the data.
		/// </summary>
		/// <param name="size">the initial capacity of the queue, it can be resized</param>
		private void init(int size) 
		{
			capacity = size;
			data = new Object[capacity + 1];
			value = new double[capacity + 1];
			value[0] = maxPriority;
			data[0] = null;
		}
		/**
		 * This function adds the given object into the <code>PriorityQueue</code>,
		 * its priority is the long priority. The way in which priority can be
		 * associated with the elements of the queue is by keeping the priority
		 * and the elements array entrys parallel.
		 *
		 * @param 
		 * @param 
		 * <code>PriorityQueue</code>
		 */
		/// <summary>
		/// This function adds the given object into the PriorityQueue its priority
		/// is the long priority. The way in which priority can be associated with the
		/// elements of the queue is by keeping the priority and the elements array entrys parallel.
		/// </summary>
		/// <param name="element">element is the object that is to be entered into this PriorityQueue</param>
		/// <param name="priority">priority this is the priority that the object holds in the PriorityQueue</param>
		internal void add(Object element, double priority) 
		{
			if (count++ >= capacity) 
			{
				expandCapacity();
			}
			/* put this as the last element */
			value[count] = priority;
			data[count] = element;
			bubbleUp(count);
		}

		/**
		 * Remove is a function to remove the element in the queue with the
		 * maximum priority. Once the element is removed then it can never be
		 * recovered from the queue with further calls. The lowest priority
		 * object will leave last.
		 *
		 * @return the object with the highest priority or if it's empty
		 * null
		 */
		/// <summary>
		/// Remove is a function to remove the element in the queue with the
		///  maximum priority. Once the element is removed then it can never be 
		///  recovered from the queue with further calls. The lowest priority object will leave last.
		/// </summary>
		/// <returns>the object with the highest priority or if it's empty null</returns>
		internal Object remove() 
		{
			if (count == 0)
				return null;
			Object element = data[1];
			/* swap the last element into the first */
			data[1] = data[count];
			value[1] = value[count];
			/* let the GC clean up */
			data[count] = null;
			value[count] = 0L;
			count--;
			bubbleDown(1);
			return element;
		}

		internal Object front() 
		{
			return data[1];
		}

		internal double getMaxPriority() 
		{
			return value[1];
		}

	
		/// <summary>
		/// Bubble down is used to put the element at subscript 'pos' 
		/// into it's rightful place in the heap (i.e heap is another name 
		/// for PriorityQueue). If the priority of an element at subscript 
		/// 'pos' is less than it's children then it must be put under one 
		/// of these children, i.e the ones with the maximum priority must come first.
		/// </summary>
		/// <param name="pos">is the position within the arrays of the element and priority</param>
		private void bubbleDown(int pos) 
		{
			Object element = data[pos];
			double priority = value[pos];
			int child;
			/* hole is position '1' */
			for (; pos * 2 <= count; pos = child) 
			{
				child = pos * 2;
				/* if 'child' equals 'count' then there
				   is only one leaf for this parent */
				if (child != count)

					/* left_child > right_child */
					if (value[child] < value[child + 1])
						child++; /* choose the biggest child */
				/* percolate down the data at 'pos', one level
					  i.e biggest child becomes the parent */
				if (priority < value[child]) 
				{
					value[pos] = value[child];
					data[pos] = data[child];
				}
				else 
				{
					break;
				}
			}
			value[pos] = priority;
			data[pos] = element;
		}

		
		
		/// <summary>
		/// Bubble up is used to place an element relatively low in the queue 
		/// to it's rightful place higher in the queue, but only if it's priority 
		/// allows it to do so, similar to bubbleDown only in the other direction 
		/// this swaps out its parents.
		/// </summary>
		/// <param name="pos">the position in the arrays of the object to be bubbled up</param>
		private void bubbleUp(int pos) 
		{
			Object element = data[pos];
			double priority = value[pos];
			/* when the parent is not less than the child, end*/
			while (value[pos / 2] < priority) 
			{
				/* overwrite the child with the parent */
				value[pos] = value[pos / 2];
				data[pos] = data[pos / 2];
				pos /= 2;
			}
			value[pos] = priority;
			data[pos] = element;
		}

		/// <summary>
		/// This ensures that there is enough space to keep adding elements 
		/// to the priority queue. It is however advised to make the capacity 
		/// of the queue large enough so that this will not be used as it is 
		/// an expensive method. This will copy across from 0 as 'off' equals 
		/// 0 is contains some important data.
		/// </summary>	
		private void expandCapacity() 
		{
			capacity = count * 2;
			Object[] elements = new Object[capacity + 1];
			double[] prioritys = new double[capacity + 1];
			System.Array.Copy(data, 0, elements, 0, data.Length);
			System.Array.Copy(value, 0, prioritys, 0, data.Length);
			data = elements;
			value = prioritys;
		}

		/// <summary>
		/// This method will empty the queue. This also helps garbage 
		/// collection by releasing any reference it has to the elements 
		/// in the queue. This starts from offset 1 as off equals 0 for the elements array.
		/// </summary>
		internal void clear() 
		{
			for (int i = 1; i < count; i++) 
			{
				data[i] = null; /* help gc */
			}
			count = 0;
		}
		/// <summary>
		/// The number of elements in the queue. 
		/// The length indicates the number of elements 
		/// that are currently in the queue.
		/// </summary>
		/// <returns>the number of elements in the queue</returns>
		internal int length() 
		{
			return count;
		}
	}
}
