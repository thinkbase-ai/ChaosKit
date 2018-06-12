using System;
using System.Collections.Generic;

namespace ConceptStrings.ChaosTools
{
	/// <summary>
	/// 
	/// </summary>
	internal class KDNode
	{
		internal HyperPoint k;
		internal Object v;
		internal KDNode left, right;
		internal bool deleted;
		internal List<object> duplicates = null;
		private KDNode(HyperPoint key, Object val) 
		{
	
			k = key;
			v = val;
			left = null;
			right = null;
			deleted = false;
		}
		internal static KDNode ins(HyperPoint key, Object val, KDNode t, int lev, int K)
		{
			if (t == null) 
			{
				t = new KDNode(key, val);
			}
	
			else if (key.equals(t.k)) 
			{
				if(t.deleted)
				{
					t.deleted = false;
					t.v = val;
				}
				else
				{
					if(t.duplicates == null)
						t.duplicates = new List<object>();
					t.duplicates.Add(val);
				}
			}

			else if (key.coord[lev] > t.k.coord[lev]) 
			{
				t.right = ins(key, val, t.right, (lev+1)%K, K);
			}
			else 
			{
				t.left = ins(key, val, t.left, (lev+1)%K, K);
			}
	
			return t;
		}

		internal static KDNode srch(HyperPoint key, KDNode t, int K) 
		{

			for (int lev=0; t!=null; lev=(lev+1)%K) 
			{

				if (!t.deleted && key.equals(t.k)) 
				{
					return t;
				}
				else if (key.coord[lev] > t.k.coord[lev]) 
				{
					t = t.right;
				}
				else 
				{
					t = t.left;
				}
			}

			return null;
		}
        internal static void rsearch(HyperPoint lowk, HyperPoint uppk, KDNode t, int lev, int K, List<KDNode> v) 
		{

			if (t == null) return;
			if (lowk.coord[lev] <= t.k.coord[lev]) 
			{
				rsearch(lowk, uppk, t.left, (lev+1)%K, K, v);
			}
			int j;
			for (j=0; j<K && lowk.coord[j]<=t.k.coord[j] && 
				uppk.coord[j]>=t.k.coord[j]; j++) 
				;
			if (j==K) v.Add(t);
			if (uppk.coord[lev] > t.k.coord[lev]) 
			{
				rsearch(lowk, uppk, t.right, (lev+1)%K, K, v);
			}
		}
		internal static void nnbr(KDNode kd, HyperPoint target, HyperRect hr,double max_dist_sqd, int lev, int K,	NearestNeighborList nnl) 
		{

			// 1. if kd is empty then set dist-sqd to infinity and exit.
			if (kd == null) 
			{
				return;
			}

			// 2. s := split field of kd
			int s = lev % K;

			// 3. pivot := dom-elt field of kd
			HyperPoint pivot = kd.k;
			double pivot_to_target = HyperPoint.sqrdist(pivot, target);

			// 4. Cut hr into to sub-hyperrectangles left-hr and right-hr.
			//    The cut plane is through pivot and perpendicular to the s
			//    dimension.
			HyperRect left_hr = hr; // optimize by not cloning
			HyperRect right_hr = (HyperRect) hr.clone();
			left_hr.max.coord[s] = pivot.coord[s];
			right_hr.min.coord[s] = pivot.coord[s];

			// 5. target-in-left := target_s <= pivot_s
			bool target_in_left = target.coord[s] < pivot.coord[s];

			KDNode nearer_kd;
			HyperRect nearer_hr;
			KDNode further_kd;
			HyperRect further_hr;

			// 6. if target-in-left then
			//    6.1. nearer-kd := left field of kd and nearer-hr := left-hr
			//    6.2. further-kd := right field of kd and further-hr := right-hr
			if (target_in_left) 
			{
				nearer_kd = kd.left;
				nearer_hr = left_hr;
				further_kd = kd.right;
				further_hr = right_hr;
			}
				//
				// 7. if not target-in-left then
				//    7.1. nearer-kd := right field of kd and nearer-hr := right-hr
				//    7.2. further-kd := left field of kd and further-hr := left-hr
			else 
			{
				nearer_kd = kd.right;
				nearer_hr = right_hr;
				further_kd = kd.left;
				further_hr = left_hr;
			}

			// 8. Recursively call Nearest Neighbor with paramters
			//    (nearer-kd, target, nearer-hr, max-dist-sqd), storing the
			//    results in nearest and dist-sqd
			nnbr(nearer_kd, target, nearer_hr, max_dist_sqd, lev + 1, K, nnl);

			KDNode nearest = (KDNode) nnl.getHighest();
			double dist_sqd;

			if (!nnl.isCapacityReached()) 
			{
				dist_sqd = Double.MaxValue;
			}
			else 
			{
				dist_sqd = nnl.getMaxPriority();
			}

			// 9. max-dist-sqd := minimum of max-dist-sqd and dist-sqd
			max_dist_sqd = Math.Min(max_dist_sqd, dist_sqd);

			// 10. A nearer point could only lie in further-kd if there were some
			//     part of further-hr within distance sqrt(max-dist-sqd) of
			//     target.  If this is the case then
			HyperPoint closest = further_hr.closest(target);
			if (HyperPoint.eucdist(closest, target) < Math.Sqrt(max_dist_sqd)) 
			{

				// 10.1 if (pivot-target)^2 < dist-sqd then
				if (pivot_to_target < dist_sqd) 
				{

					// 10.1.1 nearest := (pivot, range-elt field of kd)
					nearest = kd;

					// 10.1.2 dist-sqd = (pivot-target)^2
					dist_sqd = pivot_to_target;

					// add to nnl
					if (!kd.deleted) 
					{
						nnl.insert(kd, dist_sqd);
					}

					// 10.1.3 max-dist-sqd = dist-sqd
					// max_dist_sqd = dist_sqd;
					if (nnl.isCapacityReached()) 
					{
						max_dist_sqd = nnl.getMaxPriority();
					}
					else 
					{
						max_dist_sqd = Double.MaxValue;
					}
				}

				// 10.2 Recursively call Nearest Neighbor with parameters
				//      (further-kd, target, further-hr, max-dist_sqd),
				//      storing results in temp-nearest and temp-dist-sqd
				nnbr(further_kd, target, further_hr, max_dist_sqd, lev + 1, K, nnl);
				KDNode temp_nearest = (KDNode) nnl.getHighest();
				double temp_dist_sqd = nnl.getMaxPriority();

				// 10.3 If tmp-dist-sqd < dist-sqd then
				if (temp_dist_sqd < dist_sqd) 
				{

					// 10.3.1 nearest := temp_nearest and dist_sqd := temp_dist_sqd
					nearest = temp_nearest;
					dist_sqd = temp_dist_sqd;
				}
			}

				// SDL: otherwise, current point is nearest
			else if (pivot_to_target < max_dist_sqd) 
			{
				nearest = kd;
				dist_sqd = pivot_to_target;
			}
		}

	}
}
