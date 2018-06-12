using System;
using System.Collections.Generic;


namespace ConceptStrings.ChaosTools
{
    /// <summary>
    /// Holds a sorted list of events and provides the tools to return data organized for training and various measures of complexity
    /// </summary>
    [Serializable]
	internal class EventList
    {
        #region constants
        internal const int ITERATE_GAP_LENGTH = 1;	// gap between validation sets when iterated
		internal const double RTOL = 3.0;      // ratio tolerance for acceptance as an fnn, was 3.0
		internal const double ATOL = 0.3;       // Absolute tolerance for acceptance as an fnn was 0.3 ANE 23/08/05
		internal const double FNN_THRESHOLD = 0.1;
		internal const int MAXDIM = 10;
		internal const int IMPORTLINE = 512;  // length of a line for importing text files 
		internal const int MAXDELAY = 17;     // search up to this separation
		internal const double BANDWIDTH = 0.1;  // initial Value for search radii in MMI estimation 
		internal const double BANDWIDTHSQUARED = BANDWIDTH * BANDWIDTH;
		internal const int FNN_CUTOFF = 20; // min number of points for FNN algorithm
		internal const double MISUCCESSRATIO = 10.0;
		internal const int MAXPOINTS = 10; // number of neighbours to be considered in any search
		internal const int TIMEDIFF  = 10; // distance in time required of a candidate next point
		internal const double SCALEMN = 0.0;
		internal const double SCALEMX = 1.0;
		internal const double ANGLEMAX = 0.3;
        #endregion

        #region constructors
        internal EventList()
		{
			useTradingTimes = false; 
			useStoredOldVal = false;
			_embedDimension = 2;
			_embedSeparation = 1;
			m_bIsFifoFull = false;
			fifoIndex = 0;
			differencing = false;
			m_bParamSet = false;
			sampleTime = new TimeSpan(0,0,5,0);
			m_fLyapunov = 0.0;
			m_fHurst = 0.0;
			granularity = 0.0001;
			useGranularity = false;
			m_fKaplanGlass = 0.0;
			sorted = false;
			events = new List<Event>();
			_average = 0.0;
			_stdDev = 0.0;
		}
   	
		internal EventList(int wEmbDim, int wEmbSep)
		{ 
			useStoredOldVal = false;
			useTradingTimes = false;
			_embedDimension = wEmbDim;
			_embedSeparation = wEmbSep;
			fifo = new double[_embedDimension * _embedSeparation + 1];
			m_bIsFifoFull = false;
			fifoIndex = 0;
			differencing = false;
			m_bParamSet = false;
			sampleTime = new TimeSpan(0,0,5,0);
			m_fLyapunov = 0.0;
			m_fHurst = 0.0;
			granularity = 0.0001;
			useGranularity = false;
			m_fKaplanGlass = 0.0;
			sorted = false;
            events = new List<Event>();
			_average = 0.0;
			_stdDev = 0.0;
		}

		internal EventList(EventList a)
		{
			differencing = a.differencing;
			_embedDimension = a._embedDimension;
			_embedSeparation = a._embedSeparation;
			fifo = new double[_embedDimension * _embedSeparation + 1];
			maxVal = a.maxVal;
			minVal = a.minVal;
			maxChange = a.maxChange;
			m_bIsFifoFull = a.m_bIsFifoFull;
			fifoIndex = a.fifoIndex;
			m_pos = a.m_pos;
			sampleTime = a.sampleTime; 
			m_bParamSet = a.m_bParamSet;
			tradingTimes = a.tradingTimes;
			useTradingTimes = a.useTradingTimes; 
			useStoredOldVal = a.useStoredOldVal;	
			useGranularity = a.useGranularity;
			oldVal = a.oldVal;
			m_fLyapunov = a.m_fLyapunov;
			m_fHurst = a.m_fHurst;
            events = new List<Event>();
			for(int n = 0; n < a.events.Count; n++)
			{
				events.Add(a.events[n]);
			}
			events.Sort();
			sorted = true;
			_average = 0.0;
			_stdDev = 0.0;
		}

        #endregion

        #region fields
        protected int _embedDimension;
        protected int _embedSeparation;
        protected double maxVal;
        protected double minVal;
        protected double maxChange;
        protected double oldVal;                // for differencing
        protected bool useStoredOldVal;
        protected bool m_bParamSet;
        internal TradingTimes tradingTimes;
        // not archived
        protected bool m_bIsFifoFull;
        protected int fifoIndex;
        internal double[] fifo;
        protected int m_pos;
        protected double m_fLyapunov;
        protected double m_fHurst;
        protected double m_fKaplanGlass;
        protected List<Event> events;
        protected bool sorted;
        protected double _average;
        protected double _stdDev;
        internal TimeSpan deadPeriodThreshold = TimeSpan.Zero;
        internal bool deadPeriodDetectionEnabled = false;
        private bool inDeadPeriod = false;
        private bool interpolated = false;


        #endregion

        #region properties

        public int embedDimension
        {
            set
            {
                if (value >= 16 || value < 1)
                    throw new Exception("Embed dimension out of range");
                _embedDimension = value;
                fifo = new double[_embedDimension * _embedSeparation + 1];
                m_bIsFifoFull = false;
                fifoIndex = 0;
                m_bParamSet = true;
            }
            get
            {
                return _embedDimension;
            }
        }

        public int embedSeparation
		{
			get
			{
				return _embedSeparation;
			}
			set
			{
				_embedSeparation = value;
			}
		}

        public double average
        {
            get
            {
                return _average;
            }
        }

        public double stdDev
        {
            get
            {
                return _stdDev;
            }
        }

        public TimeSpan sampleTime { get; set; }

        public bool useTradingTimes { get; set; }

        public bool differencing { get; set; }

        public double granularity { get; set; }

        public bool useGranularity { get; set; }

        public bool MovingAverage { get; set; }
        /// <summary>
        /// Access the sampled values in the temporal database by sample time index
        /// </summary>
        /// <remarks>0 = first element in the database, Count - 1 = last, and values between these points are interpolated</remarks>
        public double this[int index]
        {
            get
            {
                if (index == 0)
                    return ((Event)events[index]).GetEventValue();
                int count = Count;
                if (index >= count - 1)
                    return ((Event)events[events.Count - 1]).GetEventValue();
                long ticks = (long)index * (((Event)events[events.Count - 1]).GetEventTime().Ticks - ((Event)events[0]).GetEventTime().Ticks) / (long)(count - 1);
                DateTime SampleTime = this.GetEarliestEvent() + new TimeSpan(ticks);
                double res;
                GetValue(out res, SampleTime);
                return res;
            }
        }

        /// <summary>
        /// Returns the number of samples in the temporal database.
        /// </summary>
        /// <exception cref="Exception">if sample time not set</exception>
        public int Count
        {
            get
            {
                if (events.Count == 0)
                    return 0;
                if (sampleTime.Ticks == 0)
                    throw new Exception("Sample time not set");
                return (int)(new TimeSpan(((Event)events[events.Count - 1]).GetEventTime().Ticks - ((Event)events[0]).GetEventTime().Ticks).Ticks / sampleTime.Ticks) + 1;
            }
        }

        #endregion

        #region static methods

        internal static double MIFromProbs(int length, double[] pX, double[] pY, double[] pXY)
        {
            double dMISum = 0.0;
            for (int i = 0; i < length; i++)
            {
                double dTemp = pX[i] * pY[i];
                if (dTemp > 0.0)
                {
                    double dLocMI = pXY[i] / dTemp;

                    if (dLocMI > 0.0)
                        dMISum += Math.Log(dLocMI) * 1.442695040889;
                }
            }// equates to log to the base 2 

            return (dMISum / (double)length);

        }

        /// <summary>
        ///  The kernel density function
        ///  </summary>  
        ///  <remarks>Approximates a delta function:
        ///  approx = delta( x_i - x_j ).
        ///  This is the multidimensional Epanechnikov Kernel:
        ///  (1 - |X/h|^2 )  * (d+2)/(Pi^(d/2) * Gamma( d/2 + 1) *h^d)
        /// </remarks>
        /// <param name="distSqred"></param>
        /// <param name="hSqr"></param>
        /// <param name="h">bandwidth</param>
        /// <param name="dim"></param>
        /// <returns></returns>
        internal static double KernelFunction(double distSqred, double hSqr, double h, int dim)
        {
            if (dim > 10)
                throw new Exception("Kernel Function: dim cannot exceed 10");
            double dX = distSqred / hSqr;
            if (dX > 1.0)
                return 0.0;
            double dTop = (dim + 2) * (1.0 - dX);
            double dBottom = 1.0;
            switch (dim)
            {
                case 1:
                    dBottom = 4.0 * h;
                    break;
                case 2:
                    dBottom = 6.283185308 * h * h;
                    break;
                case 3:
                    dBottom = 8.377580412 * Math.Pow(h, dim);
                    break;
                case 4:
                    dBottom = 9.869604404 * Math.Pow(h, dim);
                    break;
                case 5:
                    dBottom = 10.52757803 * Math.Pow(h, dim);
                    break;
                case 6:
                    dBottom = 10.33542556 * Math.Pow(h, dim);
                    break;
                case 7:
                    dBottom = 9.449531945 * Math.Pow(h, dim);
                    break;
                case 8:
                    dBottom = 8.117424257 * Math.Pow(h, dim);
                    break;
                case 9:
                    dBottom = 6.597017809 * Math.Pow(h, dim);
                    break;
                default:
                    dBottom = 5.100328084 * Math.Pow(h, dim);
                    break;
            }
            return (dTop / dBottom);
        }

        internal static void GeneralKernelEstimator(int wLength, int dimension, double[,] X, KDTree pTree, double[] pHSqreds, ref double[] P)
        {
            List<KDTreePoint> distances = new List<KDTreePoint>();
            for (int i = 0; i < wLength; i++)
            {
                double h2 = pHSqreds[i];
                double h = Math.Sqrt(h2);
                int foundCount;
                double[] sourceVector = new double[dimension];
                for (int p = 0; p < dimension; p++)
                    sourceVector[p] = X[i, p];
                pTree.SearchByDistance(sourceVector, h, out foundCount, ref distances, wLength);
                for (int j = 0; j < foundCount; j++)
                {
                    int idx = (int)((KDTreePoint)distances[j]).index;
                    if (idx < wLength)
                        P[idx] += KernelFunction(((KDTreePoint)distances[j]).distSquared, h2, h, dimension);
                }
            }
            for (int i = 0; i < wLength; i++)
                P[i] = P[i] / (double)wLength;
        }

        /// <summary>
        /// Calculates the mutual information in bits between two vectors of doubles of equal length
        /// </summary>
        /// <param name="x">one array</param>
        /// <param name="y">the other</param>
        internal static double MutualInformation(double[] x, double[] y)
        {
            int length = Math.Min(x.Length, y.Length);
            if (length == 0)
                throw new Exception("Mutual information: overlap between two vectors is 0");
            double[,] pVector1 = Normalise(x);
            double[,] pVector2 = Normalise(y);

            KDTree pTree1 = new KDTree(pVector1);
            KDTree pTree2 = new KDTree(pVector2);
            double[,] pVector3 = new double[length, 2];

            for (int i = 0; i < length; i++)
            {
                pVector3[i, 0] = x[i];
                pVector3[i, 1] = y[i];
            }
            KDTree pTree3 = new KDTree(pVector3);

            double[] pHSqred1 = new double[length];
            double[] pHSqred2 = new double[length];

            // For a first cut 

            double[] pXY = new double[length];
            double[] pX = new double[length];
            double[] pY = new double[length];

            for (int i = 0; i < length; i++)
            {
                pXY[i] = 0.0;
                pX[i] = 0.0;
                pY[i] = 0.0;
                pHSqred1[i] = BANDWIDTHSQUARED;
            }
            GeneralKernelEstimator(length, 2, pVector3, pTree3, pHSqred1, ref pXY);
            GeneralKernelEstimator(length, 1, pVector1, pTree1, pHSqred1, ref pX);
            GeneralKernelEstimator(length, 1, pVector2, pTree2, pHSqred1, ref pY);
            double res1 = MIFromProbs(length, pX, pY, pXY);
            double[] pXY2 = new double[length];
            double[] pX2 = new double[length];
            double[] pY2 = new double[length];

            AdaptBandwidths(length, pXY, pHSqred1, ref pHSqred2);
            GeneralKernelEstimator(length, 1, pVector1, pTree1, pHSqred2, ref pX2);
            GeneralKernelEstimator(length, 1, pVector2, pTree2, pHSqred2, ref pY2);
            GeneralKernelEstimator(length, 2, pVector3, pTree3, pHSqred2, ref pXY2);
            return MIFromProbs(length, pX2, pY2, pXY2);
        }

        private static double[,] Normalise(double[] vector)
        {
            if (vector.Rank > 1)
                throw new Exception("Cannot normalise multidimensional arrays");
            double fMean = 0.0;
            for (int i = 0; i < vector.Length; i++)
                fMean += vector[i];
            fMean /= (double)vector.Length;
            // calculate the std deviation.
            double fStdDev = 0.0;
            for (int i = 0; i < vector.Length; i++)
                fStdDev += (vector[i] - fMean) * (vector[i] - fMean);

            fStdDev = Math.Sqrt(fStdDev / (double)vector.Length);
            double[,] vectorOut = new double[vector.Length, 1];
            if (fStdDev != 0.0)
            {
                // rescale the data to have zero mean and unity std dev.
                for (int i = 0; i < vector.Length; i++)
                    vectorOut[i, 0] = (vector[i] - fMean) / fStdDev;
            }
            else
            {
                for (int i = 0; i < vector.Length; i++)
                    vectorOut[i, 0] = vector[i];
            }
            return vectorOut;
        }

        private static void AdaptBandwidths(int N, double[] P, double[] hin, ref double[] hout)
        {
            double logsum = 0.0;
            for (int i = 0; i < N; i++)
                logsum += Math.Log(P[i]);
            double g = Math.Exp(logsum / N);
            /* g is the geometric mean of all the densities */
            for (int i = 0; i < N; i++)
            {
                double newhsqred = hin[i] * g / P[i];
                hout[i] = newhsqred;
            }  /* FOR i */
        }

        /// <summary>
        /// Plots the best straight line between the points
        /// </summary>
        /// <remarks>Numerical recipes in C, page 524
        /// pData is arranged as y,x pairs of points	in that order. wDataCount is the no of pairs.
        /// The offset and gradient is calculated of the regression line, the offset is returned by reference.</remarks>
        /// <param name="data"></param>
        /// <param name="offset">offset of line</param>
        /// <returns>Slope of line</returns>
        internal static double LinearRegression(double[] data, out double offset)
        {
            double SumX = 0.0;
            double SumY = 0.0;
            double SumXX = 0.0;
            double SumXY = 0.0;
            if (data.Length == 0)
            {
                offset = 0.0;
                return 0.0;
            }
            for (int n = 0; n < data.Length / 2; n++)
            {
                SumX += data[2 * n + 1];
                SumY += data[2 * n];
                SumXX += data[2 * n + 1] * data[2 * n + 1];
                SumXY += data[2 * n + 1] * data[2 * n];
            }
            double Delta = (double)data.Length / 2 * SumXX - SumX * SumX;
            if (Delta == 0.0)
            {
                offset = 0.0;
                return 0.0;
            }
            offset = (SumXX * SumY - SumX * SumXY) / Delta;
            return ((double)data.Length / 2 * SumXY - SumX * SumY) / Delta;
        }

        #endregion

        #region public methods

        #endregion

        #region internal methods

        internal void AddEvent(double Value, DateTime t)
        {
            Event pNewEvent = new Event(Value, t);
            if (events.Count == 0)
            {
                maxVal = Value;
                minVal = Value;
            }
            else
            {
                maxVal = Math.Max(maxVal, Value);
                minVal = Math.Min(minVal, Value);
            }
            events.Add(pNewEvent);
            sorted = false;
        }
        /// <summary>
        /// finds the Value of the series at the time given,
        /// </summary>
        /// <remarks>random access</remarks>
        /// <param name="Value"></param>
        /// <param name="t"></param>
        /// <returns>false if empty or time is earlier than range held</returns>
        internal bool GetValue(out double Value, DateTime t)
        {
            Value = 0.0;
            if (!sorted)
            {
                events.Sort();
                sorted = true;
            }
            if (events.Count == 0)
                return false;
            else if (((Event)events[0]).GetEventTime() > t)
                return false; // t is earlier than the range of data.
            else if (((Event)events[events.Count - 1]).GetEventTime() < t)// t is later than the range of data
            {
                Value = ((Event)events[events.Count - 1]).GetEventValue();
                m_pos = events.Count - 1;
            }
            else // Use BinarySearch
            {
                Event tempEvent = new Event(0.0, t);
                int index = events.BinarySearch(tempEvent);
                if (index < 0)
                {
                    index = ~index;
                    index--;
                    // Check for dead period
                    Event previous = (Event)events[index];
                    Event next = (Event)events[index + 1];
                    if (this.deadPeriodDetectionEnabled)
                    {
                        TimeSpan difference = next.GetEventTime() - previous.GetEventTime();
                        if (difference > this.deadPeriodThreshold)
                            inDeadPeriod = true;
                        else
                            inDeadPeriod = false;
                    }
                    else
                        inDeadPeriod = false;

                    //					inDeadPeriod = this.deadPeriodDetectionEnabled && (next.GetEventTime() - previous.GetEventTime() >this.deadPeriodThreshold);

                    if (interpolated)
                        Value = Interpolate(previous, next, t);
                    else
                        Value = previous.GetEventValue();
                }
                else
                    Value = ((Event)events[index]).GetEventValue();
                m_pos = index;
            }
            return true;
        }
        /// <summary>
        /// finds the Value of the series at the time given,
        /// </summary>
        /// <remarks>Should be called after GetValue, starts search at previous location</remarks>
        /// <param name="Value">The required value by reference</param>
        /// <param name="t">Time to look for</param>
        /// <returns>false if empty or time precedes database contents</returns>
        internal bool GetNextValue(out double Value, DateTime t)
        {
            Value = 0.0;
            if (!sorted)
            {
                events.Sort();
                sorted = true;
            }
            // finds the Value of the series at the time given,
            // must be called after GetValue
            if (events.Count == 0)
                return false;
            else if (((Event)events[0]).GetEventTime() > t)
                return false; // t is earlier than the range of data.
            else if (((Event)events[events.Count - 1]).GetEventTime() < t)// t is later than the range of data
            {
                Value = ((Event)events[events.Count - 1]).GetEventValue();
                m_pos = events.Count - 1;
            }
            else // Use BinarySearch
            {
                Event tempEvent = new Event(0.0, t);
                //				int index = events.BinarySearch(m_pos,events.Count - m_pos -1,tempEvent,null);
                int index = events.BinarySearch(tempEvent);
                if (index < 0)
                {
                    index = ~index;
                    index--;
                    // Check for dead period
                    Event previous = (Event)events[index];
                    Event next = (Event)events[index + 1];
                    if (this.deadPeriodDetectionEnabled)
                    {
                        TimeSpan difference = next.GetEventTime() - previous.GetEventTime();
                        if (difference > this.deadPeriodThreshold)
                            inDeadPeriod = true;
                        else
                            inDeadPeriod = false;
                    }
                    else
                        inDeadPeriod = false;

                    //					inDeadPeriod = this.deadPeriodDetectionEnabled && (next.GetEventTime() - previous.GetEventTime() >this.deadPeriodThreshold);
                    if (interpolated)
                        Value = Interpolate(previous, next, t);
                    else
                        Value = previous.GetEventValue();
                }
                else
                    Value = ((Event)events[index]).GetEventValue();
                m_pos = index;
            }
            return true;
        }

        internal bool IsDifferenced()
        {
            return differencing;
        }
        internal void SetDifferencing(bool diff)
        {
            differencing = diff;
        }


        internal void SetTradingTimes(TradingTimes t)
        {
            if (t != null)
            {
                tradingTimes = t;
                useTradingTimes = true;
                if (!tradingTimes.twentyFourHrs)  // not 24hrs so use trading times
                    return;
                for (int i = 0; i < 7; i++)
                    if (tradingTimes.tradingDays[i] == false)  // at least one day not traded, so ditto
                        return;
            }
            useTradingTimes = false;
        }

        internal bool IsTradingTime(DateTime time)
        {
            if (!tradingTimes.IsTradingTime(time, sampleTime))
                return false;
            return !DeadBigger();
        }


        // This function checks to see if the intra-day dead period is bigger than the sample time             
        internal bool DeadBigger()
        {
            TimeSpan DeadSize = (tradingTimes.tradingEnd - tradingTimes.tradingStart);

            if (tradingTimes.tradingStart > tradingTimes.tradingEnd) // reverse orientation
                DeadSize += new TimeSpan(1, 0, 0, 0);

            return (DeadSize > sampleTime);
        }



        internal void ClearAll()
        {
            events.Clear();
        }

        internal void ResetMaxMin()
        {
            maxChange = 0.0;
        }

        /// <summary>
        /// Loads the embedding fifo
        /// </summary>
        /// <remarks>The use of trading times gives rise to a number of changes to LoadFifo
        ///
        /// Firstly, If we're using trading times, we return false if the SampleTime turns out to be
        /// an invalid trading time (this Value is not loaded into the fifo).
        ///
        /// Secondly, if we're differencing, we have to be sure that both the SampleTime Value and the
        /// previous Value used to produce the difference are valid. We used to define the previous Value
        /// as the one at time (SampleTiime - ctsStep), but if trading times are used, we cannot 
        /// guarantee that this is valid. To overcome this problem, we now access a previous valid Value 
        /// stored as event list member data. This Value is kept up to date by LoadFifo and its 
        /// sister function, EmbedForChaosMetrics.  However, the first time these functions are
        /// called (eg. at the start of training/validation pattern production) we have to search 
        /// backwards for the previous Value since we do not have an appropriate stored Value to use. 
        /// </remarks>
        /// <param name="SampleTime"></param>
        /// <param name="ctsStep"></param>
        /// <returns></returns>
        internal bool LoadFifo(DateTime SampleTime, TimeSpan ctsStep)
        {
            double val = 0.0;
            double oldVal = 0.0;

            if (useTradingTimes && !IsTradingTime(SampleTime)) // returns false if trading times 	
                return false;                                     // selected and this SampleTime invalid

            if (!m_bIsFifoFull)
            {
                if (!GetValue(out val, SampleTime))
                {							// if the request is outside the time scale for this item
                    m_bIsFifoFull = false; 		// means patterns can't be generated
                    return true;
                }
                if (differencing)
                {
                    if (!useStoredOldVal)
                    {
                        if (useTradingTimes)
                        {
                            bool ok = false;

                            for (SampleTime -= ctsStep; GetValue(out oldVal, SampleTime); SampleTime -= ctsStep)
                                if (IsTradingTime(SampleTime)) { ok = true; break; }

                            if (!ok) return false;
                        }
                        else
                        {
                            if (!GetValue(out oldVal, SampleTime - ctsStep))
                                return true;
                        }
                    }

                    double fDiff = val - oldVal;
                    oldVal = val;
                    val = ScaleToNet(fDiff, -maxChange, maxChange);
                }
                else
                {
                    oldVal = val;
                    val = ScaleToNet(val, minVal, maxVal);
                }
            }
            else
            {
                if (!GetNextValue(out val, SampleTime))
                {// if the request is outside the time scale for this item
                    m_bIsFifoFull = false; // means patterns can't be generated
                    return true;
                }
                if (differencing)
                {
                    if (!useStoredOldVal)
                    {
                        if (useTradingTimes)
                        {
                            bool ok = false;

                            for (SampleTime -= ctsStep; GetValue(out oldVal, SampleTime); SampleTime -= ctsStep)
                                if (IsTradingTime(SampleTime)) { ok = true; break; }

                            if (!ok) return false;
                        }
                        else
                        {
                            if (!GetValue(out oldVal, SampleTime - ctsStep))
                                return true;
                        }
                    }

                    double fDiff = val - oldVal;
                    oldVal = val;
                    val = ScaleToNet(fDiff, -maxChange, maxChange);
                }
                else
                {
                    oldVal = val;
                    val = ScaleToNet(val, minVal, maxVal);
                }
            }

            fifo[fifoIndex] = val;
            useStoredOldVal = true;

            if (fifoIndex == 0)
                fifoIndex = fifo.Length - 1; // wrap array
            else
            {
                if (fifoIndex == 1 || fifo.Length == 1)
                    m_bIsFifoFull = true;
                fifoIndex--;
            }
            return true;
        }


        // This is the backwards-fill version of load fifo which is used to fill up all of the validation
        // fifos by stepping back over dead periods if necessary.
        /// <summary>
        /// Fills the embedding Fifo
        /// </summary>
        /// <param name="SampleTime"></param>
        /// <param name="ctsStep"></param>
        /// <returns></returns>

        internal bool FillFifo(DateTime SampleTime, TimeSpan ctsStep)
        {
            double val = 0.0;
            double oldVal = 0.0;

            if (useTradingTimes && !IsTradingTime(SampleTime)) // returns false if trading times 	
                return false;                                     // selected and this SampleTime invalid

            if (!GetValue(out val, SampleTime))
            {	// if the request is outside the time scale for this item
                // means patterns can't be generated
                throw new Exception("Chaos Analysis: Can't perform embedding for PredictNext with data and parameters as set.");
            }

            if (differencing)
            {
                if (useTradingTimes)
                {
                    bool ok = false;

                    for (SampleTime -= ctsStep; GetValue(out oldVal, SampleTime); SampleTime -= ctsStep)
                        if (IsTradingTime(SampleTime)) { ok = true; break; }

                    if (!ok)
                        throw new Exception("Chaos Analysis: Can't perform embedding for PredictNext with data and parameters as set.");
                }
                else
                    if (!GetValue(out oldVal, SampleTime - ctsStep))
                        throw new Exception("Chaos Analysis: Can't perform embedding for PredictNext with data and parameters as set.");

                double fDiff = val - oldVal;
                val = ScaleToNet(fDiff, -maxChange, maxChange);
            }
            else
                val = ScaleToNet(val, minVal, maxVal);

            fifo[fifoIndex] = val;

            if (fifoIndex == fifo.Length - 1 || fifo.Length == 1) // now full
                return true;
            else
                fifoIndex++;

            return false;
        }



        internal bool CalcEmbeddingDimension()
        {
            // Uses the technique of false Nearest Neighbours to calculate
            // the embedding  dimension. It is quite possible not to get a
            // definite, or in-range answer. If so the function returns false.
            // A definite result is when The percentage of FNNs is lower than 1%
            // for a valid embedding dimension. Noise, or higher dimensional chaos
            // than we can handle, can cause non-convergence. In this case we take
            // the lowest answer if below 10%, and return false otherwise.


            double[] fRes = new double[MAXDIM];
            fRes[0] = 1.0;
            for (int d = 1; d < MAXDIM; d++)
            {
                fRes[d] = CalcFalseNeighbours(d, _embedSeparation);
            }
            //look for rapid drop
            for (int d = 1; d < MAXDIM; d++)
            {
                if (fRes[d - 1] - fRes[d] > 0.5 || fRes[d] < 0.4/*FNN_THRESHOLD*/)
                {
                    _embedDimension = Math.Max(d, 2);
                    return true;
                }
            }
            //look for small gradient
            for (int d = 1; d < MAXDIM; d++)
            {
                if (fRes[d - 1] - fRes[d] < 0.02 && fRes[d] < 0.8)
                {
                    _embedDimension = Math.Max(d, 2);
                    return true;
                }
            }
            // If not, look for the first minimum 
            double fSmallest = 1.0;
            int wBestIndex = 0;
            for (int d = 1; d < MAXDIM; d++)
            {
                if (fSmallest > fRes[d])
                {
                    fSmallest = fRes[d];
                    wBestIndex = d;
                }
                else if (fSmallest < fRes[d])
                {
                    break;
                }
            }
            if (fSmallest < 0.8)
            {
                _embedDimension = Math.Max(wBestIndex, 2);
                return true;
            }
            return false;
        }


        internal double[,] EmbedForChaosMetrics(int dimension, int separation, out int patterns)
        {
            if (!sorted)
            {
                events.Sort();
                sorted = true;
            }
            if (events.Count == 0)
            {
                patterns = 0;
                return null;
            }
            DateTime EndTime = GetLatestEvent();
            DateTime StartTime = GetEarliestEvent();
            patterns = 0;
            maxChange = 0.0;
            bool isFifoFull = false;
            // indicates data is present for time cursor position.
            bool valid = false;
            double val;
            double previousVal;

            TimeSpan ctsTemp = EndTime - StartTime;
            int fifoSize = (dimension - 1) * separation + 1;
            double[] fifo = new double[fifoSize];
            int estimated_patterns = (int)(new TimeSpan(EndTime.Ticks - StartTime.Ticks).Ticks / sampleTime.Ticks) + 1;
            double[,] embeddedData = new double[estimated_patterns, dimension];
            int dataIndex = 0;
            fifoIndex = 0;
            DateTime ctTemp;

            useStoredOldVal = false;

            for (ctTemp = StartTime; ctTemp <= EndTime; ctTemp += sampleTime)
            {
                if (!useTradingTimes || IsTradingTime(ctTemp))
                {
                    if (!isFifoFull)
                        valid = GetValue(out val, ctTemp);
                    else
                        valid = GetNextValue(out val, ctTemp);

                    if (differencing)
                    {
                        DateTime Time = ctTemp;

                        if (useStoredOldVal)
                            previousVal = oldVal;
                        else
                        {
                            if (useTradingTimes)
                            {
                                bool ok = false;

                                for (Time -= sampleTime; GetValue(out previousVal, Time); Time -= sampleTime)
                                {
                                    if (IsTradingTime(Time))
                                    {
                                        ok = true;
                                        break;
                                    }
                                }
                                if (!ok)
                                    valid = false;
                            }
                            else
                            {
                                if (!GetValue(out previousVal, Time - sampleTime))
                                    valid = false;
                            }
                        }
                        double fDiff = val - previousVal;
                        oldVal = val;
                        val = fDiff;
                        if (fDiff < 0.0f)
                            fDiff = -fDiff;
                        maxChange = Math.Max(fDiff, maxChange);
                    }
                    else //not differencing
                    {
                        oldVal = val;
                    }
                    if (valid)//i.e. if current time point is within the stored series
                    {
                        useStoredOldVal = true;

                        if (fifoSize > 1)
                        {
                            if (this.inDeadPeriod)//reset the fifo. New 22/08/05
                            {
                                isFifoFull = false;
                                fifoIndex = 0;
                            }
                            fifo[fifoIndex] = val;

                            if (fifoIndex == 0)
                                fifoIndex = fifoSize - 1; // wrap array
                            else
                            {
                                if (fifoIndex == 1)
                                    isFifoFull = true;
                                fifoIndex--;
                            }

                            if (isFifoFull)
                            {
                                patterns++;
                                int p = (fifoIndex + 1) % fifoSize;

                                for (int q = dimension; q > 0; q--)
                                {
                                    embeddedData[dataIndex / dimension, q - 1] = fifo[(p + (q - 1) * separation) % fifoSize];
                                    dataIndex++;
                                }
                            }
                        }
                        else // Fifo is length 1
                        {
                            if (!this.inDeadPeriod)//new 22/08/05
                            {
                                patterns++;
                                embeddedData[dataIndex++, 0] = val;
                                isFifoFull = true;
                            }
                        }
                    }
                }
            }
            return embeddedData;

        }


        internal double[,] EmbedForPrediction(int dimension, int separation, out int patterns, out double[] target, bool bAttention, TimeSpan attentionSpan, DateTime endTime, out double[] lastStim)
        {
            // returns stimulus and target arrays. 
            // both arrays are larger than patterns elements. 
            // The patterns'th vector in the returned stimulus array is the validation / test pattern
            lastStim = new double[dimension];
            target = null;
            if (!sorted)
            {
                events.Sort();
                sorted = true;
            }
            if (events.Count == 0)
            {
                throw new Exception("Chaos Analysis: Temporal database empty");
            }
            DateTime startTime = GetEarliestEvent();
            patterns = 0;
            bool isFifoFull = false;
            bool valid = false;
            double val;
            double previousVal;
            if (bAttention)
            {
                if (endTime - attentionSpan > GetEarliestEvent())
                    startTime = endTime - attentionSpan;
            }
            // ensure start time is an integral number of sampleTimes before endTime.
            TimeSpan temp = endTime - startTime;
            startTime += new TimeSpan(temp.Ticks % sampleTime.Ticks);
            //ensure at least 5 samples can be created
            if (temp.Ticks / sampleTime.Ticks < 5)
                throw new Exception("Chaos Analysis: With the current contents of the temporal database, and the choice of sample time and start time there is too little data to create predictions.");
            //set up the embedding fifo.
            int fifoSize = (dimension - 1) * separation + 1;
            fifo = new double[fifoSize];
            // arrays are set at largest possible size to save time calculating the required size.
            int estimated_patterns = (int)(new TimeSpan(endTime.Ticks - startTime.Ticks).Ticks / sampleTime.Ticks);
            double[,] embeddedData = new double[estimated_patterns, dimension];
            int dataIndex = 0;
            //			target = new double [events.Count]; // ANE 30/06/05
            target = new double[estimated_patterns];
            fifoIndex = 0;
            DateTime ctTemp;

            useStoredOldVal = false;

            for (ctTemp = startTime; ctTemp <= endTime; ctTemp += sampleTime)
            {
                if (!useTradingTimes || IsTradingTime(ctTemp))
                {
                    if (!isFifoFull)
                        valid = GetValue(out val, ctTemp);
                    else
                        valid = GetNextValue(out val, ctTemp);

                    if (differencing)
                    {
                        DateTime Time = ctTemp;

                        if (useStoredOldVal)
                            previousVal = oldVal;
                        else
                        {
                            if (useTradingTimes)
                            {
                                bool ok = false;

                                for (Time -= sampleTime; GetValue(out oldVal, Time); Time -= sampleTime)
                                {
                                    if (IsTradingTime(Time))
                                    {
                                        ok = true;
                                        break;
                                    }
                                }
                                if (!ok)
                                    valid = false;
                            }
                            else
                            {
                                if (!GetValue(out oldVal, Time - sampleTime))
                                    valid = false;
                            }
                        }
                        double fDiff = val - oldVal;
                        oldVal = val;
                        val = fDiff;
                        if (fDiff < 0.0f)
                            fDiff = -fDiff;
                        maxChange = Math.Max(fDiff, maxChange);
                    }
                    else
                    {
                        oldVal = val;
                    }

                    if (valid)
                    {
                        useStoredOldVal = true;

                        if (fifoSize > 1)
                        {
                            if (this.inDeadPeriod)//reset the fifo. New 22/08/05
                            {
                                isFifoFull = false;
                                fifoIndex = 0;
                            }
                            fifo[fifoIndex] = val;
                            if (fifoIndex == 0)
                                fifoIndex = fifoSize - 1; // wrap array
                            else
                            {
                                if (fifoIndex == 1)
                                    isFifoFull = true;
                                fifoIndex--;
                            }
                            if (isFifoFull)
                            {
                                patterns++;
                                int p = (fifoIndex + 1) % fifoSize;
                                for (int q = dimension; q > 0; q--)
                                {
                                    int row = dataIndex / dimension;
                                    int column = dataIndex % dimension;
                                    double Value = fifo[(p + (q - 1) * separation) % fifoSize];
                                    if (differencing)
                                        embeddedData[row, column] = ScaleToNet(Value, -maxChange, maxChange);
                                    else
                                        embeddedData[row, column] = ScaleToNet(Value, minVal, maxVal);
                                    dataIndex++;
                                }
                                if (patterns > 1)
                                {
                                    if (differencing)
                                        val = ScaleToNet(val, -maxChange, maxChange);
                                    else
                                        val = ScaleToNet(val, minVal, maxVal);
                                    target[patterns - 2] = val;
                                }
                            }
                        }
                        else        // Fifo is length 1
                        {
                            if (this.inDeadPeriod)// New 22/08/05
                            {
                                patterns++;
                                if (differencing)
                                    embeddedData[dataIndex++, 0] = ScaleToNet(val, 0.0f, maxChange);
                                else
                                    embeddedData[dataIndex++, 0] = ScaleToNet(val, minVal, maxVal);
                                isFifoFull = true;
                                if (patterns > 1)
                                    target[patterns - 2] = val;
                            }
                        }
                    }
                }
            }
            patterns--; // because we can't generate a target for the last pattern we decrease the count
            // however the last pattern in embeddedData, at patterns, is used for predictions if required.
            //For C# we can't resize arrays, so in the case of financial data with trading times
            //we may have less data than the maximum generatable between the start and finish times
            // we need to copy the arrays to remove the waste at the end.
            double[,] compactedData = new double[patterns, dimension];
            Array.Copy(embeddedData, compactedData, patterns * dimension);
            //now grab the last pattern as the stimulus for the first prediction
            for (int n = 0; n < dimension; n++)
                lastStim[n] = embeddedData[patterns, n];
            return compactedData;
        }

        internal double CalcFalseNeighbours(int dim, int sep)
        {
            int patterns;
            int falseCount = 0;
            double range = GetLargestEvent() - GetSmallestEvent();

            if (range == 0.0f || _stdDev == 0.0)
                return 0.0; // no variation in the data

            double[,] pEmbedMatrix = EmbedForChaosMetrics(dim, sep, out patterns);

            if (patterns < FNN_CUTOFF)  // prevents crashing when sep > patterns etc. 
                return 0.0;

            KDTree cTree = new KDTree(pEmbedMatrix);

            double[] searchPattern = new double[dim];

            for (int i = 0; i < patterns - sep; i++)
            {
                List<KDTreePoint> indexes = new List<KDTreePoint>();
                for (int j = 0; j < dim; j++)
                    searchPattern[j] = pEmbedMatrix[i, j];
                cTree.Search(searchPattern, 2, ref indexes);
                KDTreePoint nearestPoint = indexes[1];
                double nearestDist = Math.Sqrt(nearestPoint.distSquared); // The zeroth result is the search pattern
                if (i + sep >= patterns || (int)nearestPoint.index + sep >= patterns)
                    continue;
                double extra = pEmbedMatrix[i + sep, 0] - pEmbedMatrix[(int)nearestPoint.index + sep, 0]; // extra distance for next dim.
                if (extra < 0.0f)
                    extra = -extra; // take the magnitude  
                if (nearestDist != 0.0)    // financial values are discrete, thus identical values are common, the distance can be zero.
                {
                    if (extra / nearestDist > RTOL)   	  // we have a false neighbor 
                        falseCount++;
                    else if (extra / _stdDev > ATOL) // and again
                        falseCount++;
                }
                else
                {
                    if (extra > 0.0)
                        falseCount++;
                    else if ((nearestDist + extra) / _stdDev > ATOL) // and again
                        falseCount++;
                }
            }
            return (double)falseCount / (double)(patterns - sep);
        }
        /// <summary>
        /// Calculate the optimum number of samples between each embedding point
        /// </summary>
        /// <remarks>Uses the first minimum in the Auto-Mutual-Information to calculate
        /// the separation. The logic is that each dimension of the embedded data
        /// should hold as much new information as possible. Separating the dimensions
        /// by the above achieves this. The first minimum is selected, because in practice,
        /// given finite data, the smaller the separation, the more training patterns can be generated.</remarks>
        /// <returns>false for failure</returns>
        internal bool CalcEmbeddingSeparation()
        {
            double[,] pVector1;     //variables now given function scope so they can be tidied
            double[,] pVector2;     //by the exception handling
            double[,] pVector3;
            double[] pHSqred1;
            double[] pHSqred2;

            double[] pXY;
            double[] pXY2;
            double[] pX2;
            double[] pY2;

            KDTree pTree1;
            KDTree pTree2;
            KDTree pTree3;
            //----------------------------------------------------------------------------------------------
            double[] MIScores = new double[18];
            double fMI0 = 0.0;

            int wLength1;
            pVector1 = EmbedForChaosMetrics(1, 1, out wLength1);

            if (wLength1 == 0)
                return false;

            //--------------------------------------------------------------
            // Calculate the mean
            double fMean = 0.0;

            for (int i = 0; i < wLength1; i++)
                fMean += pVector1[i, 0];
            fMean /= (double)wLength1;
            //--------------------------------------------------------------
            // calculate the std deviation.
            double fStdDev = 0.0;

            for (int i = 0; i < wLength1; i++)
                fStdDev += (pVector1[i, 0] - fMean) * (pVector1[i, 0] - fMean);

            fStdDev = Math.Sqrt(fStdDev / (double)wLength1);
            //--------------------------------------------------------------
            if (fStdDev == 0.0)
                return false; // no variance in data - avoid div by zero
            //--------------------------------------------------------------
            // rescale the data to have zero mean and unity std dev.
            for (int i = 0; i < wLength1; i++)
                pVector1[i, 0] = (pVector1[i, 0] - fMean) / fStdDev;
            //--------------------------------------------------------------

            pTree1 = new KDTree(pVector1);
            int wDeltaT;
            //-------------------------------------------------------------------------------------------    
            // M A I N   L O O P    
            try
            {
                for (wDeltaT = 0; wDeltaT <= MAXDELAY; wDeltaT++)	// calculate for different time delays 
                {
                    pVector2 = null;   // Vars initialized to 0 so we can tidy up properly
                    pVector3 = null;   // in the case of an exception
                    pHSqred1 = null;
                    pHSqred2 = null;

                    pXY = null;
                    pXY2 = null;
                    pX2 = null;
                    pY2 = null;

                    pTree2 = null;
                    pTree3 = null;
                    //--------------------------------------------------------------------------------------
                    int wLength2 = wLength1 - wDeltaT;

                    pVector2 = new double[wLength2, 1];
                    pVector3 = new double[wLength2, 2];

                    for (int i = 0; i < wLength2; i++)
                    {
                        pVector2[i, 0] = pVector1[i + wDeltaT, 0];
                        pVector3[i, 0] = pVector1[i, 0];
                        pVector3[i, 1] = pVector2[i, 0];
                    }
                    pTree2 = new KDTree(pVector2);

                    int wLength = Math.Min(wLength1, wLength2);
                    pTree3 = new KDTree(pVector3);

                    double fMIAdaptive;

                    pHSqred1 = new double[wLength];
                    pHSqred2 = new double[wLength];

                    // For a first cut 

                    pXY = new double[wLength];
                    double[] pX = new double[wLength];
                    double[] pY = new double[wLength];


                    for (int i = 0; i < wLength; i++)
                    {
                        pXY[i] = 0.0;
                        pX[i] = 0.0;
                        pY[i] = 0.0;
                        pHSqred1[i] = BANDWIDTHSQUARED;
                    }
                    GeneralKernelEstimator(wLength, 2, pVector3, pTree3, pHSqred1, ref pXY);
                    GeneralKernelEstimator(wLength, 1, pVector1, pTree1, pHSqred1, ref pX);
                    GeneralKernelEstimator(wLength, 1, pVector2, pTree2, pHSqred1, ref pY);
                    double res1 = MIFromProbs(wLength, pX, pY, pXY);
                    double dLogSum = 0.0;

                    for (int i = 0; i < wLength; i++)
                    {
                        if (pXY[i] > 0.0)
                            dLogSum += Math.Log(pXY[i]);
                    }

                    double g = BANDWIDTHSQUARED * Math.Exp(dLogSum / wLength);
                    pXY2 = new double[wLength];
                    pX2 = new double[wLength];
                    pY2 = new double[wLength];

                    for (int i = 0; i < wLength; i++)
                    {
                        if (pXY[i] > 0.0)
                            pHSqred2[i] = (double)(g / pXY[i]);
                        else
                            pHSqred2[1] = BANDWIDTHSQUARED;
                        pXY2[i] = pX2[i] = pY2[i] = 0.0;
                    }
                    GeneralKernelEstimator(wLength, 1, pVector1, pTree1, pHSqred2, ref pX2);
                    GeneralKernelEstimator(wLength, 1, pVector2, pTree2, pHSqred2, ref pY2);
                    GeneralKernelEstimator(wLength, 2, pVector3, pTree3, pHSqred2, ref pXY2);
                    fMIAdaptive = MIFromProbs(wLength, pX2, pY2, pXY2);

                    MIScores[wDeltaT] = fMIAdaptive;

                    if (wDeltaT == 0)
                        fMI0 = fMIAdaptive; // interim success heuristic
                    else if (fMIAdaptive <= 0.0f || fMIAdaptive * MISUCCESSRATIO <= fMI0)
                    {
                        _embedSeparation = wDeltaT;
                        return true;
                    }
                }
            }
            catch // clean up after an exception
            {
                return false;
            }
            double fMin = fMI0;
            int Sep = 0;
            double fMax = 0.0; // this is the first maximum after the deepest minimum

            for (int Count = 1; Count < 18; Count++)
            {
                if (fMin > MIScores[Count])
                {
                    Sep = Count;
                    fMin = MIScores[Count];
                    fMax = MIScores[Count];
                }
                if (fMax < MIScores[Count])
                {
                    fMax = MIScores[Count];
                }
            }

            if (Sep == 0 || Sep == 17) // either unlagged was smallest or the MI was continually decreasing
                return false;
            if (fMax < 1.1f * fMin)
                return false;
            _embedSeparation = Sep;
            return true;
        }

        internal bool SetupEmbedding(bool lyapunov, bool kaplan, bool hurst)
        {
            bool status = true;
            if (!m_bParamSet) // Don't bother if we're not using it.
            {
                // now's the time to calculate embedding dimension and separation
                int length;
                double[,] pVector1 = EmbedForChaosMetrics(1, 1, out length);
                _average = 0.0;
                for (int n = 0; n < length; n++)//changed 22/08/05, was n < pVector1.Length
                    _average += pVector1[n, 0];
                _average /= length;//changed 22/08/05, was _average /= pVector1.Length
                _stdDev = 0.0;
                for (int n = 0; n < length; n++)//changed 22/08/05, was n < pVector1.Length
                {
                    double q = pVector1[n, 0] - _average;
                    _stdDev += q * q;
                }
                _stdDev = Math.Sqrt(_stdDev / length);//changed 22/08/05, was _stdDev /pVector1.Length
                if (length == 0)
                {
                    _embedSeparation = 1;
                    status = false;
                }
                _embedSeparation = 1;
                if (!CalcEmbeddingDimension())
                {
                    _embedDimension = 2;
                    _embedSeparation = 1;
                    status = false;
                }
                /*				FastMutualMI mi = new FastMutualMI(pVector1);
                                if(mi.separation >= 20)//max
                                {
                                    _embedSeparation = 1;
                                    status = false;
                                }
                                else
                                    _embedSeparation = Math.Max(1,mi.separation);*/
                CalculateSeparationUsingPrediction();
                m_bParamSet = true;
            }
            if (lyapunov)
            {
                Lyapunov();

            }
            if (kaplan)
                KaplanGlass();
            if (hurst)
            {
                double[] fTuples = RSAnalysis(false);
                if (fTuples == null)
                    throw new Exception("No variation found in in data - check your source.");
                double fOffset;
                m_fHurst = Math.Min(Math.Max(LinearRegression(fTuples, out fOffset), 0.0), 1.0);
            }
            return status;
        }


        /// <summary>
        /// Calculates RescaledRange
        /// </summary>
        /// <remarks>Based on Chaos and Order in the Capital Markets, page 63
        /// RSArray points to an array of pairs of values. Odd values are Y values, Even are X values
        /// bLogReturns is usually false.</remarks>
        /// <param name="bLogReturns">Use the log of the differences</param>
        /// <returns>RSS Series</returns>
        internal double[] RSAnalysis(bool bLogReturns)
        {

            double[] RSArray;
            if (events.Count == 0 || events.Count < 30)
                return null;
            TimeSpan ctsTemp = new TimeSpan(GetLatestEvent().Ticks - GetEarliestEvent().Ticks);
            int lSampleCount = (int)(ctsTemp.TotalSeconds / sampleTime.TotalSeconds);
            if (lSampleCount < 50)
                return null;
            DateTime ctTemp = GetEarliestEvent();
            double[] pTemp = new double[lSampleCount];
            int tempIndex = 0;
            double oldVal = 1.0;
            if (maxVal - minVal == 0.0)
            {// no variation in data set.
                for (int n = 0; n < lSampleCount; n++)
                    pTemp[n] = 0.0;
            }
            else
            {
                for (int n = 0; n < lSampleCount; n++)
                {
                    double val;
                    if (n == 0)
                        GetValue(out val, ctTemp);
                    else
                        GetNextValue(out val, ctTemp);
                    if (bLogReturns)
                    {
                        double fNewVal = ScaleToNet(val, minVal, maxVal);
                        pTemp[tempIndex++] = Math.Log(fNewVal / oldVal);
                        if (fNewVal != 0.0)
                            oldVal = fNewVal;
                        else
                            oldVal = 1.0;
                    }
                    else
                        pTemp[tempIndex++] = ScaleToNet(val, minVal, maxVal);
                    ctTemp += sampleTime;
                }
            }
            int wNoOfMeasures = (int)(lSampleCount / 10);
            // allocate receiving array
            RSArray = new double[(wNoOfMeasures - 1) * 2];
            int RSCount = 0;
            for (int p = wNoOfMeasures; p >= 2; p--)
            {// all the possible ways of splitting the data
                int q = (int)(lSampleCount / (int)p); // size of clump
                int Offset = 0;
                double fRSAve = 0.0;
                for (int r = 0; r < p; r++)
                {// all the split clumps
                    double fCum = 0.0;
                    double fMin = 0.0;
                    double fMax = 0.0;
                    double fAve = 0.0;
                    for (int s = 0; s < q; s++)
                    { // first pass for average
                        double val = pTemp[Offset + s];
                        fAve += val;
                    }
                    fAve /= (double)q;
                    double fSDev = 0.0;
                    for (int s = 0; s < q; s++)
                    {// second to calculate standard deviation & cumulative stuff
                        double fDiff = pTemp[Offset + s] - fAve;
                        fSDev += fDiff * fDiff;
                        fCum += fDiff;
                        fMin = Math.Min(fCum, fMin);
                        fMax = Math.Max(fCum, fMax);
                    }// Peters is confusing, I assume an average over the clump
                    // not the cumulative average.
                    if (fSDev != 0.0f)
                    {
                        fSDev = Math.Sqrt(fSDev / (double)q);
                        fRSAve += (fMax - fMin) / fSDev;
                    }
                    Offset += q;
                }
                fRSAve /= (double)p;
                if (fRSAve > 0.0f)
                    RSArray[RSCount++] = Math.Log(fRSAve);
                else
                    RSArray[RSCount++] = 0.0;
                RSArray[RSCount++] = Math.Log((double)((int)q * sampleTime.TotalSeconds));
            }
            return RSArray;
        }


        internal double Lyapunov()
        {
            //based on the Algorithm of Wolf, Swift, Swinney & Vastano,
            // Physica 16D 1985, 285-317
            // This assumes that embedding dimension & Separation have already been calculated.
            int patterns;
            int evolve = 1;
            double range = GetLargestEvent() - GetSmallestEvent();
            if (range == 0.0f)
                return (m_fLyapunov = 0.0f); // no variation in the data
            // embed the data apropriately
            double[,] pEmbedMatrix = EmbedForChaosMetrics(_embedDimension, _embedSeparation, out patterns);
            // normalise the data so that distance heuristics are coherent
            for (int k = 0; k < patterns; k++)
            {
                for (int p = 0; p < _embedDimension; p++)
                    pEmbedMatrix[k, p] = ScaleToNet(pEmbedMatrix[k, p], minVal, maxVal);
            }

            KDTree cTree = new KDTree(pEmbedMatrix);

            double[] searchPattern = new double[_embedDimension];
            double[] nearNeighbour = new double[_embedDimension];
            double[] nextNeighbour = new double[_embedDimension];
            int foundCount = 0;
            double lyapunovSum = 0.0;
            for (int i = 0; i < patterns - evolve; i += evolve)
            {
                List<KDTreePoint> indexes = new List<KDTreePoint>();
                double neighbourDist;
                int index;
                // copy the ith vector into searchPattern;    	
                for (int j = 0; j < _embedDimension; j++)
                    searchPattern[j] = pEmbedMatrix[i, j];
                cTree.Search(searchPattern, MAXPOINTS, ref indexes);
                if (foundCount == 0)
                {// first time round don't be choosy
                    KDTreePoint nearestPoint = indexes[1];
                    neighbourDist = Math.Sqrt(nearestPoint.distSquared); // The zeroth result is the search pattern itself
                    index = (int)nearestPoint.index;
                    if (index >= patterns - 1)	 // if its the last in the set we can't track where it goes...
                    {// so use the next nearest. (this can't be the last!)
                        nearestPoint = (KDTreePoint)indexes[2];
                        index = (int)nearestPoint.index;
                        neighbourDist = Math.Sqrt(nearestPoint.distSquared);
                    }
                    foundCount++;
                }
                else
                {
                    // This time choose a neighbour that meets certain criteria.
                    // These criteria arent quite the same as Wolf et. al. We have a presorted set of
                    // nearest neighbours so we search through those for the smallest angle or the first
                    // to be less than ANGLEMAX. So SCALEMX should be 5 times that in the Wolf paper. There
                    // is no point in zooming out, we have the MAXPOINTS best points returned by the KDTREE.
                    double fBestAngle = 3.14;
                    int wBestIndex = MAXPOINTS;
                    for (int p = 1; p < MAXPOINTS; p++)
                    {
                        neighbourDist = Math.Sqrt(((KDTreePoint)indexes[p]).distSquared);
                        index = (int)((KDTreePoint)indexes[p]).index;
                        if (Math.Abs((int)index - (int)i) > TIMEDIFF && index < patterns - 1 && neighbourDist > 0.0)
                        {
                            if (neighbourDist > SCALEMN && neighbourDist < SCALEMX)
                            {// calculate the angle subtended
                                double fAngle;
                                for (int j = 0; j < _embedDimension; j++)
                                    nearNeighbour[j] = pEmbedMatrix[index, j] - pEmbedMatrix[i + evolve, j];
                                // vector is relative to fiducial iterate
                                if ((fAngle = AngleBetweenVecs(nextNeighbour, nearNeighbour, _embedDimension)) < ANGLEMAX)
                                {
                                    wBestIndex = p;
                                    break;
                                }
                                else if (fAngle < fBestAngle)
                                {
                                    fBestAngle = fAngle;
                                    wBestIndex = p;
                                }
                            }
                        }
                    }
                    if (wBestIndex == MAXPOINTS)
                    {
                        continue; // don't calculate this one, carry on.
                    }
                    index = (int)((KDTreePoint)indexes[wBestIndex]).index;
                    neighbourDist = Math.Sqrt(((KDTreePoint)indexes[wBestIndex]).distSquared);
                    foundCount++;
                }
                // now caculate the distance between the iterates of the search pattern and the nearest neighbours.
                double fIterateDist = 0.0;
                for (int j = 0; j < _embedDimension; j++)
                {
                    double fTemp = pEmbedMatrix[(index + evolve), j] - pEmbedMatrix[(i + evolve), j];
                    fIterateDist += fTemp * fTemp;
                    nextNeighbour[j] = fTemp; //vector is relative to fiducial iterate
                }// we need this for angle calculations
                if ((fIterateDist = Math.Sqrt(fIterateDist)) == 0.0f || neighbourDist == 0.0)	// would cause instability
                    foundCount--;	  // don't include this one
                else
                {
                    double fTemp = Math.Log(fIterateDist / neighbourDist) * 1.442695040889; // log2
                    lyapunovSum += fTemp;
                }
            }
            if (foundCount == 0)
                return (m_fLyapunov = 0.0f); // should never happen
            m_fLyapunov = lyapunovSum / (double)foundCount;
            return m_fLyapunov; 	// result is the average.
        }

        internal double KaplanGlass()
        {
            int patterns;
            double range = GetLargestEvent() - GetSmallestEvent();
            if (range == 0.0f)
            {
                this.m_fKaplanGlass = 0.0;
                return m_fKaplanGlass; // no variation in the data
            }
            // embed the data apropriately
            double[,] pEmbedMatrix = EmbedForChaosMetrics(_embedDimension, _embedSeparation, out patterns);
            // normalise the data so that distance heuristics are coherent
            for (int k = 0; k < patterns; k++)
            {
                for (int p = 0; p < _embedDimension; p++)
                    pEmbedMatrix[k, p] = ScaleToNet(pEmbedMatrix[k, p], -maxChange, maxChange);
            }
            HyperCubeList cCubeList = new HyperCubeList(_embedDimension, 0.0625f);
            int wStart = 0;
            for (int k = 0; k < patterns; k++)
            {
                if (k >= 1)
                {// use the kth and k-1th vector, determine which hypercubes they visit
                    double[] startVector = new double[_embedDimension];
                    double[] kthVector = new double[_embedDimension];
                    for (int n = 0; n < _embedDimension; n++)
                    {
                        startVector[n] = pEmbedMatrix[wStart, n];
                        kthVector[n] = pEmbedMatrix[k, n];
                    }
                    if (cCubeList.AddCubes(startVector, kthVector))
                        wStart = k; // returns false if trajectory slice too short
                }
            }
            return (m_fKaplanGlass = cCubeList.EvaluateVectors());
        }

        internal void SetSampleTime(TimeSpan SampleTime)
        {
            sampleTime = SampleTime;
        }
        internal TimeSpan GetSampleTime()
        {
            return sampleTime;
        }
        internal void ChangeSampleTime(TimeSpan SampleTime)
        {
            if (sampleTime != SampleTime)
            {// if the time changes, the parameters aren't valid
                sampleTime = SampleTime;
                m_bParamSet = false;
            }
        }

        internal double GetLyapunov()
        {
            return m_fLyapunov;
        }
        internal double GetHurst()
        {
            return m_fHurst;
        }
        internal double GetKaplanGlass()
        {
            return (m_fKaplanGlass);
        }

        internal double AngleBetweenVecs(double[] fVec1, double[] fVec2, int wWidth)
        {
            double fNorm1 = 0.0;
            double fNorm2 = 0.0;
            double fResult = 0.0;
            for (int p = 0; p < wWidth; p++)
            {
                fNorm1 += fVec1[p] * fVec1[p];
                fNorm2 += fVec2[p] * fVec2[p];
                fResult += fVec1[p] * fVec2[p];
            }
            fNorm1 = Math.Sqrt(fNorm1);
            fNorm2 = Math.Sqrt(fNorm2);
            if (fNorm2 == 0.0f || fNorm1 == 0.0f)
                return 3.14159;
            fResult /= fNorm1 * fNorm2;
            fResult = Math.Min(fResult, 1.0f);
            fResult = Math.Max(fResult, -1.0f);
            return Math.Acos(fResult);
        }

        internal double ScaleFromNet(double val)
        {
            double fTemp = (val) * (maxVal - minVal) + minVal;
            if (useGranularity)
            {
                fTemp = fTemp / granularity;
                if (fTemp - Math.Floor(fTemp) > 0.5) fTemp += 1.0;
                return (double)((int)fTemp) * granularity;
            }
            return fTemp;
        }
        double ScaleToNet(double val, double fMin, double fMax)
        {
            double fTemp = fMax - fMin;
            return val / fTemp - fMin / fTemp;
        }

        internal double ScaleFromNetDifferenced(double val)
        {
            double fTemp = val * maxChange * 2.0f - maxChange;
            if (useGranularity)
            {
                fTemp = fTemp / granularity;
                if (fTemp - Math.Floor(fTemp) > 0.5f) fTemp += 1.0;
                return (double)((int)fTemp) * granularity;
            }
            return fTemp;
        }

        internal DateTime GetEarliestEvent()
        {
            if (events.Count == 0)
                return DateTime.MinValue;
            return ((Event)events[0]).GetEventTime();
        }
        internal void SetGranularity(double granularity)
        {
            this.granularity = granularity;
            useGranularity = true;

        }
        internal void ClearGranularity()
        {
            useGranularity = false;
        }
        internal DateTime GetLatestEvent()
        {
            if (events.Count == 0)
                return DateTime.MaxValue;
            return ((Event)events[events.Count - 1]).GetEventTime();
        }
        double GetLargestEvent()
        {
            return maxVal;
        }
        double GetSmallestEvent()
        {
            return minVal;
        }
        double GetMaxChange()
        {
            return maxChange;
        }
        internal void RemoveBefore(DateTime startTime)
        {
            if (!sorted)
            {
                events.Sort();
                sorted = true;
            }
            if (events.Count == 0)
                return;
            else if (((Event)events[0]).GetEventTime() >= startTime)
                return; // t is earlier than the range of data.
            else if (((Event)events[events.Count - 1]).GetEventTime() < startTime)// t is later than the range of data
            {
                events.Clear();
                return;
            }
            else // Use BinarySearch
            {
                Event tempEvent = new Event(0.0, startTime);
                int index = events.BinarySearch(tempEvent);
                if (index < 0)
                {
                    index = ~index;
                    index--;
                }
                events.RemoveRange(0, index);
            }
        }

        /// <summary>
        /// Load the fifo with a new prediction and create a stimulus pattern
        /// </summary>
        /// <param name="prediction"></param>
        /// <returns></returns>
        internal double[] UpdateFifoWithPrediction(double prediction)
        {
            fifo[fifoIndex] = prediction;
            if (fifoIndex == 0)
                fifoIndex = fifo.Length - 1; // wrap array
            else
                fifoIndex--;
            return ReadoutFifo();
        }
        /// <summary>
        /// Create a stimulus pattern from the fifo
        /// </summary>
        /// <returns>the pattern</returns>
        internal double[] ReadoutFifo()
        {
            double[] stimulus = new double[_embedDimension];
            int p = (fifoIndex + 1) % fifo.Length;
            int r = 0;
            for (int q = _embedDimension; q > 0; q--)
            {
                double Value = fifo[(p + (q - 1) * _embedSeparation) % fifo.Length];
                if (differencing)
                    stimulus[r++] = ScaleToNet(Value, -maxChange, maxChange);
                else
                    stimulus[r++] = ScaleToNet(Value, minVal, maxVal);
            }
            return stimulus;
        }
        /// <summary>
        /// The best separation is that which generates the best predictions
        /// </summary>
        private void CalculateSeparationUsingPrediction()
        {
            int dimension = _embedDimension;
            double[] errors = new double[MAXDELAY];
            double bestError = double.MaxValue;
            int bestSeparation = 0;
            //search over neighbour embeddings
            /*			int bestDimension = dimension;
                        int bestSeparation = 1;
                        double lowestError = double.MaxValue;
                        for(int p = 0; p < 3; p++)
                        {
                            if(_embedDimension == 2) 
                                dimension = _embedDimension + p;
                            else
                                dimension = _embedDimension + p - 2;
                            */

            for (int separation = 1; separation < MAXDELAY; separation++)
            {
                int patterns;
                double[,] stim;
                double[] targ;
                double[] validStim;
                stim = EmbedForPrediction(dimension, separation, out patterns, out targ, false, new TimeSpan(1, 0, 0, 0), this.GetLatestEvent(), out validStim);
                //Convert from row to Column vector, because that's what SugMay expects
                double[,] columnTarg = new double[patterns, 1];
                for (int s = 0; s < patterns; s++)
                {
                    columnTarg[s, 0] = targ[s];
                }
                SugMay sugMay = new SugMay(stim, columnTarg);
                for (int s = 0; s < patterns; s++)
                {
                    for (int n = 0; n < dimension; n++)
                        validStim[n] = stim[s, n];
                    double confidence;
                    double prediction = sugMay.Predict(validStim, 0, out confidence);
                    /*					if(IsDifferenced())
                                            prediction = ScaleFromNetDifferenced(prediction);
                                        else
                                            prediction = ScaleFromNet(prediction);*/
                    double error = targ[s] - prediction;
                    errors[separation - 1] += error * error;
                }
                errors[separation - 1] = Math.Sqrt(errors[separation - 1] / patterns);
                if (separation <= 5)
                {
                    if (errors[separation - 1] < bestError)
                    {
                        bestError = errors[separation - 1];
                        bestSeparation = separation;
                    }
                }
                else if (separation > 5 && errors[separation - 1] > bestError)
                {
                    this._embedSeparation = bestSeparation;
                    break;
                }
                else if (separation > 5 && errors[separation - 1] >= errors[separation - 2])
                {
                    this._embedSeparation = separation - 1;
                    break;
                }
            }
        }
        /// <summary>
        /// Performs linear interpolation between the values of events
        /// </summary>
        /// <param name="first">the earlier event</param>
        /// <param name="second">the later event</param>
        /// <param name="time">the point in time for which to generate an interpolated value</param>
        /// <returns>the interpolated value.</returns>
        private double Interpolate(Event first, Event second, DateTime time)
        {
            TimeSpan interval = second.GetEventTime() - first.GetEventTime();
            TimeSpan offset = time - first.GetEventTime();
            if (interval.TotalSeconds <= 0.0 || offset > interval)
                return (first.GetEventValue() + second.GetEventValue()) / 2.0;
            double fraction = offset.TotalSeconds / interval.TotalSeconds;
            return first.GetEventValue() + (second.GetEventValue() - first.GetEventValue()) * fraction;
        }
        internal bool IsInterpolated()
        {
            return interpolated;
        }
        internal void SetInterpolation(bool val)
        {
            interpolated = val;
        }

        internal void Sort()
        {
            if (!sorted)
            {
                events.Sort();
                sorted = true;
            }
        }
        #endregion



    }


	internal class FastMutualMI
	{
		#region constructor
		internal FastMutualMI(double [,] data)
		{
			s = new int[data.Length]; 
			q = new int[data.Length]; 
			pop = new double[data.Length]; 
			pindex = new int[data.Length]; 
			taumax = 20;	//local	
			// n0 = nearest power of 2 less than data size
			n0 = 1; 
			while ((n0 + taumax) <= data.Length)  
				n0 *= 2;
			n0 /= 2;
			bin = -1; 
			int j = n0;
			pow2 = new int[25];
			for (int i=0; i <25; i++)  //KMAX
			{
				pow2[i] = j;
				j /= 2;
			}
			for (int i=0; i < n0; i++)  
			{
				pop[i] = data[i,0];
				pindex[i] = i;
			}
			SwinneyAndVastano(data);
		}
		#endregion
		#region fields
		private int [] s;
		private int [] q;
		private double [] pop;
		private int [] pindex;
		private int [] pow2; 
		private int n0;
		private int bin;
		private int taumax; 
		internal int separation;
		#endregion
		#region helpers
		private void SwinneyAndVastano(double[,] data)
		{	
			double smallestInfo = double.MaxValue;
			double [] results = new double[taumax+1];
			SortElements(0,n0-1);
			for (int i=0; i<n0; i++)  
				s[i] = pindex[i]; 
			for(int tau=0; tau <= taumax; tau++)  
			{
				/* now do tau offset for q[] data */
				for (int i=tau; i<tau+n0; i++)  
				{
					pop[i-tau] = data[i,0];
					pindex[i-tau] = i-tau;
				}
				SortElements(0,n0-1);
				for (int i=0; i<n0; i++)  
					q[pindex[i]] = i;
				/* now find I(S,Q) according to formula (19) */
				int [] kmarray = new int[25];
				double x,y;

				kmarray[0] = 0;
				x = (double) n0;
				y = ffunct(kmarray,0);
				results[tau] = (1.0/x)*y - Math.Log(x)*1.4426950408889634073599246810019;
				if(results[tau] > smallestInfo)
				{
					separation = tau - 1;
					break;
				}
				else
					smallestInfo = results[tau];
			}
			
		}


		private double ffunct(int [] kmarray,int m)
		{
			double Value;
			int [] temparray = new int[25];
			for (int j=0; j <= m; j++)  
				temparray[j] = kmarray[j];
			int n = number(temparray,m);		
			Value = ((double)n);
			if (n<=1)  
			{
				Value = 0.0;
			} 
			else if (n == 2)  
			{
				Value = 4.0;
			} 
			else if (m == bin)  
			{
				/* no substructure */
				Value = Value * Math.Log(Value)*1.4426950408889634073599246810019;
			} 
			else 
			{
				/* assume substructure exists */
				Value = Value * 2.0;
				for (int j=0; j<=3; j++)  
				{
					temparray[m+1] = j;
					Value += ffunct(temparray,m+1);
				}
			}
			return Value;
		}

		private int number(int [] karray2,int m2)
		{
			int ivalue = n0;
			if (m2 > 0)  
			{
				int los = 0;
				int loq = 0;
				int his = n0; 
				int hiq = n0;
				for (int i=1; i <= m2; i++)  
				{
					if (karray2[i]%2==0)  
						his -= pow2[i];
					else                  
						los += pow2[i];
					if (karray2[i]<2)     
						hiq -= pow2[i];
					else                  
						loq += pow2[i];
				}
				ivalue = 0;
				for (int i=los; i < his; i++)  
				{
					int j = q[s[i]];
					if (j >= loq  && j < hiq)
						ivalue++;
				}
			} 
			return ivalue;
		}

		private void SortElements(int p,int r)
		{
			if (p < r) 
			{
				int q = PartitionElements(p,r);
				SortElements(p,q);
				SortElements(q+1,r);
			}
		}


		private int PartitionElements(int p,int r)
		{

			double x = pop[p];
			int i = p - 1;
			int j = r + 1;
			while(true) 
			{
				do
				{
					--j;
				}while(pop[j] < x);
				do
				{
					++i;
				}while(pop[i] > x);

				if (i < j) 
				{
					/* here's where the in-place swap takes place */
					/* fortunately pindex[] stores the *original* location */
					double temp = pop[i]; 
					pop[i] = pop[j]; 
					pop[j] = temp;
					int tempi = pindex[i]; 
					pindex[i] = pindex[j]; 
					pindex[j] = tempi;
				}
				else
					return j;
			}
		}
		#endregion
	}
}
