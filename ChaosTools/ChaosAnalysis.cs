using System;
using System.Collections.Generic;

namespace ConceptStrings.ChaosTools
{
    /// <summary>
    /// Control class for chaos tools library
    /// </summary>
    /// <remarks>The Chaos tools library is concerned with the analysis of time series data of a single variate.
    /// Data values supplied to the library are time tagged double length floating point values, and methods exist to load arrays, 
    /// single values at a fixed point in time, or values at a given time offset from previous values.
    /// An internal temporal database stores and orders the values. In order to unravel the internal dynamics of the
    /// time series embedding is used according to the method of Floris Takens, and two tools are provided to automatically determine 
    /// the optimal embedding dimension and separation. Predictions are generated 1 and 'n' steps ahead of a given point in the temporal database
    /// using a memory-based learning algorithm. Various measures of chaotic behavior are calculated, including Lyapunov exponent and Hurst Exponent. 
    /// More information can be found at: <see href="https://www.academia.edu/466858/Time_series_prediction_using_supervised_learning_and_tools_from_chaos_theory"/> (Acrobat reader required.)>
    ///  </remarks>
    [Serializable]
	public class ChaosAnalysis
    {
        #region constructor
        /// <summary>
		/// Constructor
		/// </summary>
		/// <remarks>Creates a new ChaosTools instance.</remarks>
		public ChaosAnalysis()
		{
			eventList = new EventList();
			m_bAttentionLimit = false;
		}

        #endregion

        #region fields

        internal EventList eventList;
        private bool m_bAttentionLimit;
        private TimeSpan m_ctsAttentionSpan;
        private SugMay sugMay;

        #endregion

        #region properties
        /// <summary>
		/// The sample time of the analyses
		/// </summary>
		public TimeSpan sampleTime
		{
			set
			{
				if(value == new TimeSpan(0,0,0,0,0))
					throw new Exception("sample time cannot be zero");
				eventList.SetSampleTime(value);
			}
			get
			{
				return eventList.GetSampleTime();
			}
		}

        /// <summary>
        /// The timestamp of the earliest data item in the temporal database
        /// </summary>
        /// <remarks>gives <see cref="DateTime.MinValue"/>if database is empty</remarks>
        public DateTime earliest
        {
            get
            {
                return eventList.GetEarliestEvent();
            }
        }
        /// <summary>
        /// supplies the earliest in time value.
        /// </summary>
        public double earliestValue
        {
            get
            {
                double val;
                eventList.GetValue(out val, eventList.GetEarliestEvent());
                return val;
            }
        }
        /// <summary>
        /// The timestamp of the latest data item in the temporal database
        /// </summary>
        /// <remarks>gives <see cref="DateTime.MaxValue"/>if database is empty</remarks>
        public DateTime latest
        {
            get
            {
                return eventList.GetLatestEvent();
            }
        }
        /// <summary>
        /// Supplies the farthest along in time value;
        /// </summary>
        public double latestValue
        {
            get
            {
                double val;
                eventList.GetValue(out val, eventList.GetLatestEvent());
                return val;
            }
        }

        /// <summary>
        /// The number of sampled intervals in the temporal database.
        /// </summary>
        /// <exception cref="Exception">if sample time not set</exception>
        public int Count
        {
            get
            {
                return eventList.Count;
            }
        }
        /// <summary>
        /// Access the sampled values in the temporal database by sample time index
        /// </summary>
        /// <remarks>0 = first element in the database, Count - 1 = last, and values betweenthese points are interpolated</remarks>
        public double this[int index]
        {
            get
            {
                return eventList[index];
            }
        }
        /// <summary>
        /// Access the sampled values in the temporal database by sample time
        /// </summary>
        /// <remarks>For sample times before the earliest value NaN is returned. 
        /// For sample times after the last value, the last value is returned. </remarks>
        public double this[DateTime sample]
        {
            get
            {
                double res;
                if (!eventList.GetValue(out res, sample))
                    res = double.NaN;
                return res;
            }
            set
            {
                eventList.AddEvent(value, sample);
            }
        }

        /// <summary>
        /// Sets differencing for predictions
        /// </summary>
        /// <remarks>Where data are constrained to a known range this is not required. If data are unconstrained, as in financial time series, first differences are taken of the 
        /// time series to create a difference series, and this is used in training/prediction.
        /// i.e. without differencing <code>S[0], S[1], ... S[n]</code> are used. 
        /// With differencing <code>S[1] - S[0], S[2] - S[1], ... S[n] - S[n-1]</code> are used.</remarks>
        public bool differencing
        {
            get
            {
                return eventList.IsDifferenced();
            }
            set
            {
                eventList.SetDifferencing(value);
            }
        }

        public bool MovingAverage
        {
            get
            {
                return eventList.MovingAverage;
            }
            set
            {
                eventList.MovingAverage = value;
            }

        }
        /// <summary>
        /// If true, interpolation is used.
        /// </summary>
        /// <remarks>When a requested sample is between two samples in the database the value used will be linearly interpolated between the two values if true, or the earlier value is used if false. </remarks>
        public bool interpolation
        {
            get
            {
                return eventList.IsInterpolated();
            }
            set
            {
                eventList.SetInterpolation(value);
            }
        }
        ///<summary>Time period for detecting breaks in the stored series</summary>
        /// <remarks>
        /// Some real world data has breaks in it.
        ///	ChaosKit will by default interpolate accross these breaks, 
        ///	and whether linear or spline interpolation is used, 
        ///	bad data is liable to get into the training set. 
        ///	Setting this period greater than <see cref="TimeSpan.Zero"/>
        ///	turns on dead period detection, and any gaps greater than this 
        ///	threshold will cause Chaoskit to treat the first subsequent data point
        ///	as if it were the start of a time series. 
        /// </remarks>
        public TimeSpan deadPeriodThreshold
        {
            get
            {
                if (eventList.deadPeriodDetectionEnabled)
                    return eventList.deadPeriodThreshold;
                return TimeSpan.Zero;
            }
            set
            {
                if (value != TimeSpan.Zero)
                {
                    eventList.deadPeriodDetectionEnabled = true;
                    eventList.deadPeriodThreshold = value;
                }
            }
        }

        /// <summary>
        /// Get the lyapunov exponent measure of this time series
        /// </summary>
        /// <remarks>This measures the loss of predictability inherent 
        /// in the time series in bits per sample step, and thus the 
        /// "chaoticness" of the time series. The calculation is averaged 
        /// over the time series.</remarks>
        public double Lyapunov
        {
            get
            {
                return eventList.GetLyapunov();
            }
        }

        /// <summary>
        /// Get the result of the RSS analysis
        /// </summary>
        /// <returns>The Hurst exponent</returns>
        public double Hurst
        {
            get
            {
                return eventList.GetHurst();
            }
        }

        /// <summary>
        /// Get the fractal dimension of this series.
        /// </summary>
        /// <remarks>This is not the dimension of the attractor generating the time series, but the graph of the series itself.
        /// It thus ranges from 1 for a straight line, to 2.0 for a series so active it fills the graphical space.
        /// Directly related to the Hurst exponent, it is a useful measure akin to volatility, but not dependant on scale.</remarks>
        public double FractalDimension
        {
            get
            {
                return 2.0 - eventList.GetHurst();
            }
        }

        /// <summary>
        /// Get the calculated optimal embedding separation for this time series.
        /// </summary>
        /// <remarks>The method used is based on auto-mutual information and is described in
        /// <see href="http://www.scientio.com/resources/thesis.pdf"/> page 79</remarks>
        /// <returns>the calculated optimal embedding separation</returns>
        public int optimalEmbeddingSeparation
        {
            get
            {
                return (int)eventList.embedSeparation;
            }
            set
            {
                eventList.embedSeparation = value;
            }
        }
        /// <summary>
        /// Get the calculated optimal embedding dimension for this time series.
        /// </summary>
        /// <remarks>Uses the method of "False Nearest Neighbours" to calculate the optimum embedding dimension.
        /// See <see href="http://www.scientio.com/resources/thesis.pdf"/> Section 4.3.2 page 81 for more details.</remarks>
        /// <returns>The calculated optimal embedding dimension</returns>
        public int optimalEmbeddingDimension
        {
            get
            {
                return (int)eventList.embedDimension;
            }
            set
            {
                eventList.embedDimension = value;
            }
        }
        /// <summary>
        /// The associated trading times object
        /// </summary>
        public TradingTimes tradingTimes
        {
            get
            {
                if (!eventList.useTradingTimes)
                    return null;
                return eventList.tradingTimes;
            }
            set
            {
                eventList.SetTradingTimes(value);
            }
        }

        /// <summary>
        /// Determines during prediction how much of the past history is used.
        /// </summary>
        public TimeSpan AttentionSpan
        {
            get
            {
                return m_ctsAttentionSpan;
            }
            set
            {
                if (value == TimeSpan.Zero)
                    ResetAttentionSpan();
                SetAttentionSpan(value);
            }
        }

        /// <summary>
        /// defines the minimum step of a prediction. I.e. if 0.1, predictions will be made to the nearest tenth.
        /// </summary>
        public double Granularity
        {
            set
            {
                eventList.SetGranularity(value);
            }
            get
            {
                return eventList.granularity;
            }
        }


        #endregion

        #region public methods

        /// <summary>
		/// Empty events and analyses.
		/// </summary>
		public void ClearData()
		{
			eventList.ClearAll();
			eventList.ResetMaxMin();
		}
		/// <summary>
		/// Removes the elements of the stored series before the given date and time.
		/// </summary>
		/// <param name="startDate">Start date and time of the data to be retained.</param>
		public void RemoveBefore(DateTime startDate)
		{
			eventList.RemoveBefore(startDate);
		}
		/// <summary>
		/// Load data items consisting of double, time pairs
		/// </summary>
		/// <remarks>Length of two arrays must be equal and > 0.</remarks>
		/// <param name="data">An array of data items</param>
		/// <param name="times">An array of corresponding times</param>
		/// <exception cref="Exception">thrown if lengths of arrays not equal.</exception>
		public void LoadTaggedData(double [] data, DateTime [] times)
		{
			if(data.Length != times.Length)
				throw new Exception("LoadTaggedData data and times arrays not equal lengths.");
			for(int n = 0; n < data.Length; n++)
				eventList.AddEvent(data[n], times[n]);
		}

        /// <summary>
        /// Add a data point at a fixed point in time
        /// </summary>
        /// <param name="data">data point</param>
        /// <param name="time">time of the data point measurement/event</param>
        public void AddDataPointFixed(double data, DateTime time)
        {
            eventList.AddEvent(data, time);
        }
        /// <summary>
        /// Adds a data point at a time offset from the last entered.
        /// </summary>
        /// <param name="data">the data point value</param>
        /// <param name="interval">the interval from the last.</param>
        public void AddDataPointRelative(double data, TimeSpan interval)
        {
            DateTime last = eventList.GetLatestEvent();
            last += interval;
            eventList.AddEvent(data, last);
        }

        /// <summary>
        /// Calculate the embedding parameters and the various measures.
        /// </summary>
        /// <remarks>Not all time series can be embedded. 
        /// Specifically time series that are truly random or exhibit 
        /// hyperchaos or chaotic behaviour beyond the range of search are rejected.
        /// Since the Lyapunov exponent calculation relies on the embedding process 
        /// the calculated value should not be trusted if the function returns false.</remarks>
        /// <param name="hurst">If true Hurst exponent is calculated</param>
        /// <param name="lyapunov">If true Lyapunov exponent is measured</param>
        /// <param name="kaplan">If true Kaplan-Glass measure is calculated.</param>
        /// <returns>false if embedding parameters could not be inferred.</returns>
        public bool CalculateMeasures(bool lyapunov, bool kaplan, bool hurst)
        {
            return eventList.SetupEmbedding(lyapunov, kaplan, hurst);
        }
        /// <summary>
        /// Load a series of sampled data values.
        /// </summary>
        /// <remarks>this function assumes that there are no breaks in the sample sequence.
        /// Usage:<code>
        ///double [] logistic = new double[10000];
        ///logistic[0] = 0.1;
        ///for(int n = 1; n &lt; 10000; n++)
        ///	logistic[n] = 4.0 * logistic[n-1] * (1.0 - logistic[n-1]);
        ///ChaosAnalysis analysis = new ChaosAnalysis();
        ///analysis.LoadSampledData(logistic,DateTime.Now, new TimeSpan(0,0,1));
        ///analysis.sampleTime = new TimeSpan(0,0,1);
        ///</code></remarks>
        /// <param name="data">An array of data values</param>
        /// <param name="startTime">the start of the sampling</param>
        /// <param name="sampleStep">the sampling interval</param>
        public void LoadSampledData(double[] data, DateTime startTime, TimeSpan sampleStep)
        {
            if (sampleStep == new TimeSpan(0, 0, 0, 0, 0))
                throw new Exception("LoadSampledData SampleStep cannot be zero");
            eventList.SetSampleTime(sampleStep);
            for (int n = 0; n < data.Length; n++)
            {
                eventList.AddEvent(data[n], startTime);
                startTime += sampleStep;
            }
        }

        /// <summary>
        /// Predict one or more values beyond 
        /// the most recent in the temporal database.
        /// </summary>
        /// <remarks>The predictions are "iterated". This means that 
        /// predictions are themselves used as the basis for further predictions.
        /// Predictions, as they are generated, are treated as genuine time 
        /// samples and injected into the embedding process. They are not
        /// added to the temporal database, however.
        /// This function is relatively slow because "training" - i.e. assembly and embedding of
        /// all the existing training examples - is performed before predictions are generated.
        /// if start is not within the range of valid data in the temporal database an exception is thrown.
        /// otherwise predictions are generated at 
        /// <code>
        ///start + <see cref="sampleTime"/>, 
        ///start + 2 * <see cref="sampleTime"/>
        /// </code> and so on.
        /// Example usage:
        /// <code>		
        ///double [] logistic = new double[1000];
        ///logistic[0] = 0.1;
        ///for(int n = 1; n &lt; 1000; n++)
        ///logistic[n] = 4.0 * logistic[n-1] * (1.0 - logistic[n-1]);
        ///ChaosAnalysis analysis = new ChaosAnalysis();
        ///TimeSpan sampleTime = new TimeSpan(0,0,1);
        ///analysis.LoadSampledData(logistic,DateTime.Now, sampleTime );
        ///Assert.IsTrue(analysis.CalculateMeasures(true,false,true));
        ///logistic[0] = logistic[logistic.Length - 1];
        ///double error = 0.0;
        ///for(int n = 1; n &lt; 1000; n++)
        ///{
        ///	logistic[n] = 4.0 * logistic[n-1] * (1.0 - logistic[n-1]);
        ///	double [] result = analysis.Predict(1,analysis.latest);
        ///	analysis.AddDataPointRelative(logistic[n],sampleTime);
        ///	double diff = logistic[n] - result[0];
        ///	error += diff * diff;
        ///}
        ///error = Math.Sqrt(error/1000);
        ///	</code></remarks>
        /// <param name="predictions">Number of predictions to create <see cref="sampleTime"/> apart.</param>
        /// <param name="start">The time to create predictions from</param>
        /// <exception cref="Exception">Thrown if start is outside bounds of temporal database, or choice of parameters and data make predictions impossible.</exception>
        /// <returns>The predicted values</returns>
        public double[] Predict(int predictions, DateTime start)
        {
            int patterns;
            return Predict(predictions, start, out patterns);
        }

        /// <summary>
        /// Predict one or more values beyond 
        /// the most recent in the temporal database.
        /// </summary>
        /// <remarks>The predictions are "iterated". This means that 
        /// predictions are themselves used as the basis for further predictions.
        /// Predictions, as they are generated, are treated as genuine time 
        /// samples and injected into the embedding process. They are not
        /// added to the temporal database, however.
        /// This function is relatively slow because "training" - i.e. assembly and embedding of
        /// all the existing training examples - is performed before predictions are generated.
        /// if start is not within the range of valid data in the temporal database an exception is thrown.
        /// otherwise predictions are generated at 
        /// <code>
        ///start + <see cref="sampleTime"/>, 
        ///start + 2 * <see cref="sampleTime"/>
        /// </code> and so on.
        /// Example usage:
        /// <code>		
        ///double [] logistic = new double[1000];
        ///logistic[0] = 0.1;
        ///for(int n = 1; n &lt; 1000; n++)
        ///logistic[n] = 4.0 * logistic[n-1] * (1.0 - logistic[n-1]);
        ///ChaosAnalysis analysis = new ChaosAnalysis();
        ///TimeSpan sampleTime = new TimeSpan(0,0,1);
        ///analysis.LoadSampledData(logistic,DateTime.Now, sampleTime );
        ///Assert.IsTrue(analysis.CalculateMeasures(true,false,true));
        ///logistic[0] = logistic[logistic.Length - 1];
        ///double error = 0.0;
        ///for(int n = 1; n &lt; 1000; n++)
        ///{
        ///	logistic[n] = 4.0 * logistic[n-1] * (1.0 - logistic[n-1]);
        ///	double [] result = analysis.Predict(1,analysis.latest);
        ///	analysis.AddDataPointRelative(logistic[n],sampleTime);
        ///	double diff = logistic[n] - result[0];
        ///	error += diff * diff;
        ///}
        ///error = Math.Sqrt(error/1000);
        ///	</code></remarks>
        /// <param name="predictions">Number of predictions to create <see cref="sampleTime"/> apart.</param>
        /// <param name="start">The time to create predictions from</param>
        /// <param name="patterns">The number of training paterns generated given the dimension, separation and sample time</param>
        /// <exception cref="Exception">Thrown if start is outside bounds of temporal database, or choice of parameters and data make predictions impossible.</exception>
        /// <returns>The predicted values</returns>
        public double[] Predict(int predictions, DateTime start, out int patterns)
        {
            int dimension = eventList.embedDimension;
            int separation = eventList.embedSeparation;
            if (start > latest)
                start = latest;
            if (start < earliest)
                throw new Exception("Chaos Analysis: start time in predict out of bounds of temporal database.");
            double[,] stim;
            double[] targ;
            double[] validStim;
            stim = eventList.EmbedForPrediction(dimension, separation, out patterns, out targ, m_bAttentionLimit, m_ctsAttentionSpan, start, out validStim);
            //Convert from row to Column vector, because that's what SugMay expects
            double[,] columnTarg = new double[patterns, 1];
            for (int s = 0; s < patterns; s++)
            {
                columnTarg[s, 0] = targ[s];
            }
            sugMay = new SugMay(stim, columnTarg);
            return IteratePredictions(predictions, validStim);
        }

        /// <summary>
        /// Predicts without retraining.
        /// </summary>
        /// <remarks><see cref="Predict"/> should be called before using this method. Because no training 
        /// is employed the method is faster than <see cref="Predict"/> but the internal model is not updated.
        /// Intersperse calls to this method with calls to <see cref="Predict"/> to keep the model up to date.
        /// Example usage:
        /// <code>
        ///double [] logistic = new double[1000];
        ///logistic[0] = 0.1;
        ///for(int n = 1; n &lt; 1000; n++)
        ///	logistic[n] = 4.0 * logistic[n-1] * (1.0 - logistic[n-1]);
        ///ChaosAnalysis analysis = new ChaosAnalysis();
        ///TimeSpan sampleTime = new TimeSpan(0,0,1);
        ///analysis.LoadSampledData(logistic,DateTime.Now, sampleTime );
        ///analysis.CalculateMeasures(true,false,true);
        ///logistic[0] = logistic[logistic.Length - 1];
        ///double error = 0.0;
        ///for(int n = 1; n &lt; 1000; n++)
        ///{
        ///	logistic[n] = 4.0 * logistic[n-1] * (1.0 - logistic[n-1]);
        ///	double [] result;
        ///	if(n == 1)
        ///		result = analysis.Predict(1,analysis.latest);
        ///	else
        ///		result = analysis.PredictNext(1,analysis.latest);
        ///	analysis.AddDataPointRelative(logistic[n],sampleTime);
        ///	double diff = logistic[n] - result[0];
        ///	error += diff * diff;
        ///}
        ///error = Math.Sqrt(error/1000);
        ///</code></remarks>
        /// <returns>The predicted values</returns>
        public double[] PredictNext(int predictions, DateTime start)
        {
            if (start > eventList.GetLatestEvent() || start < eventList.GetEarliestEvent())
                throw new Exception("Chaos Analysis: start time in predict out of bounds of temporal database.");
            double[] validStim;
            // fill the fifo from the start time backwards
            DateTime datum = start;
            while (!eventList.FillFifo(datum, this.sampleTime))
            {
                datum -= this.sampleTime;
            }
            validStim = eventList.ReadoutFifo();
            return IteratePredictions(predictions, validStim);
        }
        /// <summary>
        /// Helper function for iterating through making the predictions.
        /// </summary>
        /// <param name="predictions"></param>
        /// <param name="validStim"></param>
        /// <returns></returns>
        private double[] IteratePredictions(int predictions, double[] validStim)
        {
            double confidence;
            int n = 0;
            double[] result = new double[predictions];
            do
            {
                double prediction = sugMay.Predict(validStim, 0, out confidence);
                if (eventList.IsDifferenced())
                    result[n] = eventList.ScaleFromNetDifferenced(prediction);
                else
                    result[n] = eventList.ScaleFromNet(prediction);
                validStim = eventList.UpdateFifoWithPrediction(prediction);
                n++;
            } while (n < predictions);
            return result;
        }

        /// <summary>
        /// Gets the version of the data seen by the inference system.
        /// </summary>
        /// <param name="AttentionSpan">if not null the length of time preceding the end to use in processing</param>
        /// <returns>the samples created</returns>
        public List<double> GetSampledData(TimeSpan AttentionSpan)
        {
            double[,] stim;
            int patterns;
            int countToReturn = (int)(AttentionSpan.Ticks / sampleTime.Ticks);
            stim = eventList.EmbedForChaosMetrics(1, 1, out patterns);
            int patStart = countToReturn == 0 ? 0 : Math.Min(0, patterns - countToReturn);
            List<double> results = new List<double>();
            //return smaller of countToReturn or patterns points
            for (int i = 0; i < patterns; i++)
            {
                if(i >= patStart)
                    results.Add(stim[i, 0]);
            }
            return results;
        }

        /// <summary>
        /// Used to get an interpolated point
        /// </summary>
        /// <param name="point">Time of point</param>
        /// <returns>the value</returns>
        public double GetSampledPoint(DateTime point)
        {
            double res = 0.0;
            if (eventList.GetValue(out res, point))
                return res;
            return 0.0;
        }

        /// <summary>
        /// Sorts the stored events
        /// </summary>
        public void SortEvents()
        {
            eventList.Sort();
        }

        /// <summary>
        /// For predictions made off the end of the time series finds next valid times
        /// </summary>
        /// <param name="predictions">number of steps to predict</param>
        /// <remarks>Where trading times are used non-trading times are removed from the series. Predictions are made in this "collapsed time", 
        /// but when the time stamps of predictions are converted back into real time the times need to be processed to put non-trading times back in. 
        /// So for instance, if predicting one day ahead for traded shares on a Friday, the time of the next prediction should be Monday.</remarks>
        /// <returns>Array of times</returns>
        public DateTime[] GetNextValidTradingTimes(int predictions)
        {
            List<DateTime> times = new List<DateTime>();
            DateTime start = this.latest;
            for (int n = 0; n < predictions; n++)
            {
                do
                {
                    start += this.sampleTime;
                }
                while (this.tradingTimes != null ? !this.tradingTimes.IsTradingTime(start, this.sampleTime) : false);
                times.Add(start);
            }
            return times.ToArray();
        }

        /// <summary>
        /// For back testing find next valid trading times
        /// </summary>
        /// <param name="predictions">number of steps to predict</param>
        /// <param name="time">Time of start of predictions</param>
        /// <remarks>Where trading times are used non-trading times are removed from the series. Predictions are made in this "collapsed time", 
        /// but when the time stamps of predictions are converted back into real time the times need to be processed to put non-trading times back in. 
        /// So for instance, if predicting one day ahead for traded shares on a Friday, the time of the next prediction should be Monday.</remarks>
        /// <returns>Array of times</returns>
        public DateTime[] GetNextValidTradingTimes(int predictions, DateTime time)
        {
            List<DateTime> times = new List<DateTime>();
            DateTime start = time;
            for (int n = 0; n < predictions; n++)
            {
                do
                {
                    start += this.sampleTime;
                }
                while (this.tradingTimes != null ? !this.tradingTimes.IsTradingTime(start, this.sampleTime) : false);
                times.Add(start);
            }
            return times.ToArray();
        }


        #endregion

        #region internal methods
        /// <summary>
		/// If financial data used, set up the periods for which it is valid
		/// </summary>
		/// <param name="trading"></param>
		internal void SetTradingTimes(TradingTimes trading)
		{
			eventList.SetTradingTimes(trading);
		}

        /// <summary>
        /// Get the result of the RSS analysis
        /// </summary>
        /// <returns>The Hurst exponent</returns>
        internal double GetHurstExponent()
        {
            return eventList.GetHurst();
        }
        /// <summary>
        /// Get the Kaplan Glass statistic for this time series.
        /// </summary>
        /// <returns></returns>
        internal double GetKaplanGlassMeasure()
        {
            return eventList.GetKaplanGlass();
        }
        /// <summary>
        /// Get the fractal dimension of this series.
        /// </summary>
        /// <remarks>This is not the dimension of the attractor generating the time series, but the graph of the series itself.
        /// It thus ranges from 1 for a straight line, to 2.0 for a series so active it fills the graphical space.
        /// Directly related to the Hurst exponent, it is a useful measure akin to volatility, but not dependant on scale.</remarks>
        /// <returns>Number between 1.0 and 2.0</returns>
        internal double GetFractalDimension()
        {
            return 2.0 - eventList.GetHurst();
        }

        /// <summary>
        /// Get the lyapunov exponent measure of this time series
        /// </summary>
        /// <remarks>This measures the loss of predictability inherent 
        /// in the time series in bits per sample step, and thus the 
        /// "chaoticness" of the time series. The calculation is averaged 
        /// over the time series.</remarks>
        /// <returns>The calculated exponent</returns>
        internal double GetLyapunovExponent()
        {
            return eventList.GetLyapunov();
        }
 

        /// <summary>
        /// Control differencing
        /// </summary>
        /// <remarks>If differencing is true, predictions will be based on the 
        /// differences between succeding samples, rather than their absolute values.
        /// This is most usefull with financial series, where the absolute trading value
        /// is of less interest than the changes (which are the source of profit), and where 
        /// price levels are not often repeated.</remarks>
        /// <param name="bDifference">true for differencing on</param>
        internal void SetDifferencing(bool bDifference)
        {
            eventList.SetDifferencing(bDifference);
        }

        /// <summary>
        /// Sets the length of history used in analysis
        /// </summary>
        /// <remarks>For many social and economic time series external factors can change the implied generating function of the series, if any exists.
        /// These of course occur unpredictably, but some defense can be gained by ignoring data that is older than a given period.
        /// By setting an attention span, the prediction will ignore data items that might be predictive that are older than the attention span. </remarks>
        /// <param name="period"></param>
        internal void SetAttentionSpan(TimeSpan period)
        {
            if (period > new TimeSpan(0, 0, 0, 0, 0))
                m_bAttentionLimit = true;
            m_ctsAttentionSpan = period;
        }
        /// <summary>
        /// Stop using only a limited history.
        /// </summary>
        internal void ResetAttentionSpan()
        {
            m_bAttentionLimit = false;
        }

        /// <summary>Set the granularity of predictions.</summary>
        /// <remarks>
        /// Many financial time series are quantised, 
        /// It makes little sense to predict values to arbitrary resolution 
        /// when they must be constrained to some set of quantised values. 
        /// Setting the prediction granularity to the natural quantisation 
        /// of the time series, for instance 0.0001 for the $/£ exchange rate,
        /// will result in the predictions being rounded to the nearest quantum.
        /// </remarks>
        /// <param name="predGran">The size of the smallest possible increment</param>
        internal void SetPredictionGranularity(double predGran)
        {
            eventList.SetGranularity(predGran);
        }


        #endregion

        #region static methods

        /// <summary>
        /// Calculates the mutual information in bits between two vectors
        /// </summary>
        public static double MutualInformation(double[] x, double[] y)
        {
            return EventList.MutualInformation(x, y);
        }

        #endregion
 
    }
}
