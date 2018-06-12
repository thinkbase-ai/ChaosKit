using System;

namespace ConceptStrings.ChaosTools
{
    /// <summary>
    /// Represents periods when a given market or markets are open
    /// </summary>
    [Serializable]
	public class TradingTimes
	{
		/// <summary>
		/// constructor
		/// </summary>
		public TradingTimes()
		{
			twentyFourHrs = false;
			tradingDays = new bool[7];
			tradingDays[0] = false; //0 = Sunday
			tradingDays[1] = true;
			tradingDays[2] = true;
			tradingDays[3] = true;
			tradingDays[4] = true;
			tradingDays[5] = true;
			tradingDays[6] = false;
			tradingStart = new TimeSpan(0,9,0,0,0);
			tradingEnd = new TimeSpan(0,17,0,0,0);
			holidays  = null;
		}
		/// <summary>
		/// true if 24 hour trading
		/// </summary>
		public bool twentyFourHrs;
		/// <summary>
		/// flags for trading days, [0] = Sunday, default Mon-Fri.
		/// </summary>
		public bool [] tradingDays; 
		/// <summary>
		/// trading start time, offset from midnight
		/// </summary>
		public TimeSpan tradingStart;
		/// <summary>
		/// trading end, offset from midnight
		/// </summary>
		public TimeSpan tradingEnd;
		/// <summary>
		/// list of holiday dates sorted in time order.
		/// </summary>
		public DateTime [] holidays;
		/// <summary>
		/// Determine if a given time is a trading time
		/// </summary>
		/// <param name="time">The time to examine</param>
		/// <param name="sampleTime">the current sample time</param>
		/// <returns>true if a trading time, false otherwise</returns>
		public  bool IsTradingTime(DateTime time, TimeSpan sampleTime)
		{
			if (!tradingDays[(int)time.DayOfWeek])    // invalid day
				return false;

            if (holidays != null)
            {
                foreach (DateTime day in holidays)
                {
                    if (day.Date == time.Date)//could be optimized
                        return false;
                }
            }
		
			if(sampleTime	>= new TimeSpan(1,0,0,0)) //ignore times if daily or greater predictions
				return true;
			TimeSpan offsetTime = new TimeSpan(time.Hour,time.Minute,time.Second);
			if (tradingStart < tradingEnd) // normal orientation
			{
				if((tradingStart <= offsetTime)&& (offsetTime <= tradingEnd))
					return true;
			}
			else // reverse orientation, trading start is later in the day than trading end
			{
				if(!((tradingEnd < offsetTime) && (offsetTime < tradingStart)))
					return true;
			}
			return false;
		}

	}
}
