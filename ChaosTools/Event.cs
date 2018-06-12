using System;

namespace ConceptStrings.ChaosTools
{
	/// <summary>
	/// Represents a single time tagged data value
	/// </summary>
	[Serializable]
	internal class Event : IComparable
	{
        #region constructors
        public Event()
		{
   			dataValue = 0.0;
    		timeTag = DateTime.Now;
		}
        public Event(double Value, DateTime t)
		{
			dataValue = Value;
			timeTag = t;
		}

        public Event(double Value)
		{
			dataValue = Value;
			timeTag = DateTime.Now;
		}

        public Event(Event a)
		{
			timeTag = a.timeTag;
			dataValue = a.dataValue;    
		}

        #endregion

        #region properties

        public double dataValue { get; set; }
        public DateTime timeTag { get; set; }

        #endregion

        #region methods
        internal double GetEventValue()
		{ 
			return dataValue;
		}
		internal DateTime GetEventTime()
		{ 
			return timeTag;
		}
		internal void SetEventValue(double Value)
		{
			dataValue = Value;
		}
		public override string ToString()
		{
			return "Value: " + dataValue.ToString() + " Time: " + timeTag.ToString();
		}

        #endregion

        #region IComparable Members

        public int CompareTo(object obj)
		{
			if(!(obj is Event))
				throw new Exception("Event.CompareTo passed non Event object");
			Event evnt = obj as Event;
			return timeTag.CompareTo(evnt.timeTag);
		}

		#endregion
	}
}
