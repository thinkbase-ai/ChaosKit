using System;
using System.Collections.Generic;
using System.Diagnostics;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using ConceptStrings.ChaosTools;

namespace ChaosToolsTest
{
    [TestClass]
    public class ChaosAnalysisTest
    {
        static private int samples = 513;

        public ChaosAnalysisTest()
        {
            //
            // TODO: Add constructor logic here
            //
        }

        private TestContext testContextInstance;

        /// <summary>
        ///Gets or sets the test context which provides
        ///information about and functionality for the current test run.
        ///</summary>
        public TestContext TestContext
        {
            get
            {
                return testContextInstance;
            }
            set
            {
                testContextInstance = value;
            }
        }

        #region Additional test attributes
        //
        // You can use the following additional attributes as you write your tests:
        //
        // Use ClassInitialize to run code before running the first test in the class
        // [ClassInitialize()]
        // public static void MyClassInitialize(TestContext testContext) { }
        //
        // Use ClassCleanup to run code after all tests in a class have run
        // [ClassCleanup()]
        // public static void MyClassCleanup() { }
        //
        // Use TestInitialize to run code before running each test 
        // [TestInitialize()]
        // public void MyTestInitialize() { }
        //
        // Use TestCleanup to run code after each test has run
        // [TestCleanup()]
        // public void MyTestCleanup() { }
        //
        #endregion


        [TestMethod]
        public void TestLogistic()
        {
            Trace.WriteLine("logistic");
            double[] logistic = new double[10000];
            logistic[0] = 0.1;
            for (int n = 1; n < 10000; n++)
                logistic[n] = 4.0 * logistic[n - 1] * (1.0 - logistic[n - 1]);
            ChaosAnalysis analysis = new ChaosAnalysis();
            analysis.LoadSampledData(logistic, DateTime.Now, new TimeSpan(0, 0, 1));
            analysis.sampleTime = new TimeSpan(0, 0, 1);
            Assert.IsTrue(analysis.CalculateMeasures(true, false, true));
            Assert.AreEqual(2, analysis.optimalEmbeddingDimension);
            Assert.AreEqual(1, analysis.optimalEmbeddingSeparation);
            Assert.AreEqual(0.55934, analysis.Hurst, 0.001);
            Assert.AreEqual(0.98410, analysis.Lyapunov, 0.001);
        }
        /// <summary>
        /// Test using Tent Map
        /// </summary>
        [TestMethod]
        public void TestTent()
        {
            Trace.WriteLine("Tent");
            double[] tent = new double[samples];
            tent[0] = 0.1;
            for (int n = 1; n < samples; n++)
            {
                if (tent[n - 1] < 0.5)
                    tent[n] = 1.9 * tent[n - 1];
                else
                    tent[n] = 1.9 * (1 - tent[n - 1]);
            }
            ChaosAnalysis analysis = new ChaosAnalysis();
            analysis.LoadSampledData(tent, DateTime.Now, new TimeSpan(0, 0, 1));
            analysis.sampleTime = new TimeSpan(0, 0, 1);
            Assert.IsTrue(analysis.CalculateMeasures(true, false, true));
            Assert.AreEqual(2, analysis.optimalEmbeddingDimension);
            Assert.AreEqual(1, analysis.optimalEmbeddingSeparation);
            Assert.AreEqual(0.4570, analysis.Hurst, 0.001);
            Assert.AreEqual(0.9018, analysis.Lyapunov, 0.001);
        }

        /// <summary>
        /// Test using Sine
        /// </summary>
        [TestMethod]
        public void TestSine()
        {
            Trace.WriteLine("Sine");
            double[] sine = new double[samples];
            for (int n = 0; n < samples; n++)
            {
                sine[n] = Math.Sin(n * Math.PI / 180);
            }
            ChaosAnalysis analysis = new ChaosAnalysis();
            analysis.LoadSampledData(sine, DateTime.Now, new TimeSpan(0, 0, 1));
            analysis.sampleTime = new TimeSpan(0, 0, 1);
            Assert.IsTrue(analysis.CalculateMeasures(true, false, true));
            Assert.AreEqual(2, analysis.optimalEmbeddingDimension);
            Assert.AreEqual(2, analysis.optimalEmbeddingSeparation);
            Assert.AreEqual(0.997, analysis.Hurst, 0.001);
            Assert.AreEqual(0.384, analysis.Lyapunov, 0.001);
        }

        [TestMethod]
        public void TestPartCrazy()
        {
            Trace.WriteLine("PartCrazy");
            double[] partCrazy = new double[samples];
            for (int n = 0; n < samples; n++)
            {
                partCrazy[n] = (n % 3) / 2;
            }
            ChaosAnalysis analysis = new ChaosAnalysis();
            analysis.LoadSampledData(partCrazy, DateTime.Now, new TimeSpan(0, 0, 1));
            analysis.sampleTime = new TimeSpan(0, 0, 1);
            Assert.IsTrue(analysis.CalculateMeasures(true, false, true));
            Assert.AreEqual(2, analysis.optimalEmbeddingDimension);
            Assert.AreEqual(1, analysis.optimalEmbeddingSeparation);
            Assert.AreEqual(0.0217, analysis.Hurst, 0.001);
            Assert.AreEqual(0.5, analysis.Lyapunov, 0.001);
        }


        /// <summary>
        /// Test with a straight line. Hurst should be 1.0, Lyapunov 0.0 since no information is lost.
        /// </summary>
        [TestMethod]
        public void TestLinear()
        {
            Trace.WriteLine("Linear");
            double[] linear = new double[samples];
            linear[0] = 0.1;
            for (int n = 1; n < samples; n++)
            {
                linear[n] = linear[n - 1] + 1.0;
            }
            ChaosAnalysis analysis = new ChaosAnalysis();
            analysis.LoadSampledData(linear, DateTime.Now, new TimeSpan(0, 0, 1));
            analysis.sampleTime = new TimeSpan(0, 0, 1);
            Assert.IsTrue(analysis.CalculateMeasures(true, false, true));
            Assert.AreEqual(2, analysis.optimalEmbeddingDimension);
            Assert.AreEqual(1, analysis.optimalEmbeddingSeparation);
            Assert.AreEqual(1.000, analysis.Hurst, 0.001);
            Assert.AreEqual(0.000, analysis.Lyapunov, 0.001);
        }
        /// <summary>
        /// Test with continually reverting waveform. Fractal dimension should be ~ 2.0.
        /// Hurst should be 0.0;
        /// </summary>
        [TestMethod]
        public void TestCrazy()
        {
            Trace.WriteLine("Crazy");
            double[] linear = new double[samples];
            linear[0] = 0.1;
            for (int n = 1; n < samples; n++)
            {
                if (n % 2 == 0)
                    linear[n] = 0.1;
                else
                    linear[n] = 0.5;
            }
            ChaosAnalysis analysis = new ChaosAnalysis();
            analysis.LoadSampledData(linear, DateTime.Now, new TimeSpan(0, 0, 1));
            analysis.sampleTime = new TimeSpan(0, 0, 1);
            Assert.IsTrue(analysis.CalculateMeasures(true, false, true));
            Assert.AreEqual(0.0000, analysis.Hurst, 0.001);
            Assert.AreEqual(0.000, analysis.Lyapunov, 0.001);
        }
        /// <summary>
        /// Test with the Mackay Glass equation
        /// </summary>
        [TestMethod]
        public void TestMackayGlass()
        {
            Trace.WriteLine("MackayGlass");
            double[] mackayGlass = new double[10100];
            for (int n = 0; n < 18; n++)
                mackayGlass[n] = 0.7;
            for (int n = 19; n < 10100; n++)
            {
                mackayGlass[n] = 0.9 * mackayGlass[n - 1] + 0.2 * (mackayGlass[n - 18] / (1.0 + Math.Pow(mackayGlass[n - 18], 10)));
            }
            double[] mackayGlass2 = new double[10000];
            for (int n = 100; n < 10100; n++)
            {
                mackayGlass2[n - 100] = mackayGlass[n];
            }
            ChaosAnalysis analysis = new ChaosAnalysis();
            analysis.LoadSampledData(mackayGlass2, DateTime.Now, new TimeSpan(0, 0, 1));
            analysis.sampleTime = new TimeSpan(0, 0, 1);
            Assert.IsTrue(analysis.CalculateMeasures(true, false, true));
            Assert.AreEqual(2, analysis.optimalEmbeddingDimension);
            Assert.AreEqual(1, analysis.optimalEmbeddingSeparation);
            Assert.AreEqual(0.65412, analysis.Hurst, 0.005);
            Assert.AreEqual(1.07493, analysis.Lyapunov, 0.01);
        }
        /// <summary>
        /// Test with random numbers
        /// </summary>
        [TestMethod]
        public void TestRandom()
        {
            Trace.WriteLine("Random");
            double[] random = new double[samples];
            Random rand = new Random(1234);
            for (int n = 0; n < samples; n++)
            {
                random[n] = rand.NextDouble();
            }
            ChaosAnalysis analysis = new ChaosAnalysis();
            analysis.LoadSampledData(random, DateTime.Now, new TimeSpan(0, 0, 1));
            analysis.sampleTime = new TimeSpan(0, 0, 1);
            Assert.IsFalse(analysis.CalculateMeasures(true, false, true));
            Assert.AreEqual(2, analysis.optimalEmbeddingDimension);
            Assert.AreEqual(2, analysis.optimalEmbeddingSeparation);
            Assert.AreEqual(0.6283, analysis.Hurst, 0.001);
            Assert.AreEqual(2.94264, analysis.Lyapunov, 0.01);
        }
        /*	/// <summary>
            /// Test the MI calcs with random data
            /// </summary>
            public void TestMIRandom()
            {
                double [] random1 = new double[samples];
                double [] random2 = new double[samples];
                Random rand = new Random(3567);
                for(int n = 0; n < samples; n++)
                {
                    random1[n] = rand.NextDouble();
                    random2[n] = rand.NextDouble();
                }
                Assert.AreEqual(0.0000,ChaosAnalysis.MutualInformation(random1,random2),0.001);
                for(int n = 0; n < samples; n++)
                {
                    random1[n] = (double)n;
                    random2[n] = (double)n;
                }
                Assert.AreEqual(0.0000,ChaosAnalysis.MutualInformation(random1,random2),0.001);
            }*/
        /// <summary>
        /// Test the search by distance method
        /// </summary>
        [TestMethod]
        public void TestKDTreeSearchByDistance()
        {
            double[,] testVectors = new double[50, 3];
            Random rand = new Random(3567);
            bool finished = false;
            int n = 0;
            while (!finished)
            {
                double x = rand.NextDouble();
                double y = rand.NextDouble();
                double z = rand.NextDouble();
                //create vector, reject if not within 0.25 of 0,0,0
                if (x * x + y * y + z * z < 0.25)
                {
                    testVectors[n, 0] = x;
                    testVectors[n, 1] = y;
                    testVectors[n, 2] = z;
                    n++;
                    if (n == 25)
                        finished = true;
                }
            }
            finished = false;
            while (!finished)
            {
                double x = rand.NextDouble();
                double y = rand.NextDouble();
                double z = rand.NextDouble();
                //create vector, reject if not within 0.25 of 0,0,0
                if (x * x + y * y + z * z > 0.25)
                {
                    testVectors[n, 0] = x;
                    testVectors[n, 1] = y;
                    testVectors[n, 2] = z;
                    n++;
                    if (n == 50)
                        finished = true;
                }
            }
            //now search by distance should return the first 25 values.
            KDTree tree = new KDTree(testVectors);
            double[] zero = new double[3];
            zero[0] = 0.0;
            zero[1] = 0.0;
            zero[2] = 0.0;
            int count;
            List<KDTreePoint> distances = new List<KDTreePoint>();
            tree.SearchByDistance(zero, 0.5, out count, ref distances, 30);
            Assert.AreEqual(25, count);
            Assert.IsTrue(((KDTreePoint)distances[count - 1]).distSquared < 0.25);
            zero[0] = testVectors[0, 0];
            zero[1] = testVectors[0, 1];
            zero[2] = testVectors[0, 2];
            tree.SearchByDistance(zero, 0.5, out count, ref distances, 30);
            // check that the target, which is also part of the dataset is the first
            Assert.AreEqual(0, ((KDTreePoint)distances[0]).index);
            Assert.AreEqual(0.0, ((KDTreePoint)distances[0]).distSquared);
        }
        /// <summary>
        /// Test the search buy distance method
        /// </summary>
        [TestMethod]
        public void TestKDTreeSearchByDistanceObject()
        {
            TreeReferencableObject[,] testVectors = new TreeReferencableObject[50, 3];
            Random rand = new Random(3567);
            bool finished = false;
            int n = 0;
            while (!finished)
            {
                double x = rand.NextDouble();
                double y = rand.NextDouble();
                double z = rand.NextDouble();
                //create vector, reject if not within 0.25 of 0,0,0
                if (x * x + y * y + z * z < 0.25)
                {
                    testVectors[n, 0] = new TreeReferencableObject(x);
                    testVectors[n, 1] = new TreeReferencableObject(y);
                    testVectors[n, 2] = new TreeReferencableObject(z);
                    n++;
                    if (n == 25)
                        finished = true;
                }
            }
            finished = false;
            while (!finished)
            {
                double x = rand.NextDouble();
                double y = rand.NextDouble();
                double z = rand.NextDouble();
                //create vector, reject if not within 0.25 of 0,0,0
                if (x * x + y * y + z * z > 0.25)
                {
                    testVectors[n, 0] = new TreeReferencableObject(x);
                    testVectors[n, 1] = new TreeReferencableObject(y);
                    testVectors[n, 2] = new TreeReferencableObject(z);
                    n++;
                    if (n == 50)
                        finished = true;
                }
            }
            //now search by distance should return the first 25 values.
            KDTree tree = new KDTree(testVectors);
            double[] zero = new double[3];
            zero[0] = 0.0;
            zero[1] = 0.0;
            zero[2] = 0.0;
            int count;
            List<KDTreePoint> distances = new List<KDTreePoint>();
            tree.SearchByDistance(zero, 0.5, out count, ref distances, 30);
            Assert.AreEqual(25, count);
            Assert.IsTrue(((KDTreePoint)distances[count - 1]).distSquared < 0.25);
            zero[0] = testVectors[0, 0].GetValue();
            zero[1] = testVectors[0, 1].GetValue();
            zero[2] = testVectors[0, 2].GetValue();
            tree.SearchByDistance(zero, 0.5, out count, ref distances, 30);
            // check that the target, which is also part of the dataset is the first
            Assert.AreEqual(0, ((KDTreePoint)distances[0]).index);
            Assert.AreEqual(0.0, ((KDTreePoint)distances[0]).distSquared);
        }
        /// <summary>
        /// Test the Search method
        /// </summary>
        [TestMethod]
        public void TestKDTreeSearch()
        {
            double[,] testVectors = new double[50, 3];
            Random rand = new Random(3567);
            for (int n = 0; n < 50; n++)
            {
                testVectors[n, 0] = rand.NextDouble();
                testVectors[n, 1] = rand.NextDouble();
                testVectors[n, 2] = rand.NextDouble();
            }
            KDTree tree = new KDTree(testVectors);
            double[] zero = new double[3];
            zero[0] = 0.0;
            zero[1] = 0.0;
            zero[2] = 0.0;
            List<KDTreePoint> distances = new List<KDTreePoint>();
            tree.Search(zero, 5, ref distances);
            Assert.AreEqual(5, distances.Count);
            zero[0] = 0.6;
            zero[1] = 0.6;
            zero[2] = 0.6;
            tree.Search(zero, 5, ref distances);
            Assert.AreEqual(5, distances.Count);
            zero[0] = testVectors[0, 0];
            zero[1] = testVectors[0, 1];
            zero[2] = testVectors[0, 2];
            tree.Search(zero, 5, ref distances);
            // check that the target, which is also part of the dataset is the first
            Assert.AreEqual(0, ((KDTreePoint)distances[0]).index);
            Assert.AreEqual(0.0, ((KDTreePoint)distances[0]).distSquared);
        }
        /// <summary>
        /// Test the Search method for objects
        /// </summary>
        [TestMethod]
        public void TestKDTreeSearchObject()
        {
            TreeReferencableObject[,] testVectors = new TreeReferencableObject[50, 3];
            Random rand = new Random(3567);
            for (int n = 0; n < 50; n++)
            {
                testVectors[n, 0] = new TreeReferencableObject(rand.NextDouble());
                testVectors[n, 1] = new TreeReferencableObject(rand.NextDouble());
                testVectors[n, 2] = new TreeReferencableObject(rand.NextDouble());
            }
            KDTree tree = new KDTree(testVectors);
            TreeReferencableObject[] zero = new TreeReferencableObject[3];
            zero[0] = new TreeReferencableObject(0.0);
            zero[1] = new TreeReferencableObject(0.0);
            zero[2] = new TreeReferencableObject(0.0);
            List<KDTreePoint> distances = new List<KDTreePoint>();
            tree.Search(zero, 5, ref distances);
            Assert.AreEqual(5, distances.Count);
            zero[0] = new TreeReferencableObject(0.6);
            zero[1] = new TreeReferencableObject(0.6);
            zero[2] = new TreeReferencableObject(0.6);
            tree.Search(zero, 5, ref distances);
            Assert.AreEqual(5, distances.Count);
            zero[0] = testVectors[0, 0];
            zero[1] = testVectors[0, 1];
            zero[2] = testVectors[0, 2];
            tree.Search(zero, 5, ref distances);
            // check that the target, which is also part of the dataset is the first
            Assert.AreEqual(0, ((KDTreePoint)distances[0]).index);
            Assert.AreEqual(0.0, ((KDTreePoint)distances[0]).distSquared);
        }
        /// <summary>
        /// Test the Unique Search facility
        /// </summary>
        [TestMethod]
        public void TestKDTreeUniqueSearch()
        {
            double[,] testVectors = new double[50, 3];
            Random rand = new Random(3567);
            for (int n = 0; n < 25; n++)
            {
                testVectors[n, 0] = rand.NextDouble();
                testVectors[n, 1] = rand.NextDouble();
                testVectors[n, 2] = rand.NextDouble();
            }
            //second part of array is duplicate of first
            for (int n = 25; n < 50; n++)
            {
                testVectors[n, 0] = testVectors[n - 25, 0];
                testVectors[n, 1] = testVectors[n - 25, 1];
                testVectors[n, 2] = testVectors[n - 25, 2];
            }
            KDTree tree = new KDTree(testVectors);
            double[] zero = new double[3];
            zero[0] = 0.0;
            zero[1] = 0.0;
            zero[2] = 0.0;
            List<KDTreePoint> distances = new List<KDTreePoint>();
            tree.UniqueSearch(zero, 5, ref distances);
            Assert.AreEqual(5, distances.Count);
            zero[0] = testVectors[0, 0];
            zero[1] = testVectors[0, 0];
            zero[2] = testVectors[0, 0];
            tree.UniqueSearch(zero, 5, ref distances);
            Assert.AreEqual(5, distances.Count);
            Assert.IsFalse((int)((KDTreePoint)distances[0]).index == 0);
            Assert.IsFalse(((KDTreePoint)distances[0]).distSquared == 0.0);
            Assert.IsFalse(((KDTreePoint)distances[0]).index == ((KDTreePoint)distances[1]).index);
            Assert.IsFalse(((KDTreePoint)distances[1]).index == ((KDTreePoint)distances[2]).index);
            Assert.IsFalse(((KDTreePoint)distances[2]).index == ((KDTreePoint)distances[3]).index);
            Assert.IsFalse(((KDTreePoint)distances[3]).index == ((KDTreePoint)distances[4]).index);
            zero[0] = 0.6;
            zero[1] = 0.6;
            zero[2] = 0.6;
            tree.UniqueSearch(zero, 5, ref distances);
            Assert.AreEqual(5, distances.Count);
        }
        /// <summary>
        /// Test the Unique Search facility
        /// </summary>
        [TestMethod]
        public void TestKDTreeUniqueSearchObject()
        {
            TreeReferencableObject[,] testVectors = new TreeReferencableObject[50, 3];
            Random rand = new Random(3567);
            for (int n = 0; n < 25; n++)
            {
                testVectors[n, 0] = new TreeReferencableObject(rand.NextDouble());
                testVectors[n, 1] = new TreeReferencableObject(rand.NextDouble());
                testVectors[n, 2] = new TreeReferencableObject(rand.NextDouble());
            }
            //second part of array is duplicate of first
            for (int n = 25; n < 50; n++)
            {
                testVectors[n, 0] = testVectors[n - 25, 0];
                testVectors[n, 1] = testVectors[n - 25, 1];
                testVectors[n, 2] = testVectors[n - 25, 2];
            }
            KDTree tree = new KDTree(testVectors);
            double[] zero = new double[3];
            zero[0] = 0.0;
            zero[1] = 0.0;
            zero[2] = 0.0;
            List<KDTreePoint> distances = new List<KDTreePoint>();
            tree.UniqueSearch(zero, 5, ref distances);
            Assert.AreEqual(5, distances.Count);
            zero[0] = testVectors[0, 0].GetValue();
            zero[1] = testVectors[0, 0].GetValue();
            zero[2] = testVectors[0, 0].GetValue();
            tree.UniqueSearch(zero, 5, ref distances);
            Assert.AreEqual(5, distances.Count);
            Assert.IsFalse((int)((KDTreePoint)distances[0]).index == 0);
            Assert.IsFalse(((KDTreePoint)distances[0]).distSquared == 0.0);
            Assert.IsFalse(((KDTreePoint)distances[0]).index == ((KDTreePoint)distances[1]).index);
            Assert.IsFalse(((KDTreePoint)distances[1]).index == ((KDTreePoint)distances[2]).index);
            Assert.IsFalse(((KDTreePoint)distances[2]).index == ((KDTreePoint)distances[3]).index);
            Assert.IsFalse(((KDTreePoint)distances[3]).index == ((KDTreePoint)distances[4]).index);
            zero[0] = 0.6;
            zero[1] = 0.6;
            zero[2] = 0.6;
            tree.UniqueSearch(zero, 5, ref distances);
            Assert.AreEqual(5, distances.Count);
        }
        /*		/// <summary>
                /// 
                /// </summary>
                [TestMethod]
                public void TestFastMutualMILogistic()
                {
                    double [,] logistic = new double[samples,1];
                    logistic[0,0] = 0.1;
                    for(int n = 1; n < samples; n++)
                        logistic[n,0] = 4.0 * logistic[n-1,0] * (1.0 - logistic[n-1,0]);
                    FastMutualMI fast = new FastMutualMI(logistic);
                }*/
        /// <summary>
        /// 
        /// </summary>
        /*		public void TestFastMutualMIMackayGlass()
                {
                    double [,] mackayGlass = new double[samples,1];
                    for(int n = 0; n < 18; n++)
                        mackayGlass[n,0] = 0.7;
                    for(int n = 19; n < samples; n++)
                    {
                        mackayGlass[n,0] = 0.9 * mackayGlass[n-1,0] + 0.2 * (mackayGlass[n-18,0]/(1.0 + Math.Pow(mackayGlass[n-18,0],10)));
                    }
                    FastMutualMI fast = new FastMutualMI(mackayGlass);
                }*/
        /// <summary>
        /// 
        /// </summary>
        [TestMethod]
        public void TestRossler()
        {
            Trace.WriteLine("Rossler");
            double h = 0.05; /* or smaller */
            double a = 0.2;
            double b = 0.2;
            double c = 5.7;
            XYZ p = new XYZ(0, 0, 0); ;
            XYZ plast = new XYZ(0.1, 0, 0);
            double[] xDimRossler = new double[10000];
            int index = 0;
            for (int i = 0; i < 10100; i++)
            {
                p.x = plast.x + h * (-plast.y - plast.z);
                p.y = plast.y + h * (plast.x + a * plast.y);
                p.z = plast.z + h * (b + plast.z * (plast.x - c));
                if (i > 100)
                {
                    xDimRossler[index++] = p.x;
                }
                plast = p;
            }
            ChaosAnalysis analysis = new ChaosAnalysis();
            analysis.LoadSampledData(xDimRossler, DateTime.Now, new TimeSpan(0, 0, 1));
            analysis.sampleTime = new TimeSpan(0, 0, 1);
            Assert.IsTrue(analysis.CalculateMeasures(true, false, true));
            Assert.AreEqual(2, analysis.optimalEmbeddingDimension);
            Assert.AreEqual(1, analysis.optimalEmbeddingSeparation);
            Assert.AreEqual(0.82081, analysis.Hurst, 0.001);
            Assert.AreEqual(0.838, analysis.Lyapunov, 0.001);

        }
        /// <summary>
        /// Tests using the Lorenz Attractor
        /// </summary>
        [TestMethod]
        public void TestLorentz()
        {
            Trace.WriteLine("Lorenz");
            double h = 0.01;
            double a = 10.0;
            double b = 28.0;
            double c = 8.0 / 3.0;
            XYZ p = new XYZ(0, 0, 0); ;
            XYZ plast = new XYZ(0.1, 0, 0);
            double[] xDimLorentz = new double[10000];
            int index = 0;
            for (int i = 0; i < 10100; i++)
            {
                p.x = plast.x + h * a * (plast.y - plast.x);
                p.y = plast.y + h * (plast.x * (b - plast.z) - plast.y);
                p.z = plast.z + h * (plast.x * plast.y - c * plast.z);
                if (i > 100)
                {
                    xDimLorentz[index++] = p.x;
                }
                plast = p;
            }
            ChaosAnalysis analysis = new ChaosAnalysis();
            analysis.LoadSampledData(xDimLorentz, DateTime.Now, new TimeSpan(0, 0, 1));
            analysis.sampleTime = new TimeSpan(0, 0, 1);
            Assert.IsTrue(analysis.CalculateMeasures(true, false, true));
            Assert.AreEqual(2, analysis.optimalEmbeddingDimension);
            Assert.AreEqual(4, analysis.optimalEmbeddingSeparation);
            Assert.AreEqual(0.90129, analysis.Hurst, 0.001);
            Assert.AreEqual(0.0789, analysis.Lyapunov, 0.001);
        }
        /// <summary>
        /// Test the predict function on the logistic series.
        /// </summary>
        [TestMethod]
        public void TestPredictLogistic()
        {
            double[] logistic = new double[1000];
            logistic[0] = 0.1;
            for (int n = 1; n < 1000; n++)
                logistic[n] = 4.0 * logistic[n - 1] * (1.0 - logistic[n - 1]);
            ChaosAnalysis analysis = new ChaosAnalysis();
            TimeSpan sampleTime = new TimeSpan(0, 0, 1);
            analysis.LoadSampledData(logistic, DateTime.Now, sampleTime);
            Assert.IsTrue(analysis.CalculateMeasures(true, false, true));
            logistic[0] = logistic[logistic.Length - 1];
            double error = 0.0;
            for (int n = 1; n < 1000; n++)
            {
                logistic[n] = 4.0 * logistic[n - 1] * (1.0 - logistic[n - 1]);
                double[] result = analysis.Predict(1, analysis.latest);
                analysis.AddDataPointRelative(logistic[n], sampleTime);
                double diff = logistic[n] - result[0];
                error += diff * diff;
            }
            error = Math.Sqrt(error / 1000);
            Assert.AreEqual(0.0017366, error, 0.001);
        }
        /// <summary>
        /// Test the predict function on the logistic series.
        /// This time use Predict for the first prediction, and PredictNext
        /// thereafter to save time - this runs 4 times faster - accuracy drops off 
        /// compared to the previous test with only Predict because 
        /// new values from the series do not become part of the training data set.
        /// </summary>
        [TestMethod]
        public void TestPredictLogisticWithPredictNext()
        {
            double[] logistic = new double[1000];
            logistic[0] = 0.1;
            for (int n = 1; n < 1000; n++)
                logistic[n] = 4.0 * logistic[n - 1] * (1.0 - logistic[n - 1]);
            ChaosAnalysis analysis = new ChaosAnalysis();
            TimeSpan sampleTime = new TimeSpan(0, 0, 1);
            analysis.LoadSampledData(logistic, DateTime.Now, sampleTime);
            Assert.IsTrue(analysis.CalculateMeasures(true, false, true));
            logistic[0] = logistic[logistic.Length - 1];
            double error = 0.0;
            for (int n = 1; n < 1000; n++)
            {
                logistic[n] = 4.0 * logistic[n - 1] * (1.0 - logistic[n - 1]);
                double[] result;
                if (n == 1)
                    result = analysis.Predict(1, analysis.latest);
                else
                    result = analysis.PredictNext(1, analysis.latest);
                analysis.AddDataPointRelative(logistic[n], sampleTime);
                double diff = logistic[n] - result[0];
                error += diff * diff;
            }
            error = Math.Sqrt(error / 1000);
            Assert.AreEqual(0.00403954, error, 0.001);
        }
        /// <summary>
        /// Test the predict function on the logistic series
        /// with multiple predictions.
        /// </summary>
        /// <remarks>This illustrates the rapid drop off in 
        /// accuracy of iterated predictions to be expected 
        /// with chaotic series.</remarks>
        [TestMethod]
        public void TestPredictLogisticMultiple()
        {
            double[] logistic = new double[1000];
            logistic[0] = 0.1;
            for (int n = 1; n < 1000; n++)
                logistic[n] = 4.0 * logistic[n - 1] * (1.0 - logistic[n - 1]);
            ChaosAnalysis analysis = new ChaosAnalysis();
            TimeSpan sampleTime = new TimeSpan(0, 0, 1);
            analysis.LoadSampledData(logistic, DateTime.Now, sampleTime);
            Assert.IsTrue(analysis.CalculateMeasures(true, false, true));
            logistic[0] = logistic[logistic.Length - 1];
            double error = 0.0;
            double[] result = analysis.Predict(10, analysis.latest);
            double[] errors = new double[10];
            for (int n = 1; n < 11; n++)
            {
                logistic[n] = 4.0 * logistic[n - 1] * (1.0 - logistic[n - 1]);
                double diff = logistic[n] - result[n - 1];
                error += diff * diff;
                errors[n - 1] = Math.Sqrt(diff * diff);
            }
            error = Math.Sqrt(error / 10);
            Assert.AreEqual(0.00631, errors[0], 0.001);
            Assert.AreEqual(0.01522, errors[1], 0.001);
            Assert.AreEqual(0.01838, errors[2], 0.001);
            Assert.AreEqual(0.05689, errors[3], 0.001);
            Assert.AreEqual(0.08324, errors[4], 0.001);
            Assert.AreEqual(0.29876, error, 0.001);
        }
        /// <summary>
        /// Test the predict function on a linear series. 
        /// </summary>
        [TestMethod]
        public void TestPredictLinear()
        {
            double[] linear = new double[1000];
            linear[0] = 0.1;
            for (int n = 1; n < 1000; n++)
            {
                linear[n] = linear[n - 1] + 1.0;
            }
            ChaosAnalysis analysis = new ChaosAnalysis();
            TimeSpan sampleTime = new TimeSpan(0, 0, 1);
            analysis.LoadSampledData(linear, DateTime.Now, sampleTime);
            Assert.IsTrue(analysis.CalculateMeasures(true, false, true));
            linear[0] = linear[linear.Length - 1];
            double error = 0.0;
            for (int n = 1; n < 1000; n++)
            {
                linear[n] = linear[n - 1] + 1.0;
                double[] result = analysis.Predict(1, analysis.latest);
                analysis.AddDataPointRelative(linear[n], sampleTime);
                double diff = linear[n] - result[0];
                error += diff * diff;
            }
            error = Math.Sqrt(error / 1000);
            Assert.AreEqual(1.5848, error, 0.01);
        }
        /// <summary>
        /// Test with a linear series differenced
        /// </summary>
        /// <remarks> This series continually extends its range. 
        /// D1fferencing produces more accurate results.
        /// This is also an interesting test because all the differences are the same.
        /// It thus represents an extreme set of cirumstances that must be handled.</remarks>
        public void TestPredictLinearDifferenced()
        {
            double[] linear = new double[1000];
            linear[0] = 0.1;
            for (int n = 1; n < 1000; n++)
            {
                linear[n] = linear[n - 1] + 1.0;
            }
            ChaosAnalysis analysis = new ChaosAnalysis();
            analysis.differencing = true;
            TimeSpan sampleTime = new TimeSpan(0, 0, 1);
            analysis.LoadSampledData(linear, DateTime.Now, sampleTime);
            Assert.IsTrue(analysis.CalculateMeasures(true, false, true));
            linear[0] = linear[linear.Length - 1];
            double error = 0.0;
            for (int n = 1; n < 1000; n++)
            {
                linear[n] = linear[n - 1] + 1.0;
                double[] result = analysis.Predict(1, analysis.latest);
                analysis.AddDataPointRelative(linear[n], sampleTime);
                //we're predicting the difference, so relative to last
                double diff = linear[n] - (result[0] + linear[n - 1]);
                error += diff * diff;
            }
            error = Math.Sqrt(error / 1000);
            Assert.AreEqual(0.0, error, 0.01);
        }
        /// <summary>
        /// Tests the temporal database
        /// </summary>
        [TestMethod]
        public void TestTemporalDatabase()
        {
            double[] logistic = new double[1000];
            logistic[0] = 0.1;
            for (int n = 1; n < 1000; n++)
                logistic[n] = 4.0 * logistic[n - 1] * (1.0 - logistic[n - 1]);
            ChaosAnalysis analysis = new ChaosAnalysis();
            DateTime start = new DateTime(1997, 1, 1, 0, 0, 0, 0);
            analysis.sampleTime = new TimeSpan(0, 0, 1);
            analysis.LoadSampledData(logistic, start, new TimeSpan(0, 0, 1));
            Assert.AreEqual(start, analysis.earliest);
            Assert.AreEqual(start + new TimeSpan(9990000000), analysis.latest);
            Assert.AreEqual(1000, analysis.Count);
            Assert.AreEqual(logistic[0], analysis[0]);
            Assert.AreEqual(logistic[999], analysis[999]);
            Assert.AreEqual(logistic[599], analysis[599]);
            Assert.AreEqual(logistic[768], analysis[768]);
            //remove before start
            analysis.RemoveBefore(new DateTime(1997, 1, 1, 0, 0, 0, 0));
            Assert.AreEqual(start, analysis.earliest);
            Assert.AreEqual(start + new TimeSpan(9990000000), analysis.latest);
            //remove before fixed point
            DateTime newStart = new DateTime(1997, 1, 1, 0, 12, 0, 0);
            analysis.RemoveBefore(newStart);
            Assert.AreEqual(newStart, analysis.earliest);
            Assert.AreEqual(start + new TimeSpan(9990000000), analysis.latest);
            newStart = new DateTime(1997, 1, 1, 1, 12, 0, 0);
            analysis.RemoveBefore(newStart);
            Assert.AreEqual(DateTime.MinValue, analysis.earliest);
            Assert.AreEqual(DateTime.MaxValue, analysis.latest);
        }
        [TestMethod]
        public void TestTemporalDatabaseWithTaggedData()
        {
            double[] logistic = new double[1000];
            DateTime[] times = new DateTime[1000];

            DateTime start = new DateTime(1997, 1, 1, 0, 0, 0, 0);
            logistic[0] = 0.1;
            times[0] = start;
            TimeSpan offset = new TimeSpan(0, 1, 0);
            DateTime iterator = start;
            for (int n = 1; n < 1000; n++)
            {
                logistic[n] = 4.0 * logistic[n - 1] * (1.0 - logistic[n - 1]);
                iterator += offset;
                times[n] = iterator;
            }
            ChaosAnalysis analysis = new ChaosAnalysis();
            analysis.LoadTaggedData(logistic, times);
            analysis.sampleTime = offset;
            analysis.interpolation = true;
            Assert.AreEqual(start, analysis.earliest);
            Assert.AreEqual(iterator, analysis.latest);
            Assert.AreEqual(1000, analysis.Count);
            Assert.AreEqual(logistic[0], analysis[0]);
            Assert.AreEqual(logistic[999], analysis[999]);
            Assert.AreEqual(logistic[599], analysis[599]);
            Assert.AreEqual(logistic[768], analysis[768]);
            //test interpolation
            TimeSpan testOfset = new TimeSpan((long)((double)offset.Ticks * 599.5));
            Assert.AreEqual((logistic[599] + logistic[600]) / 2.0, analysis[start + testOfset]);
            testOfset = new TimeSpan((long)((double)offset.Ticks * 200.5));
            Assert.AreEqual((logistic[200] + logistic[201]) / 2.0, analysis[start + testOfset]);
            //			Assert.ReferenceEquals(double.NaN,analysis[DateTime.MinValue]);
            Assert.AreEqual(logistic[999], analysis[DateTime.MaxValue]);

            //remove before start
            analysis.RemoveBefore(new DateTime(1997, 1, 1, 0, 0, 0, 0));
            Assert.AreEqual(start, analysis.earliest);
            Assert.AreEqual(iterator, analysis.latest);
            //remove before fixed point
            DateTime newStart = new DateTime(1997, 1, 1, 0, 12, 0, 0);
            analysis.RemoveBefore(newStart);
            Assert.AreEqual(newStart, analysis.earliest);
            Assert.AreEqual(iterator, analysis.latest);
            newStart = new DateTime(1997, 1, 1, 17, 0, 0, 0);
            analysis.RemoveBefore(newStart);
            Assert.AreEqual(DateTime.MinValue, analysis.earliest);
            Assert.AreEqual(DateTime.MaxValue, analysis.latest);

            //now test unequally spaced events 
            double[] data = new double[500];
            times = new DateTime[500];
            iterator = start;
            for (int n = 0; n < 500; n++)
            {
                data[n] = (double)n;
                times[n] = iterator;
                TimeSpan testOffset = new TimeSpan((long)((double)offset.Ticks * ((double)n / 100 + 1)));
                iterator += testOffset;
            }
            analysis = new ChaosAnalysis();
            analysis.LoadTaggedData(data, times);
            analysis.sampleTime = offset;
            analysis.interpolation = true;
            Assert.AreEqual(start, analysis.earliest);
            Assert.AreEqual(1742, analysis.Count);//This is not the count of events, but the count of the number of samples in the database, based on the sample time.
            //			Assert.AreEqual(data[0],analysis[times[0]]);
            Assert.AreEqual(data[199], analysis[times[199]]);
            Assert.AreEqual(data[299], analysis[times[299]]);
            Assert.AreEqual(data[468], analysis[times[468]]);
            //test interpolation
            DateTime testtime = new DateTime((times[199].Ticks + times[200].Ticks) / 2);//halfway between two samples
            Assert.AreEqual((data[199] + data[200]) / 2.0, analysis[testtime]);
            testtime = new DateTime((times[303].Ticks + times[304].Ticks) / 2);//halfway between two samples
            Assert.AreEqual((data[303] + data[304]) / 2.0, analysis[testtime], 0.000001);

        }

        [TestMethod]
        public void TestGetSampledPoint()
        {
            double[] linear = new double[1000];
            linear[0] = 0.1;
            for (int n = 1; n < 1000; n++)
            {
                linear[n] = linear[n - 1] + 1.0;
            }
            ChaosAnalysis analysis = new ChaosAnalysis();
            TimeSpan sampleTime = new TimeSpan(1, 0, 0);
            DateTime start = new DateTime(2000, 1, 1);
            analysis.LoadSampledData(linear, start, sampleTime);
            analysis.CalculateMeasures(true, false, true);
            DateTime offset1 = start + new TimeSpan(sampleTime.Ticks * 3 + sampleTime.Ticks /2 );
            Assert.AreEqual(3.1,analysis.GetSampledPoint(offset1));
            analysis.interpolation = true;
            Assert.AreEqual(3.6, analysis.GetSampledPoint(offset1),0.01);
            DateTime offset2 = start - new TimeSpan(sampleTime.Ticks * 3 + sampleTime.Ticks / 2);
            Assert.AreEqual(0.0, analysis.GetSampledPoint(offset2), 0.01);
            DateTime offset3 = analysis.latest + new TimeSpan(sampleTime.Ticks * 3 + sampleTime.Ticks / 2);
            Assert.AreEqual(999.1, analysis.GetSampledPoint(offset3), 0.01);
        }

        [TestMethod]
        public void TestGetNextValidTradingTimes()
        {
            ChaosAnalysis analysis = new ChaosAnalysis();
            TimeSpan sampleTime = new TimeSpan(1, 0, 0);
            TradingTimes tr = new TradingTimes();
            tr.tradingDays = new bool[] { false, true, true, true, true, true, false };
            double[] linear = new double[1000];
            linear[0] = 0.1;
            for (int n = 1; n < 1000; n++)
            {
                linear[n] = linear[n - 1] + 1.0;
            }
            DateTime start = new DateTime(2000, 1, 1);
            analysis.LoadSampledData(linear, start, sampleTime);
            analysis.SortEvents();
            DateTime[] offsets = analysis.GetNextValidTradingTimes(9);
            Assert.AreEqual(9, offsets.Length);

        }




    }
    internal class XYZ
    {
        internal XYZ(double _x, double _y, double _z)
        {
            x = _x;
            y = _y;
            z = _z;
        }
        internal double x;
        internal double y;
        internal double z;
    }


    internal class TreeReferencableObject : ITreeReferenceable
    {
        internal TreeReferencableObject()
        {
            internalValue = 0.0;
        }
        internal TreeReferencableObject(double val)
        {
            internalValue = val;
        }
        #region ITreeReferenceable Members

        public double GetValue()
        {
            return internalValue;
        }

        #endregion
        internal double internalValue;
    }
}
