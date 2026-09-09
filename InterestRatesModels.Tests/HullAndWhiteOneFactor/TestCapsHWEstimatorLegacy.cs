using System.Collections.Generic;
using DVPLI;
using NUnit.Framework;

namespace HullAndWhiteOneFactor
{
    [TestFixture]
    public class TestCapsHWEstimatorLegacy
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void ProvidesToReturnsStocasticProcessHW()
        {
            CapsHWEstimatorLegacy estimator = new CapsHWEstimatorLegacy();

            Assert.AreEqual(typeof(DVPLDOM.StocasticProcessHW), estimator.ProvidesTo);
        }

        [Test]
        public void EstimateDummyCalibrationIsInheritedFromBase()
        {
            CapsHWEstimatorLegacy estimator = new CapsHWEstimatorLegacy();
            InterestRateMarketData dataset = TestCommon.TestMarketDataFactory.CreateCapMarketData(
                new Vector(new double[] { 1, 2 }),
                new Matrix(new double[,] { { 0.20, 0.22 }, { 0.21, 0.23 } }));
            Fairmat.Calibration.CapVolatilityFiltering settings = new Fairmat.Calibration.CapVolatilityFiltering
            {
                DummyCalibration = true
            };

            EstimationResult result = estimator.Estimate(new List<object> { dataset }, settings);

            CollectionAssert.AreEqual(new string[] { "Alpha", "Sigma" }, result.Names);
            CollectionAssert.AreEqual(new double[] { 0.1, 0.05 }, result.Values);
        }
    }
}
