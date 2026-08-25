using System;
using System.Collections.Generic;
using DVPLI;
using NUnit.Framework;

namespace HullAndWhiteOneFactor
{
    [TestFixture]
    public class TestCapsHWEstimator
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        private static InterestRateMarketData CreateMarketData()
        {
            return new InterestRateMarketData
            {
                ZRMarketDates = new Vector(new double[] { 0, 1, 2, 5, 10 }),
                ZRMarket = new Vector(new double[] { 0.01, 0.015, 0.017, 0.02, 0.025 }),
                CapMaturity = new Vector(new double[] { 1, 2 }),
                CapRate = new Vector(new double[] { 0.01, 0.02 }),
                CapTenor = 0.5,
                CapVolatility = new Matrix(new double[,] { { 0.20, 0.22 }, { 0.21, 0.23 } })
            };
        }

        [Test]
        public void DescriptionsAndProvidesToAreFixed()
        {
            CapsHWEstimator estimator = new CapsHWEstimator();

            Assert.AreEqual("Calibrate Using Caps", estimator.ToolTipText);
            Assert.AreEqual("Calibrate HW1 using Caps", estimator.Description);
            Assert.AreEqual(typeof(HW1), estimator.ProvidesTo);
            Assert.IsInstanceOf<Fairmat.Calibration.CapVolatilityFiltering>(estimator.DefaultSettings);
        }

        [Test]
        public void GetRequirementsReturnsInterestRateMarketData()
        {
            CapsHWEstimator estimator = new CapsHWEstimator();

            EstimateRequirement[] requirements = estimator.GetRequirements(null, null);

            Assert.AreEqual(1, requirements.Length);
            Assert.AreEqual(typeof(InterestRateMarketData), requirements[0].MarketDataType);
        }

        [Test]
        public void EstimateReturnsErrorWhenCapVolatilityIsMissing()
        {
            CapsHWEstimator estimator = new CapsHWEstimator();
            InterestRateMarketData dataset = CreateMarketData();
            dataset.CapVolatility = null;

            EstimationResult result = estimator.Estimate(new List<object> { dataset });

            Assert.AreEqual("Cap not available at requested date", result.ErrorMessage);
        }

        [Test]
        public void EstimateDummyCalibrationReturnsFixedGuess()
        {
            CapsHWEstimator estimator = new CapsHWEstimator();
            InterestRateMarketData dataset = CreateMarketData();
            Fairmat.Calibration.CapVolatilityFiltering settings = new Fairmat.Calibration.CapVolatilityFiltering
            {
                DummyCalibration = true
            };

            EstimationResult result = estimator.Estimate(new List<object> { dataset }, settings);

            CollectionAssert.AreEqual(new string[] { "Alpha", "Sigma" }, result.Names);
            CollectionAssert.AreEqual(new double[] { 0.1, 0.05 }, result.Values);
            CollectionAssert.AreEqual((double[])dataset.ZRMarketDates.ToArray(), result.ZRX);
            CollectionAssert.AreEqual((double[])dataset.ZRMarket.ToArray(), result.ZRY);
        }

        [Test]
        public void EstimateDummyCalibrationThrowsWhenZRMarketDatesIsMissing()
        {
            // Unlike CapletEstimator, CapsHWEstimator does not guard the dummy
            // branch on ZRMarket being present, so it hits an unguarded null
            // dereference here.
            CapsHWEstimator estimator = new CapsHWEstimator();
            InterestRateMarketData dataset = CreateMarketData();
            dataset.ZRMarketDates = null;
            Fairmat.Calibration.CapVolatilityFiltering settings = new Fairmat.Calibration.CapVolatilityFiltering
            {
                DummyCalibration = true
            };

            Assert.Throws<NullReferenceException>(() =>
                estimator.Estimate(new List<object> { dataset }, settings));
        }

        [Test]
        public void EstimateThrowsWhenBlackCapsAreMalformed()
        {
            CapsHWEstimator estimator = new CapsHWEstimator();
            InterestRateMarketData dataset = CreateMarketData();
            // A negative strike makes the Black formula take the log of a
            // negative number, yielding NaN and triggering the guard.
            dataset.CapRate = new Vector(new double[] { -0.01, 0.02 });

            Exception ex = Assert.Throws<Exception>(() =>
                estimator.Estimate(new List<object> { dataset }));

            Assert.AreEqual("Malformed black caps", ex.Message);
        }
    }

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
            InterestRateMarketData dataset = new InterestRateMarketData
            {
                ZRMarketDates = new Vector(new double[] { 0, 1, 2, 5, 10 }),
                ZRMarket = new Vector(new double[] { 0.01, 0.015, 0.017, 0.02, 0.025 }),
                CapMaturity = new Vector(new double[] { 1, 2 }),
                CapRate = new Vector(new double[] { 0.01, 0.02 }),
                CapTenor = 0.5,
                CapVolatility = new Matrix(new double[,] { { 0.20, 0.22 }, { 0.21, 0.23 } })
            };
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
