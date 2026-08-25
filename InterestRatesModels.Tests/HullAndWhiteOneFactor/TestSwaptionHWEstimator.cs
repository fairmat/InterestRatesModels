using System;
using System.Collections.Generic;
using DVPLI;
using NUnit.Framework;

namespace HullAndWhiteOneFactor
{
    [TestFixture]
    public class TestSwaptionHWEstimator
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
                OptionMaturity = new Vector(new double[] { 1, 2 }),
                SwapDuration = new Vector(new double[] { 2, 5 }),
                SwaptionTenor = 0.5,
                SwaptionsVolatility = new Matrix(new double[,] { { 0.20, 0.22 }, { 0.21, 0.23 } })
            };
        }

        [Test]
        public void DescriptionsAndProvidesToAreFixed()
        {
            SwaptionHWEstimator estimator = new SwaptionHWEstimator();

            Assert.AreEqual("Calibrate Using Swaption", estimator.ToolTipText);
            Assert.AreEqual("Calibrate HW1 using Swaption", estimator.Description);
            Assert.AreEqual(typeof(HW1), estimator.ProvidesTo);
        }

        [Test]
        public void GetRequirementsReturnsInterestRateMarketData()
        {
            SwaptionHWEstimator estimator = new SwaptionHWEstimator();

            EstimateRequirement[] requirements = estimator.GetRequirements(null, null);

            Assert.AreEqual(1, requirements.Length);
            Assert.AreEqual(typeof(InterestRateMarketData), requirements[0].MarketDataType);
        }

        [Test]
        public void EstimateDummyCalibrationReturnsFixedGuess()
        {
            SwaptionHWEstimator estimator = new SwaptionHWEstimator();
            InterestRateMarketData dataset = CreateMarketData();
            Fairmat.Calibration.SwaptionsFiltering settings = new Fairmat.Calibration.SwaptionsFiltering
            {
                DummyCalibration = true
            };

            EstimationResult result = estimator.Estimate(new List<object> { dataset }, settings);

            CollectionAssert.AreEqual(new string[] { "Alpha", "Sigma" }, result.Names);
            CollectionAssert.AreEqual(new double[] { 0.1, 0.1 }, result.Values);
            CollectionAssert.AreEqual((double[])dataset.ZRMarketDates.ToArray(), result.ZRX);
            CollectionAssert.AreEqual((double[])dataset.ZRMarket.ToArray(), result.ZRY);
        }

        [Test]
        public void EstimateReturnsErrorWhenFilteringExcludesAllSwaptions()
        {
            SwaptionHWEstimator estimator = new SwaptionHWEstimator();
            InterestRateMarketData dataset = CreateMarketData();
            Fairmat.Calibration.SwaptionsFiltering settings = new Fairmat.Calibration.SwaptionsFiltering
            {
                MinSwaptionMaturity = 100,
                MaxSwaptionMaturity = 200,
                MinSwapDuration = 0,
                MaxSwapDuration = 30
            };

            EstimationResult result = estimator.Estimate(new List<object> { dataset }, settings);

            Assert.AreEqual("No swaptions satisfying criteria found, please relax filters", result.ErrorMessage);
        }
    }
}
