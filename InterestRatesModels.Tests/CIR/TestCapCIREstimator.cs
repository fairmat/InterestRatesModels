using System;
using System.Collections.Generic;
using DVPLI;
using NUnit.Framework;

namespace CIRProcess
{
    [TestFixture]
    public class TestCapCIREstimator
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        private static InterestRateMarketData CreateMarketData()
        {
            return TestCommon.TestMarketDataFactory.CreateCapMarketData(
                new Vector(new double[] { 1, 2, 5 }),
                new Matrix(new double[,]
                {
                    { 0.20, 0.22 },
                    { 0.21, 0.23 },
                    { 0.19, 0.20 },
                }));
        }

        [Test]
        public void DescriptionsAndProvidesToAreFixed()
        {
            CapCIREstimator estimator = new CapCIREstimator();

            Assert.AreEqual("Calibrate using Caps", estimator.ToolTipText);
            Assert.AreEqual("Calibrate from cap matrix", estimator.Description);
            Assert.AreEqual(typeof(CIR), estimator.ProvidesTo);
        }

        [Test]
        public void GetRequirementsReturnsInterestRateMarketData()
        {
            CapCIREstimator estimator = new CapCIREstimator();

            EstimateRequirement[] requirements = estimator.GetRequirements(null, null);

            Assert.AreEqual(1, requirements.Length);
            Assert.AreEqual(typeof(InterestRateMarketData), requirements[0].MarketDataType);
        }

        [Test]
        public void EstimateDummyCalibrationReturnsFixedInitialGuessPlusR0()
        {
            CapCIREstimator estimator = new CapCIREstimator();
            InterestRateMarketData dataset = CreateMarketData();
            Fairmat.Calibration.CapVolatilityFiltering settings = new Fairmat.Calibration.CapVolatilityFiltering
            {
                DummyCalibration = true
            };

            EstimationResult result = estimator.Estimate(new List<object> { dataset }, settings);

            CollectionAssert.AreEqual(CIR.parameterNames, result.Names);
            Assert.AreEqual(4, result.Values.Length);
            Assert.AreEqual(1.0, result.Values[0], 1e-12);
            Assert.AreEqual(0.01, result.Values[1], 1e-12);
            Assert.AreEqual(0.05, result.Values[2], 1e-12);
            Assert.AreEqual(dataset.ZRMarket[0], result.Values[3], 1e-12);
        }

        [Test]
        public void EstimateThrowsWhenDataIsNotInterestRateMarketData()
        {
            CapCIREstimator estimator = new CapCIREstimator();

            Assert.Throws<NullReferenceException>(() =>
                estimator.Estimate(new List<object> { "not market data" }));
        }
    }
}
