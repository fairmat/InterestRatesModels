using System;
using System.Collections.Generic;
using DVPLI;
using Fairmat.Finance;
using NUnit.Framework;

namespace Pelsser.Calibration
{
    [TestFixture]
    public class TestCapletRealWordEstimator
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void DescriptionsAndProvidesToAreFixed()
        {
            CapletRealWordEstimator estimator = new CapletRealWordEstimator();

            Assert.AreEqual("Real world calibration using caps and zero rate curve historical series", estimator.ToolTipText);
            Assert.AreEqual("Real world calibration using caps and zero rate curve historical series", estimator.Description);
            Assert.AreEqual(typeof(SquaredGaussianModel), estimator.ProvidesTo);
            Assert.IsInstanceOf<LambdaCalibrationSettings>(estimator.DefaultSettings);
        }

        [Test]
        public void GetRequirementsReturnsMarketDataAndDiscountingCurves()
        {
            CapletRealWordEstimator estimator = new CapletRealWordEstimator();

            EstimateRequirement[] requirements = estimator.GetRequirements(null, null);

            Assert.AreEqual(2, requirements.Length);
            Assert.AreEqual(typeof(InterestRateMarketData), requirements[0].MarketDataType);
            Assert.AreEqual(typeof(Fairmat.MarketData.DiscountingCurveMarketData[]), requirements[1].MarketDataType);
        }

        [Test]
        public void EstimateThrowsWhenSettingsIsNull()
        {
            CapletRealWordEstimator estimator = new CapletRealWordEstimator();
            InterestRateMarketData dataset = new InterestRateMarketData();

            Assert.Throws<NullReferenceException>(() =>
                estimator.Estimate(new List<object> { dataset, new Fairmat.MarketData.DiscountingCurveMarketData[0] }, null));
        }

        [Test]
        public void EstimateThrowsWhenYearsExceedsBondMaturity()
        {
            CapletRealWordEstimator estimator = new CapletRealWordEstimator();
            InterestRateMarketData dataset = new InterestRateMarketData();
            LambdaCalibrationSettings settings = new LambdaCalibrationSettings
            {
                Years = 10,
                BondMaturity = 5
            };

            Exception ex = Assert.Throws<Exception>(() =>
                estimator.Estimate(new List<object> { dataset, new Fairmat.MarketData.DiscountingCurveMarketData[0] }, settings));

            Assert.AreEqual("Bond maturity has to be greater of the historical series time span.", ex.Message);
        }
    }
}
