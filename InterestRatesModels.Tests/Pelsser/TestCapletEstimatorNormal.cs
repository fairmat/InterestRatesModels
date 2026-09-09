using System;
using System.Collections.Generic;
using System.Reflection;
using DVPLDOM;
using DVPLI;
using Fairmat.Finance;
using NUnit.Framework;

namespace Pelsser.Calibration
{
    [TestFixture]
    public class TestCapletEstimatorNormal
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void DescriptionAndToolTipOverridesAreFixed()
        {
            CappletEstimatorNormal estimator = new CappletEstimatorNormal();

            Assert.AreEqual("Calibrate against Caplet prices (Normal Volatilities)", estimator.Description);
            Assert.AreEqual("Caplet Estimator (Normal Volatilities)", estimator.ToolTipText);
        }

        [Test]
        public void BlackModelFactoryReturnsBachelierNormalModel()
        {
            CappletEstimatorNormal estimator = new CappletEstimatorNormal();
            Function zr = new PFunction(null);
            zr.Expr = new double[,] { { 0, 0.01 }, { 10, 0.02 } };

            MethodInfo factory = typeof(CapletEstimator).GetMethod("BlackModelFactory",
                BindingFlags.NonPublic | BindingFlags.Instance);
            object result = factory.Invoke(estimator, new object[] { zr });

            Assert.IsInstanceOf<BachelierNormalModel>(result);
        }

        [Test]
        public void GetRequirementsReturnsTwoRequirementsWithVolAtmNField()
        {
            CappletEstimatorNormal estimator = new CappletEstimatorNormal();
            EstimateQuery query = new EstimateQuery { Market = "EU" };

            EstimateRequirement[] requirements = estimator.GetRequirements(null, query);

            Assert.AreEqual(2, requirements.Length);
            Assert.AreEqual(typeof(InterestRateMarketData), requirements[0].MarketDataType);
            Assert.AreEqual("Vol-ATM-N", requirements[1].Field);
            Assert.AreEqual(typeof(Fairmat.MarketData.MatrixMarketData), requirements[1].MarketDataType);
            Assert.AreEqual("CapVolEUR", requirements[1].TickerReplacement);
        }

        [Test]
        public void EstimateDummyCalibrationIsInheritedFromBase()
        {
            CappletEstimatorNormal estimator = new CappletEstimatorNormal();
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

            CollectionAssert.AreEqual(new string[] { "alpha1", "sigma1" }, result.Names);
            CollectionAssert.AreEqual(new double[] { 0.014, 0.001 }, result.Values);
        }
    }
}
