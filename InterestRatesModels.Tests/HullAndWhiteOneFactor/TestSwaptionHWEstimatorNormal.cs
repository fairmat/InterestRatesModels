using System;
using System.Collections.Generic;
using DVPLDOM;
using DVPLI;
using Fairmat.Calibration;
using Fairmat.Finance;
using Fairmat.MarketData;
using NUnit.Framework;

namespace HullAndWhiteOneFactor
{
    /// <summary>
    /// Exposes the protected BlackModelFactory of <see cref="SwaptionHWEstimatorNormal"/> for testing.
    /// </summary>
    internal class TestableSwaptionHWEstimatorNormal : SwaptionHWEstimatorNormal
    {
        public BlackModel CallBlackModelFactory(Function zr)
        {
            return this.BlackModelFactory(zr);
        }
    }

    /// <summary>
    /// Tests for <see cref="SwaptionHWEstimatorNormal"/> and <see cref="SwaptionHWEstimatorNormalLegacy"/>:
    /// the members they override on top of the inherited <see cref="SwaptionHWEstimator"/>, plus the
    /// inherited Estimate() dummy-calibration short-circuit (deterministic, no optimizer involved).
    /// </summary>
    [TestFixture]
    public class TestSwaptionHWEstimatorNormal
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void Description_ReturnsNormalVariantText()
        {
            SwaptionHWEstimatorNormal estimator = new SwaptionHWEstimatorNormal();

            Assert.That(estimator.Description, Is.EqualTo("Calibrate HW1 using Swaption [Normal]"));
        }

        [Test]
        public void GetRequirements_ReturnsInterestRateMarketDataAndNormalVolRequirements()
        {
            SwaptionHWEstimatorNormal estimator = new SwaptionHWEstimatorNormal();
            EstimateQuery query = new EstimateQuery();
            query.Market = "EUR";

            EstimateRequirement[] requirements = estimator.GetRequirements(null, query);

            Assert.That(requirements.Length, Is.EqualTo(2));
            Assert.That(requirements[0].MarketDataType, Is.EqualTo(typeof(InterestRateMarketData)));
            Assert.That(requirements[1].Field, Is.EqualTo("Vol-ATM-N"));
            Assert.That(requirements[1].MarketDataType, Is.EqualTo(typeof(MatrixMarketData)));
            Assert.That(requirements[1].TickerReplacement, Is.EqualTo("SwaptionVolEUR"));
        }

        [Test]
        public void BlackModelFactory_ReturnsBachelierNormalModel()
        {
            TestableSwaptionHWEstimatorNormal estimator = new TestableSwaptionHWEstimatorNormal();
            PFunction zr = new PFunction(null);

            BlackModel result = estimator.CallBlackModelFactory(zr);

            Assert.That(result, Is.InstanceOf<BachelierNormalModel>());
        }

        [Test]
        public void ProvidesTo_IsInheritedHW1Type()
        {
            SwaptionHWEstimatorNormal estimator = new SwaptionHWEstimatorNormal();

            Assert.That(estimator.ProvidesTo, Is.EqualTo(typeof(HW1)));
        }

        [Test]
        public void Legacy_ProvidesTo_ReturnsStocasticProcessHW()
        {
            SwaptionHWEstimatorNormalLegacy estimator = new SwaptionHWEstimatorNormalLegacy();

            Assert.That(estimator.ProvidesTo, Is.EqualTo(typeof(StocasticProcessHW)));
        }

        [Test]
        public void Legacy_Description_IsInheritedFromNormal()
        {
            SwaptionHWEstimatorNormalLegacy estimator = new SwaptionHWEstimatorNormalLegacy();

            Assert.That(estimator.Description, Is.EqualTo("Calibrate HW1 using Swaption [Normal]"));
        }

        [Test]
        public void Estimate_DummyCalibration_ReturnsX0AndZeroRateCurve()
        {
            InterestRateMarketData irmd = new InterestRateMarketData();
            irmd.ZRMarketDates = (Vector)(new double[] { 1.0, 2.0, 5.0, 10.0 });
            irmd.ZRMarket = (Vector)(new double[] { 0.01, 0.015, 0.02, 0.025 });

            SwaptionsFiltering settings = new SwaptionsFiltering();
            settings.DummyCalibration = true;

            SwaptionHWEstimatorNormal estimator = new SwaptionHWEstimatorNormal();
            List<object> data = new List<object> { irmd };

            EstimationResult result = estimator.Estimate(data, settings);

            Assert.That(result.Names, Is.EqualTo(new string[] { "Alpha", "Sigma" }));
            Assert.That((double[])(Vector)result.Values, Is.EqualTo(new double[] { 0.1, 0.1 }).Within(1e-12));
            Assert.That(result.ZRX, Is.EqualTo(new double[] { 1.0, 2.0, 5.0, 10.0 }).Within(1e-12));
            Assert.That(result.ZRY, Is.EqualTo(new double[] { 0.01, 0.015, 0.02, 0.025 }).Within(1e-12));
        }
    }
}
