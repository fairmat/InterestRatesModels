using System;
using System.Collections.Generic;
using DVPLDOM;
using DVPLI;
using Fairmat.MarketData;
using NUnit.Framework;

namespace HullAndWhiteOneFactor
{
    /// <summary>
    /// Tests for <see cref="HW1DummyCalibator"/>: a pure stub estimator with no real
    /// calibration - Estimate() returns hardcoded Alpha/Sigma and passes the input curve
    /// through verbatim, so these are true numeric input-in/output-out assertions.
    /// </summary>
    [TestFixture]
    public class TestHW1DummyCalibator
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void ToolTipText_ReturnsDescriptiveText()
        {
            HW1DummyCalibator estimator = new HW1DummyCalibator();

            Assert.That(estimator.ToolTipText, Is.EqualTo("Provides default values of HW1 process"));
        }

        [Test]
        public void Description_ReturnsDescriptiveText()
        {
            HW1DummyCalibator estimator = new HW1DummyCalibator();

            Assert.That(estimator.Description, Is.EqualTo("Provides default values of HW1 process"));
        }

        [Test]
        public void ProvidesTo_ReturnsHW1Type()
        {
            HW1DummyCalibator estimator = new HW1DummyCalibator();

            Assert.That(estimator.ProvidesTo, Is.EqualTo(typeof(HW1)));
        }

        [Test]
        public void DefaultSettings_IsNull()
        {
            HW1DummyCalibator estimator = new HW1DummyCalibator();

            Assert.That(estimator.DefaultSettings, Is.Null);
        }

        [Test]
        public void GetRequirements_ReturnsSingleDiscountingCurveMarketDataRequirement()
        {
            HW1DummyCalibator estimator = new HW1DummyCalibator();

            EstimateRequirement[] requirements = estimator.GetRequirements(null, null);

            Assert.That(requirements.Length, Is.EqualTo(1));
            Assert.That(requirements[0].MarketDataType, Is.EqualTo(typeof(DiscountingCurveMarketData)));
        }

        [Test]
        public void Estimate_ReturnsHardcodedAlphaSigmaAndPassesThroughCurve()
        {
            DiscountingCurveMarketData curve = new DiscountingCurveMarketData();
            curve.Durations = (Vector)(new double[] { 1.0, 2.0, 5.0, 10.0 });
            curve.Values = (Vector)(new double[] { 0.99, 0.97, 0.90, 0.80 });

            HW1DummyCalibator estimator = new HW1DummyCalibator();
            EstimationResult result = estimator.Estimate(new List<object> { curve });

            Assert.That(result.Names, Is.EqualTo(new string[] { "Alpha", "Sigma" }));
            Assert.That((double[])(Vector)result.Values, Is.EqualTo(new double[] { 0.1, 0.05 }).Within(1e-12));
            Assert.That(result.ZRX, Is.EqualTo(new double[] { 1.0, 2.0, 5.0, 10.0 }).Within(1e-12));
            Assert.That(result.ZRY, Is.EqualTo(new double[] { 0.99, 0.97, 0.90, 0.80 }).Within(1e-12));
        }

        [Test]
        public void Estimate_EmptyData_ThrowsArgumentOutOfRange()
        {
            HW1DummyCalibator estimator = new HW1DummyCalibator();

            Assert.Throws<ArgumentOutOfRangeException>(() => estimator.Estimate(new List<object>()));
        }

        [Test]
        public void Estimate_WrongDataType_ThrowsNullReferenceException()
        {
            HW1DummyCalibator estimator = new HW1DummyCalibator();

            Assert.Throws<NullReferenceException>(() => estimator.Estimate(new List<object> { "not a curve" }));
        }
    }
}
