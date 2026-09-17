
using System;
using System.Collections.Generic;
using DVPLI;
using Fairmat.Calibration;
using Fairmat.Finance;
using NUnit.Framework;

namespace HullAndWhiteOneFactor
{
    /// <summary>
    /// Tests for <see cref="SwaptionHW1RealWorldEstimator"/> covering the parts of the class
    /// that don't require constructing real Fairmat market data fixtures: descriptive metadata,
    /// the guard clause in Estimate, defensive-casting behavior of the settings parameter,
    /// and the explicit IEstimator/IEstimatorEx members.
    /// </summary>
    [TestFixture]
    public class TestSwaptionHW1RealWorldEstimator
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void ToolTipText_ReturnsDescriptiveText()
        {
            SwaptionHW1RealWorldEstimator estimator = new SwaptionHW1RealWorldEstimator();

            Assert.That(estimator.ToolTipText, Is.EqualTo("Real world calibration using swaptions and zero rate curve historical series"));
        }

        [Test]
        public void Description_ReturnsDescriptiveText()
        {
            SwaptionHW1RealWorldEstimator estimator = new SwaptionHW1RealWorldEstimator();

            Assert.That(estimator.Description, Is.EqualTo("Real world calibration using swaptions and zero rate curve historical series"));
        }

        [Test]
        public void ProvidesTo_ReturnsHW1Type()
        {
            SwaptionHW1RealWorldEstimator estimator = new SwaptionHW1RealWorldEstimator();

            Assert.That(estimator.ProvidesTo, Is.EqualTo(typeof(HW1)));
        }

        [Test]
        public void DefaultSettings_ReturnsLambdaCalibrationSettings()
        {
            SwaptionHW1RealWorldEstimator estimator = new SwaptionHW1RealWorldEstimator();

            Assert.That(estimator.DefaultSettings, Is.InstanceOf<LambdaCalibrationSettings>());
        }

        [Test]
        public void Estimate_YearsGreaterThanBondMaturity_ThrowsException()
        {
            SwaptionHW1RealWorldEstimator estimator = new SwaptionHW1RealWorldEstimator();
            LambdaCalibrationSettings settings = new LambdaCalibrationSettings { Years = 10, BondMaturity = 5 };

            Exception ex = Assert.Throws<Exception>(() => estimator.Estimate(new List<object>(), settings));

            Assert.That(ex.Message, Is.EqualTo("Bond maturity has to be greater of the historical series time span."));
        }

        [Test]
        public void Estimate_NullSettings_ThrowsNullReferenceException()
        {
            SwaptionHW1RealWorldEstimator estimator = new SwaptionHW1RealWorldEstimator();

            Assert.Throws<NullReferenceException>(() => estimator.Estimate(new List<object>(), null));
        }

        [Test]
        public void Estimate_WrongSettingsType_ThrowsInvalidCastException()
        {
            SwaptionHW1RealWorldEstimator estimator = new SwaptionHW1RealWorldEstimator();
            SwaptionsFiltering settings = new SwaptionsFiltering();

            Assert.Throws<InvalidCastException>(() => estimator.Estimate(new List<object>(), settings));
        }

        [Test]
        public void IEstimator_Estimate_ThrowsNotImplementedException()
        {
            IEstimator estimator = new SwaptionHW1RealWorldEstimator();

            Assert.Throws<NotImplementedException>(() => estimator.Estimate(new List<object>(), null, null, null));
        }

        [Test]
        public void IEstimator_GetRequirements_ThrowsNotImplementedException()
        {
            IEstimator estimator = new SwaptionHW1RealWorldEstimator();

            Assert.Throws<NotImplementedException>(() => estimator.GetRequirements(null, default));
        }

        [Test]
        public void IEstimator_ProvidesTo_ThrowsNotImplementedException()
        {
            IEstimator estimator = new SwaptionHW1RealWorldEstimator();

            Assert.Throws<NotImplementedException>(() => { Type _ = estimator.ProvidesTo; });
        }

        [Test]
        public void IEstimatorEx_DefaultSettings_ThrowsNotImplementedException()
        {
            IEstimatorEx estimator = new SwaptionHW1RealWorldEstimator();

            Assert.Throws<NotImplementedException>(() => { IEstimationSettings _ = estimator.DefaultSettings; });
        }
    }
}
