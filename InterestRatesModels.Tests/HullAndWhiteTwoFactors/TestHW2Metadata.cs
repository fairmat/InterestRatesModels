using System;
using System.Collections.Generic;
using DVPLDOM;
using DVPLI;
using NUnit.Framework;

namespace HullAndWhiteTwoFactors
{
    /// <summary>
    /// Tests for HW2's metadata/plumbing surface: none of these need a live simulation.
    /// Numeric assertions read back ModelParameter values via ExportObjects(false) + .fV(),
    /// which is safe here because ModelParameter(double, string) stores a plain RightValue
    /// whose fV() returns the constant immediately, with no Parse(IProject) needed.
    /// </summary>
    [TestFixture]
    public class TestHW2Metadata
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void DefaultConstructor_SetsDocumentedDefaultParameterValues()
        {
            HW2 hw2 = new HW2();

            List<IExportable> objects = hw2.ExportObjects(false);

            Assert.That(objects.Count, Is.EqualTo(7));
            Assert.That(((ModelParameter)objects[1]).fV(), Is.EqualTo(0.1).Within(1e-12));
            Assert.That(((ModelParameter)objects[2]).fV(), Is.EqualTo(0.001).Within(1e-12));
            Assert.That(((ModelParameter)objects[3]).fV(), Is.EqualTo(0.01).Within(1e-12));
            Assert.That(((ModelParameter)objects[4]).fV(), Is.EqualTo(0.001).Within(1e-12));
            Assert.That(((ModelParameter)objects[5]).fV(), Is.EqualTo(0.001).Within(1e-12));
        }

        [Test]
        public void Populate_SetsAlphaAndSigmaAndRhoFromNamesAndValues()
        {
            HW2 hw2 = new HW2();

            string[] names = new string[] { "alpha1", "sigma1", "alpha2", "sigma2", "rho" };
            double[] values = new double[] { 0.2, 0.02, 0.03, 0.015, 0.4 };
            hw2.Populate(names, values);

            List<IExportable> objects = hw2.ExportObjects(false);

            Assert.That(((ModelParameter)objects[1]).fV(), Is.EqualTo(0.2).Within(1e-12));
            Assert.That(((ModelParameter)objects[2]).fV(), Is.EqualTo(0.02).Within(1e-12));
            Assert.That(((ModelParameter)objects[3]).fV(), Is.EqualTo(0.03).Within(1e-12));
            Assert.That(((ModelParameter)objects[4]).fV(), Is.EqualTo(0.015).Within(1e-12));
            Assert.That(((ModelParameter)objects[5]).fV(), Is.EqualTo(0.4).Within(1e-12));
        }

        [Test]
        public void Populate_AcceptsAliasNames()
        {
            HW2 hw2 = new HW2();

            string[] names = new string[] { "a1", "sigma", "a2", "sigma2", "rho" };
            double[] values = new double[] { 0.25, 0.03, 0.04, 0.02, 0.5 };
            hw2.Populate(names, values);

            List<IExportable> objects = hw2.ExportObjects(false);

            Assert.That(((ModelParameter)objects[1]).fV(), Is.EqualTo(0.25).Within(1e-12));
            Assert.That(((ModelParameter)objects[2]).fV(), Is.EqualTo(0.03).Within(1e-12));
        }

        [Test]
        public void ZeroRateReference_RoundTripsThroughSetAndGet()
        {
            HW2 hw2 = new HW2();

            hw2.SetZeroRateReference("@zr1");

            Assert.That(hw2.GetZeroRateReference(), Is.EqualTo("@zr1"));
        }

        [Test]
        public void ProcessInfo_ReturnsHW2ProcessType()
        {
            HW2 hw2 = new HW2();

            Assert.That(hw2.ProcessInfo.ProcessType, Is.EqualTo("H&W2"));
        }

        [Test]
        public void ImplementsFullSimulation_IsFalse()
        {
            HW2 hw2 = new HW2();

            Assert.That(hw2.ImplementsFullSimulation, Is.False);
        }

        [Test]
        public void ImplementsMarkovBasedSimulation_IsTrue()
        {
            HW2 hw2 = new HW2();

            Assert.That(hw2.ImplementsMarkovBasedSimulation, Is.True);
        }

        [Test]
        public void SimulationInfo_ReturnsExpectedShape()
        {
            HW2 hw2 = new HW2();

            SimulationInfo info = hw2.SimulationInfo;

            Assert.That(info.LatentSize, Is.EqualTo(1));
            Assert.That(info.NoiseSize, Is.EqualTo(2));
            Assert.That(info.StateSize, Is.EqualTo(2));
            Assert.That(info.StateDescription, Is.EqualTo(new string[] { "short rate", "latent component" }));
        }

        [Test]
        public void GetDeltaFactors_ReturnsNull()
        {
            HW2 hw2 = new HW2();

            Assert.That(hw2.GetDeltaFactors(), Is.Null);
        }

        [Test]
        public void GetVegaFactors_ReturnsSigma1AndSigma2()
        {
            HW2 hw2 = new HW2();

            IModelParameter[] vegaFactors = hw2.GetVegaFactors();

            Assert.That(vegaFactors.Length, Is.EqualTo(2));
            Assert.That(((ModelParameter)vegaFactors[0]).fV(), Is.EqualTo(0.001).Within(1e-12));
            Assert.That(((ModelParameter)vegaFactors[1]).fV(), Is.EqualTo(0.001).Within(1e-12));
        }
    }
}
