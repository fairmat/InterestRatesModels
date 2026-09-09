using System;
using System.Runtime.Serialization;
using DVPLDOM;
using DVPLI;
using NUnit.Framework;

namespace HullAndWhiteTwoFactors
{
    [TestFixture]
    public class TestHW2Process
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void DefaultConstructorSetsExpectedDefaults()
        {
            HW2 hw2 = new HW2();

            Assert.AreEqual(string.Empty, hw2._zr.Expression);
            Assert.AreEqual(0.1, hw2._a1.fV(), 1e-12);
            Assert.AreEqual(0.001, hw2._s1.fV(), 1e-12);
            Assert.AreEqual(0.01, hw2._a2.fV(), 1e-12);
            Assert.AreEqual(0.001, hw2._s2.fV(), 1e-12);
            Assert.AreEqual(0.001, hw2._rho.fV(), 1e-12);

            var exported = hw2.ExportObjects(false);
            Assert.AreEqual(0.0, ((IModelParameter)exported[6]).fV(), 1e-12);
        }

        [Test]
        public void OnDeserializedCreatesDriftAdjustmentWhenMissing()
        {
            HW2 hw2 = new HW2();
            hw2.driftAdjustment = null;

            hw2.OnDeserialized(default(StreamingContext));

            var exported = hw2.ExportObjects(false);
            Assert.IsNotNull(exported[6]);
            Assert.AreEqual(0.0, ((IModelParameter)exported[6]).fV(), 1e-12);
        }

        [Test]
        public void OnDeserializedDoesNotOverwriteExistingDriftAdjustment()
        {
            HW2 hw2 = new HW2();
            IModelParameter original = new ModelParameter(0.5);
            hw2.driftAdjustment = original;

            hw2.OnDeserialized(default(StreamingContext));

            var exported = hw2.ExportObjects(false);
            Assert.That(exported[6], Is.SameAs(original));
        }

        [Test]
        public void PopulateSetsAllFiveParameters()
        {
            HW2 hw2 = new HW2();

            hw2.Populate(
                new string[] { "alpha1", "sigma1", "alpha2", "sigma2", "rho" },
                new double[] { 0.2, 0.02, 0.03, 0.015, 0.15 });

            Assert.AreEqual(0.2, hw2._a1.fV(), 1e-12);
            Assert.AreEqual(0.02, hw2._s1.fV(), 1e-12);
            Assert.AreEqual(0.03, hw2._a2.fV(), 1e-12);
            Assert.AreEqual(0.015, hw2._s2.fV(), 1e-12);
            Assert.AreEqual(0.15, hw2._rho.fV(), 1e-12);
        }

        [Test]
        public void SetAndGetZeroRateReferenceRoundTrip()
        {
            HW2 hw2 = new HW2();

            hw2.SetZeroRateReference("@zr2");

            Assert.AreEqual("@zr2", hw2.GetZeroRateReference());
        }

        [Test]
        public void ProcessMetadataIsFixed()
        {
            HW2 hw2 = new HW2();

            Assert.AreEqual("H&W2", hw2.ProcessInfo.ProcessType);
            Assert.IsFalse(hw2.ImplementsFullSimulation);
            Assert.IsTrue(hw2.ImplementsMarkovBasedSimulation);

            DynamicInfo info = hw2.DynamicInfo;
            Assert.IsFalse(info.a_time_dependent);
            Assert.IsTrue(info.a_state_dependent);
            Assert.IsFalse(info.b_time_dependent);
            Assert.IsFalse(info.b_state_dependent);

            SimulationInfo sim = hw2.SimulationInfo;
            Assert.AreEqual(1, sim.LatentSize);
            Assert.AreEqual(2, sim.NoiseSize);
            Assert.AreEqual(2, sim.StateSize);
            CollectionAssert.AreEqual(new string[] { "short rate", "latent component" }, sim.StateDescription);

            bool[] isLog = new bool[2];
            hw2.isLog(ref isLog);
            Assert.IsFalse(isLog[0]);
            Assert.IsFalse(isLog[1]);
        }

        [Test]
        public void GetDeltaFactorsIsNullAndGetVegaFactorsReturnsBothSigmas()
        {
            HW2 hw2 = new HW2();

            Assert.IsNull(hw2.GetDeltaFactors());
            CollectionAssert.AreEqual(new IModelParameter[] { hw2._s1, hw2._s2 }, hw2.GetVegaFactors());
        }

        [Test]
        public void ExportObjectsReturnsAllSevenParametersInOrder()
        {
            HW2 hw2 = new HW2();

            var exported = hw2.ExportObjects(false);

            Assert.AreEqual(7, exported.Count);
            Assert.That(exported[0], Is.SameAs(hw2._zr));
            Assert.That(exported[1], Is.SameAs(hw2._a1));
            Assert.That(exported[2], Is.SameAs(hw2._s1));
            Assert.That(exported[3], Is.SameAs(hw2._a2));
            Assert.That(exported[4], Is.SameAs(hw2._s2));
            Assert.That(exported[5], Is.SameAs(hw2._rho));
        }

        private static ProjectROV CreateProjectWithZeroRateSymbol()
        {
            Document doc = new Document();
            ProjectROV rov = new ProjectROV(doc);
            doc.Part.Add(rov);

            AFunction zerorate = new AFunction(rov);
            zerorate.VarName = "zr";
            zerorate.m_IndependentVariables = 1;
            zerorate.m_Value = (RightValue)0.02;
            rov.Symbols.Add(zerorate);

            return rov;
        }

        [Test]
        public void ParseWithDistinctAlphasAndValidZeroRateHasNoErrors()
        {
            ProjectROV rov = CreateProjectWithZeroRateSymbol();
            HW2 hw2 = new HW2();
            hw2.SetZeroRateReference("@zr");
            hw2._a2 = (ModelParameter)0.05;

            StochasticProcessExtendible iex = new StochasticProcessExtendible(rov, hw2);
            rov.Processes.AddProcess(iex);

            bool errors = rov.Parse();

            Assert.IsFalse(errors);
            Assert.IsFalse(rov.HasErrors);
        }

        [Test]
        public void ParseFailsWhenAlpha1EqualsAlpha2()
        {
            ProjectROV rov = CreateProjectWithZeroRateSymbol();
            HW2 hw2 = new HW2();
            hw2.SetZeroRateReference("@zr");
            hw2._a2 = (ModelParameter)hw2._a1.fV();

            StochasticProcessExtendible iex = new StochasticProcessExtendible(rov, hw2);
            rov.Processes.AddProcess(iex);

            rov.Parse();

            Assert.IsTrue(rov.HasErrors);
        }

        [Test]
        public void ParseFailsWhenZeroRateReferenceIsNotAReference()
        {
            ProjectROV rov = CreateProjectWithZeroRateSymbol();
            HW2 hw2 = new HW2();
            hw2.SetZeroRateReference("not-a-reference");
            hw2._a2 = (ModelParameter)0.05;

            StochasticProcessExtendible iex = new StochasticProcessExtendible(rov, hw2);
            rov.Processes.AddProcess(iex);

            rov.Parse();

            Assert.IsTrue(rov.HasErrors);
        }

        [Test]
        public unsafe void AUsesStraightSimulationFormula()
        {
            HW2 hw2 = new HW2();
            hw2.alpha1 = 1.0;
            hw2.alpha2 = 0.05;
            hw2.theta = new double[] { 0.03, 0.04 };
            hw2.driftAdjustment = new ModelParameter(0.0);

            double[] x = new double[] { 0.02, 0.01 };
            double[] a = new double[2];

            fixed (double* px = x, pa = a)
            {
                hw2.a(0, px, pa);
            }

            Assert.AreEqual(0.03 + 0.01 - 1.0 * 0.02 + 0.0, a[0], 1e-12);
            Assert.AreEqual(-0.05 * 0.01, a[1], 1e-12);
        }

        [Test]
        public unsafe void BReturnsSigma1AndSigma2()
        {
            HW2 hw2 = new HW2();
            hw2.sigma1 = 0.02;
            hw2.sigma2 = 0.03;

            double[] x = new double[] { 0.0, 0.0 };
            double[] b = new double[2];

            fixed (double* px = x, pb = b)
            {
                hw2.b(0, px, pb);
            }

            Assert.AreEqual(0.02, b[0], 1e-12);
            Assert.AreEqual(0.03, b[1], 1e-12);
        }

        [Test]
        public unsafe void AbCombinesAAndB()
        {
            HW2 hw2 = new HW2();
            hw2.alpha1 = 1.0;
            hw2.alpha2 = 0.05;
            hw2.sigma1 = 0.02;
            hw2.sigma2 = 0.03;
            hw2.theta = new double[] { 0.03, 0.04 };
            hw2.driftAdjustment = new ModelParameter(0.0);

            double[] x = new double[] { 0.02, 0.01 };
            double[] a = new double[2];
            double[] b = new double[2];

            fixed (double* px = x, pa = a, pb = b)
            {
                hw2.ab(0, px, pa, pb);
            }

            Assert.AreEqual(0.03 + 0.01 - 1.0 * 0.02 + 0.0, a[0], 1e-12);
            Assert.AreEqual(-0.05 * 0.01, a[1], 1e-12);
            Assert.AreEqual(0.02, b[0], 1e-12);
            Assert.AreEqual(0.03, b[1], 1e-12);
        }

        [Test]
        public void X0EvaluatesZeroCurveAtSecondDateAndLeavesSecondComponentZero()
        {
            HW2 hw2 = new HW2();
            Function zr = new PFunction(null);
            zr.Expr = new double[,] { { 0, 0.01 }, { 10, 0.03 } };
            hw2.zeroRateCurve = zr;
            hw2.mDates = new double[] { 0.0, 0.5, 1.0 };

            double[] x0 = hw2.x0;

            Assert.AreEqual(2, x0.Length);
            Assert.AreEqual(zr.Evaluate(0.5), x0[0], 1e-12);
            Assert.AreEqual(0.0, x0[1], 1e-12);
        }
    }
}
