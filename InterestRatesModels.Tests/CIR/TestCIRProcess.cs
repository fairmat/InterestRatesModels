using System;
using System.Reflection;
using DVPLDOM;
using DVPLI;
using NUnit.Framework;

namespace CIRProcess
{
    [TestFixture]
    public class TestCIRProcess
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        private static double GetPrivateDouble(CIR cir, string fieldName)
        {
            FieldInfo field = typeof(CIR).GetField(fieldName, BindingFlags.NonPublic | BindingFlags.Instance);
            return (double)field.GetValue(cir);
        }

        [Test]
        public void DefaultConstructorSetsAllParametersToPoint001()
        {
            CIR cir = new CIR();

            Assert.AreEqual(0.001, cir.k.fV(), 1e-12);
            Assert.AreEqual(0.001, cir.theta.fV(), 1e-12);
            Assert.AreEqual(0.001, cir.sigma.fV(), 1e-12);
            Assert.AreEqual(0.001, cir.r0.fV(), 1e-12);
        }

        [Test]
        public void ExportObjectsReturnsAllFourParametersInOrder()
        {
            CIR cir = new CIR();

            var exported = cir.ExportObjects(false);

            Assert.AreEqual(4, exported.Count);
            Assert.That(exported[0], Is.SameAs(cir.k));
            Assert.That(exported[1], Is.SameAs(cir.theta));
            Assert.That(exported[2], Is.SameAs(cir.sigma));
            Assert.That(exported[3], Is.SameAs(cir.r0));
        }

        [Test]
        public void X0ReturnsR0()
        {
            CIR cir = new CIR();
            Assert.AreEqual(new double[] { cir.r0.fV() }, cir.x0);
        }

        [Test]
        public void GetVegaFactorsReturnsSigmaAndGetDeltaFactorsReturnsNull()
        {
            CIR cir = new CIR();

            CollectionAssert.AreEqual(new IModelParameter[] { cir.sigma }, cir.GetVegaFactors());
            Assert.IsNull(cir.GetDeltaFactors());
        }

        [Test]
        public void IsLogSetsFalse()
        {
            CIR cir = new CIR();
            bool[] isLog = new bool[1];

            cir.isLog(ref isLog);

            Assert.IsFalse(isLog[0]);
        }

        [Test]
        public void ProcessMetadataIsFixed()
        {
            CIR cir = new CIR();

            Assert.IsFalse(cir.ImplementsFullSimulation);
            Assert.IsTrue(cir.ImplementsMarkovBasedSimulation);
            Assert.AreEqual("CIR", cir.ProcessInfo.ProcessType);

            DynamicInfo info = cir.DynamicInfo;
            Assert.IsFalse(info.a_time_dependent);
            Assert.IsTrue(info.a_state_dependent);
            Assert.IsFalse(info.b_time_dependent);
            Assert.IsFalse(info.b_state_dependent);

            SimulationInfo sim = cir.SimulationInfo;
            Assert.AreEqual(0, sim.LatentSize);
            Assert.AreEqual(1, sim.NoiseSize);
            Assert.AreEqual(1, sim.StateSize);
            Assert.AreEqual(-1, sim.DefaultComponent);
            CollectionAssert.AreEqual(new string[] { "short rate" }, sim.StateDescription);
        }

        private static ProjectROV CreateProject()
        {
            Document doc = new Document();
            ProjectROV rov = new ProjectROV(doc);
            doc.Part.Add(rov);
            return rov;
        }

        [Test]
        public void ParseWithPlainNumericParametersHasNoErrors()
        {
            CIR cir = new CIR();
            ProjectROV rov = CreateProject();

            bool errors = cir.Parse(rov);

            Assert.IsFalse(errors);
        }

        [Test]
        public void SetupComputesDAndNu()
        {
            CIR cir = new CIR();
            cir.k = new ModelParameter(1.0);
            cir.theta = new ModelParameter(0.02);
            cir.sigma = new ModelParameter(0.08);
            cir.r0 = new ModelParameter(0.01);
            ProjectROV rov = CreateProject();
            cir.Parse(rov);

            cir.Setup(new double[] { 0.0, 1.0 });

            double expectedD = Math.Sqrt(1.0 * 1.0 + 2.0 * 0.08 * 0.08);
            double expectedNu = 2.0 * 1.0 * 0.02 / (0.08 * 0.08);

            Assert.AreEqual(expectedD, GetPrivateDouble(cir, "d"), 1e-12);
            Assert.AreEqual(expectedNu, GetPrivateDouble(cir, "nu"), 1e-12);
        }

        [Test]
        public void SetupWithZeroSigmaMakesNuDivideByZero()
        {
            CIR cir = new CIR();
            cir.k = new ModelParameter(1.0);
            cir.theta = new ModelParameter(0.02);
            cir.sigma = new ModelParameter(0.0);
            cir.r0 = new ModelParameter(0.01);
            ProjectROV rov = CreateProject();
            cir.Parse(rov);

            cir.Setup(new double[] { 0.0, 1.0 });

            // nu = 2 * k * theta / sigma^2, with sigma == 0 this is a division by
            // zero and yields a non-finite (infinite) value.
            Assert.IsFalse(double.IsFinite(GetPrivateDouble(cir, "nu")));
        }

        [Test]
        public void BondMatchesClosedFormFormula()
        {
            double k = 1.0;
            double theta = 0.02;
            double sigma = 0.08;
            double r0 = 0.01;

            CIR cir = new CIR();
            cir.k = new ModelParameter(k);
            cir.theta = new ModelParameter(theta);
            cir.sigma = new ModelParameter(sigma);
            cir.r0 = new ModelParameter(r0);
            ProjectROV rov = CreateProject();
            cir.Parse(rov);
            cir.Setup(new double[] { 0.0, 1.0, 2.0 });

            Matrix dynamic = new Matrix(new double[,] { { r0 } });
            double t = 1.0;
            double s = 2.0;
            double bond = cir.Bond(dynamic, new double[] { 0.0 }, 0, t, s);

            double d = Math.Sqrt(k * k + 2.0 * sigma * sigma);
            double nu = 2.0 * k * theta / (sigma * sigma);
            double T = s - t;
            double den = (k + d) * (Math.Exp(d * T) - 1.0) + 2.0 * d;
            double A = Math.Pow(2.0 * d * Math.Exp(0.5 * (k + d) * T) / den, nu);
            double B = 2.0 * (Math.Exp(d * T) - 1.0) / den;
            double expected = A * Math.Exp(-r0 * B);

            Assert.AreEqual(expected, bond, 1e-12);
        }

        [Test]
        public unsafe void AAndBFlipToZeroForNegativeState()
        {
            CIR cir = new CIR();
            cir.k = new ModelParameter(1.0);
            cir.theta = new ModelParameter(0.02);
            cir.sigma = new ModelParameter(0.08);
            cir.r0 = new ModelParameter(0.01);
            ProjectROV rov = CreateProject();
            cir.Parse(rov);

            double x = -0.5;
            double aOut = 0, bOut = 0;
            cir.a(0, &x, &aOut);
            cir.b(0, &x, &bOut);

            Assert.AreEqual(1.0 * (0.02 - 0.0), aOut, 1e-12);
            Assert.AreEqual(0.08 * Math.Sqrt(0.0), bOut, 1e-12);
        }

        [Test]
        public unsafe void AAndBUsePositiveStateDirectly()
        {
            CIR cir = new CIR();
            cir.k = new ModelParameter(1.0);
            cir.theta = new ModelParameter(0.02);
            cir.sigma = new ModelParameter(0.08);
            cir.r0 = new ModelParameter(0.01);
            ProjectROV rov = CreateProject();
            cir.Parse(rov);

            double x = 0.03;
            double aOut = 0, bOut = 0;
            cir.ab(0, &x, &aOut, &bOut);

            Assert.AreEqual(1.0 * (0.02 - 0.03), aOut, 1e-12);
            Assert.AreEqual(0.08 * Math.Sqrt(0.03), bOut, 1e-12);
        }

        [Test]
        public void PopulateUpdatesAllParametersWhenAllNamesPresent()
        {
            CIR cir = new CIR();
            EstimationResult estimate = new EstimationResult(CIR.parameterNames,
                new double[] { 0.5, 0.03, 0.09, 0.012 });

            cir.Populate(null, estimate);

            Assert.AreEqual(0.5, cir.k.fV(), 1e-12);
            Assert.AreEqual(0.03, cir.theta.fV(), 1e-12);
            Assert.AreEqual(0.09, cir.sigma.fV(), 1e-12);
            Assert.AreEqual(0.012, cir.r0.fV(), 1e-12);
        }

        [Test]
        public void PopulateDoesNotThrowWhenAParameterNameIsMissing()
        {
            CIR cir = new CIR();
            EstimationResult estimate = new EstimationResult(
                new string[] { "k", "Theta", "Sigma" },
                new double[] { 0.5, 0.03, 0.09 });

            Assert.DoesNotThrow(() => cir.Populate(null, estimate));
        }
    }
}
