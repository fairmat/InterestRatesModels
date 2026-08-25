using System;
using DVPLDOM;
using DVPLI;
using NUnit.Framework;

namespace CIRProcess
{
    [TestFixture]
    public class TestCIR2
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void DefaultConstructorSetsExpectedDefaults()
        {
            CIR2 cir2 = new CIR2();

            Assert.AreEqual(0.026758131, cir2.k1.fV(), 1e-12);
            Assert.AreEqual(0.226406137, cir2.k2.fV(), 1e-12);
            Assert.AreEqual(0.023147821, cir2.theta1.fV(), 1e-12);
            Assert.AreEqual(0.01, cir2.theta2.fV(), 1e-12);
            Assert.AreEqual(0.1, cir2.sigma1.fV(), 1e-12);
            Assert.AreEqual(0.0001, cir2.sigma2.fV(), 1e-12);
            Assert.AreEqual(0.00001, cir2.startingValue1.fV(), 1e-12);
            Assert.AreEqual(0.00001, cir2.startingValue2.fV(), 1e-12);

            Assert.AreEqual("k1", cir2.k1.Description);
            Assert.AreEqual("k2", cir2.k2.Description);
            Assert.AreEqual("Theta1", cir2.theta1.Description);
            Assert.AreEqual("Theta2", cir2.theta2.Description);
            Assert.AreEqual("Sigma1", cir2.sigma1.Description);
            Assert.AreEqual("Sigma2", cir2.sigma2.Description);
            Assert.AreEqual("starting value factor 1", cir2.startingValue1.Description);
            Assert.AreEqual("starting value factor 2", cir2.startingValue2.Description);
            Assert.AreEqual("zero rate function reference", cir2.ZRReference.Description);
        }

        [Test]
        public void ProcessTypeIsFixed()
        {
            Assert.AreEqual("Interest Rate Models/CIR (Two Factors)", CIR2.ProcessType);
        }

        [Test]
        public void ExportObjectsReturnsAllNineParametersInOrder()
        {
            CIR2 cir2 = new CIR2();

            var exported = cir2.ExportObjects(false);

            Assert.AreEqual(9, exported.Count);
            Assert.That(exported[0], Is.SameAs(cir2.startingValue1));
            Assert.That(exported[1], Is.SameAs(cir2.startingValue2));
            Assert.That(exported[2], Is.SameAs(cir2.k1));
            Assert.That(exported[3], Is.SameAs(cir2.k2));
            Assert.That(exported[4], Is.SameAs(cir2.theta1));
            Assert.That(exported[5], Is.SameAs(cir2.theta2));
            Assert.That(exported[6], Is.SameAs(cir2.sigma1));
            Assert.That(exported[7], Is.SameAs(cir2.sigma2));
            Assert.That(exported[8], Is.SameAs(cir2.ZRReference));
        }

        [Test]
        public void X0ReturnsBothStartingValues()
        {
            CIR2 cir2 = new CIR2();
            CollectionAssert.AreEqual(new double[] { cir2.startingValue1.fV(), cir2.startingValue2.fV() }, cir2.x0);
        }

        [Test]
        public void IsLogSetsBothComponentsFalse()
        {
            CIR2 cir2 = new CIR2();
            bool[] isLog = new bool[2];

            cir2.isLog(ref isLog);

            Assert.IsFalse(isLog[0]);
            Assert.IsFalse(isLog[1]);
        }

        [Test]
        public void ProcessMetadataIsFixed()
        {
            CIR2 cir2 = new CIR2();

            Assert.IsFalse(cir2.ImplementsFullSimulation);
            Assert.IsTrue(cir2.ImplementsMarkovBasedSimulation);
            Assert.AreEqual(CIR2.ProcessType, cir2.ProcessInfo.ProcessType);

            DynamicInfo info = cir2.DynamicInfo;
            Assert.IsFalse(info.a_time_dependent);
            Assert.IsTrue(info.a_state_dependent);
            Assert.IsFalse(info.b_time_dependent);
            Assert.IsTrue(info.b_state_dependent);

            SimulationInfo sim = cir2.SimulationInfo;
            Assert.AreEqual(0, sim.LatentSize);
            Assert.AreEqual(2, sim.NoiseSize);
            Assert.AreEqual(2, sim.StateSize);
            CollectionAssert.AreEqual(new string[] { "s1", "s2" }, sim.StateDescription);
        }

        [Test]
        public void SetupIsANoOp()
        {
            CIR2 cir2 = new CIR2();
            Assert.DoesNotThrow(() => cir2.Setup(new double[] { 0.0, 1.0 }));
        }

        [Test]
        public void PopulateThrowsNotImplemented()
        {
            CIR2 cir2 = new CIR2();
            Assert.Throws<NotImplementedException>(() =>
                cir2.Populate(new string[] { "k1" }, new double[] { 0.1 }));
        }

        [Test]
        public unsafe void AUsesRawStateWhileBFloorsNegativeState()
        {
            CIR2 cir2 = new CIR2();

            double[] x = new double[] { -0.5, 0.2 };
            double[] a = new double[2];
            double[] b = new double[2];

            fixed (double* px = x, pa = a, pb = b)
            {
                cir2.ab(0, px, pa, pb);
            }

            Assert.AreEqual(cir2.k1.fV() * (cir2.theta1.fV() - (-0.5)), a[0], 1e-12);
            Assert.AreEqual(cir2.k2.fV() * (cir2.theta2.fV() - 0.2), a[1], 1e-12);

            // b floors negative state to zero via AdjSqrt.
            Assert.AreEqual(0.0, b[0], 1e-12);
            Assert.AreEqual(Math.Sqrt(0.2) * cir2.sigma2.fV(), b[1], 1e-12);
        }

        [Test]
        public void TransformFloorsBothComponentsToTinyPositiveValue()
        {
            CIR2 cir2 = new CIR2();
            Matrix outDynamic = new Matrix(new double[,]
            {
                { -0.5, 0.0000000002 },
                { 0.01, -1.0 },
            });

            cir2.Transform(new double[] { 0.0, 1.0 }, outDynamic);

            Assert.AreEqual(1e-10, outDynamic[0, 0], 1e-20);
            Assert.AreEqual(0.0000000002, outDynamic[0, 1], 1e-20);
            Assert.AreEqual(0.01, outDynamic[1, 0], 1e-12);
            Assert.AreEqual(1e-10, outDynamic[1, 1], 1e-20);
        }

        [Test]
        public void CDSDefaultPThrowsWhenCalledWithoutParse()
        {
            CIR2 cir2 = new CIR2();
            Matrix dynamic = new Matrix(new double[,] { { 0.001, 0.001 } });

            Assert.Throws<NullReferenceException>(() =>
                cir2.CDSDefaultP(dynamic, new double[] { 0.0 }, 0, 1.0));
        }

        private static ProjectROV CreateProjectWithFlatZeroCurve()
        {
            Document doc = new Document();
            ProjectROV rov = new ProjectROV(doc);
            doc.Part.Add(rov);

            Function zr = new PFunction(null);
            zr.VarName = "ZR";
            zr.Expr = new double[,] { { 0.0, 0.02 }, { 10.0, 0.02 } };
            rov.Symbols.Add(zr);
            rov.Parse();

            return rov;
        }

        [Test]
        public void ParseWithZeroRateSymbolHasNoErrors()
        {
            CIR2 cir2 = new CIR2();
            ProjectROV rov = CreateProjectWithFlatZeroCurve();

            bool errors = cir2.Parse(rov);

            Assert.IsFalse(errors);
        }

        [Test]
        public void CDSDefaultPReturnsFiniteValueForWellFormedInputs()
        {
            CIR2 cir2 = new CIR2();
            ProjectROV rov = CreateProjectWithFlatZeroCurve();
            cir2.Parse(rov);

            Matrix dynamic = new Matrix(new double[,] { { 0.00001, 0.00001 } });
            double result = cir2.CDSDefaultP(dynamic, new double[] { 0.0 }, 0, 1.0);

            Assert.IsTrue(double.IsFinite(result));
        }

        [Test]
        public void CDSDefaultPIsNonFiniteWhenTauIsZero()
        {
            CIR2 cir2 = new CIR2();
            ProjectROV rov = CreateProjectWithFlatZeroCurve();
            cir2.Parse(rov);

            Matrix dynamic = new Matrix(new double[,] { { 0.00001, 0.00001 } });
            double result = cir2.CDSDefaultP(dynamic, new double[] { 1.0 }, 0, 1.0);

            Assert.IsFalse(double.IsFinite(result));
        }
    }
}
