using System;
using DVPLDOM;
using DVPLI;
using DVPLUtils;
using Fairmat.Optimization;
using NUnit.Framework;

namespace HullAndWhiteOneFactor
{
    [TestFixture]
    public class TestCapsHW1OptimizationProblem
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        private static CapHW1 CreateCapHW1()
        {
            double[,] values =
            {
                { .5, .015 },
                { 1, .0165 },
                { 2, .0168 },
                { 3, .0172 },
                { 5, .0182 },
                { 8, .0210 },
                { 10, .025 },
            };

            Function zeroratecurve = new PFunction(null);
            zeroratecurve.Expr = values;
            (zeroratecurve as PFunction).m_Function.iType = EInterpolationType.LINEAR;

            return new CapHW1(zeroratecurve);
        }

        private static CapsHW1OptimizationProblem CreateProblem(out Matrix blackCaps, out Vector capMaturity, out Vector capRate, out CapHW1 hw1Caps)
        {
            hw1Caps = CreateCapHW1();
            capMaturity = new Vector(new double[] { 1, 2 });
            capRate = new Vector(new double[] { 0.01, 0.02 });
            blackCaps = new Matrix(new double[,] { { 0.001, 0.0 }, { 0.002, 0.003 } });
            double deltaK = 0.5;
            return new CapsHW1OptimizationProblem(hw1Caps, blackCaps, capMaturity, capRate, deltaK);
        }

        [Test]
        public void BoundsMatchAlphaLowerBoundAndFixedUpperBounds()
        {
            Matrix blackCaps;
            Vector capMaturity, capRate;
            CapHW1 hw1Caps;
            CapsHW1OptimizationProblem problem = CreateProblem(out blackCaps, out capMaturity, out capRate, out hw1Caps);

            Bounds b = problem.Bounds;

            CollectionAssert.AreEqual(new double[] { HW1.alphaLowerBound, 1e-8 }, (double[])b.Lb.ToArray());
            CollectionAssert.AreEqual(new double[] { 1 - 1e-8, 0.5 }, (double[])b.Ub.ToArray());
        }

        [Test]
        public void HasNoNonLinearOrLinearConstraints()
        {
            Matrix blackCaps;
            Vector capMaturity, capRate;
            CapHW1 hw1Caps;
            CapsHW1OptimizationProblem problem = CreateProblem(out blackCaps, out capMaturity, out capRate, out hw1Caps);

            Assert.IsFalse(problem.HasNonLinearConstraints);
            Assert.IsNull(problem.LinearIneqConstraints);
        }

        [Test]
        public void ObjMatchesManualL2NormSkippingZeroCells()
        {
            Matrix blackCaps;
            Vector capMaturity, capRate;
            CapHW1 hw1Caps;
            CapsHW1OptimizationProblem problem = CreateProblem(out blackCaps, out capMaturity, out capRate, out hw1Caps);

            Vector x = new Vector(new double[] { 0.05, 0.02 });

            double actual = problem.Obj(x);

            Matrix hwCaps = hw1Caps.HWMatrixCaps(capMaturity, capRate, x[0], x[1], 0.5);
            double sum = 0;
            for (int r = 0; r < hwCaps.R; r++)
                for (int c = 0; c < hwCaps.C; c++)
                    if (blackCaps[r, c] != 0.0)
                        sum += Math.Pow(hwCaps[r, c] - blackCaps[r, c], 2);
            double expected = Math.Sqrt(sum);

            Assert.AreEqual(expected, actual, 1e-12);
        }

        [Test]
        public void ObjIsZeroWhenAllBlackCapsAreZero()
        {
            CapHW1 hw1Caps = CreateCapHW1();
            Vector capMaturity = new Vector(new double[] { 1, 2 });
            Vector capRate = new Vector(new double[] { 0.01, 0.02 });
            Matrix blackCaps = new Matrix(2, 2);

            CapsHW1OptimizationProblem problem = new CapsHW1OptimizationProblem(hw1Caps, blackCaps, capMaturity, capRate, 0.5);

            Assert.AreEqual(0.0, problem.Obj(new Vector(new double[] { 0.05, 0.02 })));
        }

        [Test]
        public void GradIsNotImplemented()
        {
            Matrix blackCaps;
            Vector capMaturity, capRate;
            CapHW1 hw1Caps;
            CapsHW1OptimizationProblem problem = CreateProblem(out blackCaps, out capMaturity, out capRate, out hw1Caps);

            Assert.Throws<NotImplementedException>(() => problem.Grad(new Vector(new double[] { 0.05, 0.02 })));
        }

        [Test]
        public void GIsNotImplemented()
        {
            Matrix blackCaps;
            Vector capMaturity, capRate;
            CapHW1 hw1Caps;
            CapsHW1OptimizationProblem problem = CreateProblem(out blackCaps, out capMaturity, out capRate, out hw1Caps);

            Assert.Throws<NotImplementedException>(() => problem.G(new Vector(new double[] { 0.05, 0.02 })));
        }
    }
}
