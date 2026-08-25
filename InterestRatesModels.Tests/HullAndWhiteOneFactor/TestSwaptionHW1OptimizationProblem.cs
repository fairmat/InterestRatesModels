using System;
using DVPLDOM;
using DVPLI;
using DVPLUtils;
using Fairmat.Optimization;
using NUnit.Framework;

namespace HullAndWhiteOneFactor
{
    [TestFixture]
    public class TestSwaptionHW1OptimizationProblem
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        private static SwaptionHW1 CreateSwaptionHW1()
        {
            Function zeroratecurve = new PFunction(null);
            zeroratecurve.Expr = new double[,] { { 0, 0.02 }, { 50, 0.02 } };
            (zeroratecurve as PFunction).m_Function.iType = EInterpolationType.LINEAR;
            return new SwaptionHW1(zeroratecurve);
        }

        private static SwaptionHW1OptimizationProblem CreateProblem(out Matrix blackSwaption, out Vector swaptionMaturity, out Vector swapDuration, out SwaptionHW1 shw1)
        {
            shw1 = CreateSwaptionHW1();
            swaptionMaturity = new Vector(new double[] { 1, 2 });
            swapDuration = new Vector(new double[] { 2, 5 });
            blackSwaption = new Matrix(new double[,] { { 10.0, 0.0 }, { 12.0, 15.0 } });
            double deltaK = 0.5;
            return new SwaptionHW1OptimizationProblem(shw1, blackSwaption, swaptionMaturity, swapDuration, deltaK);
        }

        [Test]
        public void BoundsMatchAlphaLowerBoundAndFixedUpperBounds()
        {
            Matrix blackSwaption;
            Vector swaptionMaturity, swapDuration;
            SwaptionHW1 shw1;
            SwaptionHW1OptimizationProblem problem = CreateProblem(out blackSwaption, out swaptionMaturity, out swapDuration, out shw1);

            Bounds b = problem.Bounds;

            CollectionAssert.AreEqual(new double[] { HW1.alphaLowerBound, 1e-8 }, (double[])b.Lb.ToArray());
            CollectionAssert.AreEqual(new double[] { 1 - 1e-8, 0.5 }, (double[])b.Ub.ToArray());
        }

        [Test]
        public void HasNoNonLinearOrLinearConstraints()
        {
            Matrix blackSwaption;
            Vector swaptionMaturity, swapDuration;
            SwaptionHW1 shw1;
            SwaptionHW1OptimizationProblem problem = CreateProblem(out blackSwaption, out swaptionMaturity, out swapDuration, out shw1);

            Assert.IsFalse(problem.HasNonLinearConstraints);
            Assert.IsNull(problem.LinearIneqConstraints);
        }

        [Test]
        public void ObjMatchesManualRmseSkippingZeroCells()
        {
            Matrix blackSwaption;
            Vector swaptionMaturity, swapDuration;
            SwaptionHW1 shw1;
            SwaptionHW1OptimizationProblem problem = CreateProblem(out blackSwaption, out swaptionMaturity, out swapDuration, out shw1);

            Vector x = new Vector(new double[] { 0.1, 0.01 });
            double deltaK = 0.5;

            double actual = problem.Obj(x);

            Matrix hwSwaptions = shw1.HWSwaptionMatrix(swaptionMaturity, swapDuration, x[0], x[1], deltaK);
            double sumq = 0;
            int effective = 0;
            for (int r = 0; r < swaptionMaturity.Length; r++)
                for (int c = 0; c < swapDuration.Length; c++)
                    if (blackSwaption[r, c] != 0.0)
                    {
                        sumq += Math.Pow(hwSwaptions[r, c] - blackSwaption[r, c], 2);
                        effective++;
                    }
            double expected = Math.Sqrt(sumq / effective);

            Assert.AreEqual(expected, actual, 1e-9);
        }

        [Test]
        public void ObjIsNaNWhenAllBlackSwaptionsAreZero()
        {
            SwaptionHW1 shw1 = CreateSwaptionHW1();
            Vector swaptionMaturity = new Vector(new double[] { 1, 2 });
            Vector swapDuration = new Vector(new double[] { 2 });
            Matrix blackSwaption = new Matrix(2, 1);

            SwaptionHW1OptimizationProblem problem = new SwaptionHW1OptimizationProblem(shw1, blackSwaption, swaptionMaturity, swapDuration, 0.5);

            // With no non-zero cells, "effective" stays 0 and sumq / effective is 0 / 0.
            Assert.IsTrue(double.IsNaN(problem.Obj(new Vector(new double[] { 0.1, 0.01 }))));
        }

        [Test]
        public void InternalObjWithDetailsDoesNotThrow()
        {
            Matrix blackSwaption;
            Vector swaptionMaturity, swapDuration;
            SwaptionHW1 shw1;
            SwaptionHW1OptimizationProblem problem = CreateProblem(out blackSwaption, out swaptionMaturity, out swapDuration, out shw1);

            Assert.DoesNotThrow(() => problem.Obj(new Vector(new double[] { 0.1, 0.01 }), true));
        }

        [Test]
        public void GradIsNotImplemented()
        {
            Matrix blackSwaption;
            Vector swaptionMaturity, swapDuration;
            SwaptionHW1 shw1;
            SwaptionHW1OptimizationProblem problem = CreateProblem(out blackSwaption, out swaptionMaturity, out swapDuration, out shw1);

            Assert.Throws<NotImplementedException>(() => problem.Grad(new Vector(new double[] { 0.1, 0.01 })));
        }

        [Test]
        public void GIsNotImplemented()
        {
            Matrix blackSwaption;
            Vector swaptionMaturity, swapDuration;
            SwaptionHW1 shw1;
            SwaptionHW1OptimizationProblem problem = CreateProblem(out blackSwaption, out swaptionMaturity, out swapDuration, out shw1);

            Assert.Throws<NotImplementedException>(() => problem.G(new Vector(new double[] { 0.1, 0.01 })));
        }
    }
}
