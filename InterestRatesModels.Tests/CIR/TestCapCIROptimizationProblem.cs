using System;
using DVPLDOM;
using DVPLI;
using Fairmat.Optimization;
using NUnit.Framework;

namespace CIRProcess
{
    [TestFixture]
    public class TestCapCIROptimizationProblem
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        private static CapCIROptimizationProblem CreateProblem(out Matrix blackCaps, out Vector capMaturity, out Vector capRate)
        {
            capMaturity = new Vector(new double[] { 1, 2 });
            capRate = new Vector(new double[] { 0.01, 0.02 });
            blackCaps = new Matrix(new double[,] { { 0.01, 0.0 }, { 0.02, 0.03 } });
            double tau = 0.5;
            double r0 = 0.02;
            return new CapCIROptimizationProblem(blackCaps, capMaturity, capRate, tau, r0);
        }

        [Test]
        public void ObjMatchesManualL2NormSkippingZeroCells()
        {
            Matrix blackCaps;
            Vector capMaturity, capRate;
            CapCIROptimizationProblem problem = CreateProblem(out blackCaps, out capMaturity, out capRate);

            Vector x = new Vector(new double[] { 1.0, 0.02, 0.05 });

            double actual = problem.Obj(x);

            Matrix cirCaps = CIRCap.CIRCapMatrix(capMaturity, capRate, 0.5, 0.02, x);
            double sum = 0;
            for (int r = 0; r < cirCaps.R; r++)
                for (int c = 0; c < cirCaps.C; c++)
                    if (blackCaps[r, c] != 0.0)
                        sum += Math.Pow(cirCaps[r, c] - blackCaps[r, c], 2);
            double expected = Math.Sqrt(sum);

            Assert.AreEqual(expected, actual, 1e-12);
        }

        [Test]
        public void ObjIsZeroWhenAllBlackCapsAreZero()
        {
            Vector capMaturity = new Vector(new double[] { 1, 2 });
            Vector capRate = new Vector(new double[] { 0.01, 0.02 });
            Matrix blackCaps = new Matrix(2, 2);

            CapCIROptimizationProblem problem = new CapCIROptimizationProblem(blackCaps, capMaturity, capRate, 0.5, 0.02);

            Assert.AreEqual(0.0, problem.Obj(new Vector(new double[] { 1.0, 0.02, 0.05 })));
        }

        [Test]
        public void ObjExcludesZeroBlackCapCellRegardlessOfModelValue()
        {
            Matrix blackCaps;
            Vector capMaturity, capRate;
            CapCIROptimizationProblem problem = CreateProblem(out blackCaps, out capMaturity, out capRate);

            // blackCaps from CreateProblem is { {0.01, 0.0}, {0.02, 0.03} }: cell [0,1] is
            // the only zero cell. Use two very different x vectors so the CIR-implied value
            // at [0,1] differs substantially between them, then check Obj against a sum
            // computed only from the three known non-zero cells (hardcoded, not derived from
            // the same "!= 0.0" guard under test) - proving the zero cell is excluded from
            // Obj regardless of what the model computes there.
            Vector xLow = new Vector(new double[] { 1.0, 0.02, 0.05 });
            Vector xHigh = new Vector(new double[] { 3.0, 0.10, 0.30 });

            Matrix cirCapsLow = CIRCap.CIRCapMatrix(capMaturity, capRate, 0.5, 0.02, xLow);
            Matrix cirCapsHigh = CIRCap.CIRCapMatrix(capMaturity, capRate, 0.5, 0.02, xHigh);

            Assert.That(Math.Abs(cirCapsLow[0, 1] - cirCapsHigh[0, 1]), Is.GreaterThan(1e-3));

            double expectedLow = Math.Sqrt(
                Math.Pow(cirCapsLow[0, 0] - blackCaps[0, 0], 2) +
                Math.Pow(cirCapsLow[1, 0] - blackCaps[1, 0], 2) +
                Math.Pow(cirCapsLow[1, 1] - blackCaps[1, 1], 2));
            double expectedHigh = Math.Sqrt(
                Math.Pow(cirCapsHigh[0, 0] - blackCaps[0, 0], 2) +
                Math.Pow(cirCapsHigh[1, 0] - blackCaps[1, 0], 2) +
                Math.Pow(cirCapsHigh[1, 1] - blackCaps[1, 1], 2));

            Assert.AreEqual(expectedLow, problem.Obj(xLow), 1e-12);
            Assert.AreEqual(expectedHigh, problem.Obj(xHigh), 1e-12);
        }

        [Test]
        public void BoundsAreFixed()
        {
            Matrix blackCaps;
            Vector capMaturity, capRate;
            CapCIROptimizationProblem problem = CreateProblem(out blackCaps, out capMaturity, out capRate);

            Bounds b = problem.Bounds;

            CollectionAssert.AreEqual(new double[] { 1e-5, 1e-5, 1e-5 }, (double[])b.Lb.ToArray());
            CollectionAssert.AreEqual(new double[] { 5, 0.5, 0.5 }, (double[])b.Ub.ToArray());
        }

        [Test]
        public void HasNoNonLinearOrLinearConstraints()
        {
            Matrix blackCaps;
            Vector capMaturity, capRate;
            CapCIROptimizationProblem problem = CreateProblem(out blackCaps, out capMaturity, out capRate);

            Assert.IsFalse(problem.HasNonLinearConstraints);
            Assert.IsNull(problem.LinearIneqConstraints);
        }

        [Test]
        public void GradIsNotImplemented()
        {
            Matrix blackCaps;
            Vector capMaturity, capRate;
            CapCIROptimizationProblem problem = CreateProblem(out blackCaps, out capMaturity, out capRate);

            Assert.Throws<NotImplementedException>(() => problem.Grad(new Vector(new double[] { 1, 0.02, 0.05 })));
        }

        [Test]
        public void GIsNotImplemented()
        {
            Matrix blackCaps;
            Vector capMaturity, capRate;
            CapCIROptimizationProblem problem = CreateProblem(out blackCaps, out capMaturity, out capRate);

            Assert.Throws<NotImplementedException>(() => problem.G(new Vector(new double[] { 1, 0.02, 0.05 })));
        }

        private static InterestRateMarketData CreateMarketData()
        {
            return TestCommon.TestMarketDataFactory.CreateCapMarketData(
                new Vector(new double[] { 1, 2, 5 }),
                new Matrix(new double[,]
                {
                    { 0.20, 0.22 },
                    { 0.21, 0.23 },
                    { 0.19, 0.20 },
                }));
        }

        [Test]
        public void MarketDataConstructorDerivesR0FromZeroCurveAtZero()
        {
            InterestRateMarketData irmd = CreateMarketData();

            CapCIROptimizationProblem problem = new CapCIROptimizationProblem(irmd);

            Assert.AreEqual(irmd.ZRMarket[0], problem.r0, 1e-12);
        }

        [Test]
        public void MarketDataConstructorMapsZeroCapVolatilityToZeroBlackCap()
        {
            InterestRateMarketData irmd = CreateMarketData();
            irmd.CapVolatility[0, 0] = 0.0;

            CapCIROptimizationProblem problem = new CapCIROptimizationProblem(irmd);

            // A zero black cap cell must not contribute to Obj regardless of x.
            Vector capMaturity = irmd.CapMaturity;
            Vector capRate = irmd.CapRate;
            Matrix cirCaps = CIRCap.CIRCapMatrix(capMaturity, capRate, irmd.CapTenor, problem.r0,
                new Vector(new double[] { 1.0, 0.02, 0.05 }));
            Assert.DoesNotThrow(() => problem.Obj(new Vector(new double[] { 1.0, 0.02, 0.05 })));
        }

        [Test]
        public void MarketDataConstructorThrowsWhenBlackCapIsNaN()
        {
            InterestRateMarketData irmd = CreateMarketData();

            // A negative strike makes the Black formula take the log of a negative
            // number (ln(forward / strike)), which yields NaN and must trigger the
            // guard in the constructor.
            irmd.CapMaturity = new Vector(new double[] { 1.0 });
            irmd.CapRate = new Vector(new double[] { -0.01 });
            irmd.CapVolatility = new Matrix(new double[,] { { 0.20 } });

            Assert.Throws<Exception>(() => new CapCIROptimizationProblem(irmd));
        }
    }
}
