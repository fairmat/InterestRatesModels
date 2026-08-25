using System;
using System.Collections.Generic;
using DVPLDOM;
using DVPLI;
using DVPLUtils;
using NUnit.Framework;

namespace HullAndWhiteOneFactor
{
    [TestFixture]
    public class TestHWCompactSimulator
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        private const int I = 12 * 40;

        private static Function CreateFlatZeroCurve(double rate)
        {
            Function zeroratecurve = new PFunction(null);
            zeroratecurve.Expr = new double[,] { { 0, rate }, { 50, rate } };
            (zeroratecurve as PFunction).m_Function.iType = EInterpolationType.LINEAR;
            return zeroratecurve;
        }

        // Mirrors the private Theta/F/Ft chain in HWCompactSimulator, operating on
        // the same zr/a/sigma, so ExpectedShortRate/ExpectedAverageRate can be
        // cross-checked without reflecting into the class's protected/private members.
        private static double F(IFunction zr, double t, double dt)
        {
            double zrT = zr.Evaluate(t);
            return t * (zr.Evaluate(t + dt) - zrT) / dt + zrT;
        }

        private static double Ft(IFunction zr, double t, double dt)
        {
            return (F(zr, t + dt, dt) - F(zr, t - dt, dt)) / (2 * dt);
        }

        private static double Theta(IFunction zr, double a, double sigma, double t, double dt)
        {
            return Ft(zr, t, dt) + a * F(zr, t, dt) + (1.0 - Math.Exp(-2.0 * a * t)) * sigma * sigma / (2.0 * a);
        }

        private static double ExpectedShortRateReplica(IFunction zr, double a, double sigma, double t)
        {
            double term1 = Math.Exp(-a * t) * zr.Evaluate(0);
            double ds = 0.001;
            double term2 = 0;
            for (double s = 0; s <= t; s += ds)
                term2 += Math.Exp(a * (s - t)) * Theta(zr, a, sigma, s, ds) * ds;
            return term1 + term2;
        }

        [Test]
        public void ExpectedShortRateMatchesManualReplicaOnFlatCurve()
        {
            Function zr = CreateFlatZeroCurve(0.02);
            HWCompactSimulator sim = new HWCompactSimulator();
            sim.zr = zr;
            sim.a = 0.1;
            sim.sigma = 0.01;

            double t = 0.5;
            double actual = sim.ExpectedShortRate(t);
            double expected = ExpectedShortRateReplica(zr, 0.1, 0.01, t);

            Assert.AreEqual(expected, actual, 1e-9);
        }

        [Test]
        public void ExpectedAverageRateReturnsFiniteValue()
        {
            Function zr = CreateFlatZeroCurve(0.02);
            HWCompactSimulator sim = new HWCompactSimulator();
            sim.zr = zr;
            sim.a = 0.1;
            sim.sigma = 0.01;

            double result = sim.ExpectedAverageRate(0.5);

            Assert.IsTrue(double.IsFinite(result));
        }

        [Test]
        public void SimulateProducesExpectedNumberOfDatesWithCorrectSpacing()
        {
            Function zr = CreateFlatZeroCurve(0.02);
            HWCompactSimulator sim = new HWCompactSimulator();
            sim.zr = zr;
            sim.a = 0.1;
            sim.sigma = 0.01;

            double m = 5.0;
            List<double> dates, finalRate, avgRate;
            sim.Simulate(m, out dates, out finalRate, out avgRate);

            Assert.AreEqual(I, dates.Count);
            Assert.AreEqual(I, finalRate.Count);
            Assert.AreEqual(I, avgRate.Count);

            double dt = m / I;
            for (int i = 0; i < dates.Count; i++)
                Assert.AreEqual(i * dt, dates[i], 1e-12);
        }
    }
}
