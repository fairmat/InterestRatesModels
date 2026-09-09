using System;
using DVPLDOM;
using DVPLI;
using DVPLUtils;
using Fairmat.Math;
using NUnit.Framework;

namespace HullAndWhiteOneFactor
{
    [TestFixture]
    public class TestSwaptionHW1
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        private static Function CreateFlatZeroCurve(double rate)
        {
            return TestCommon.TestCurveFactory.CreateFlatZeroCurve(rate);
        }

        [Test]
        public void ForwardSwapRateMatchesManualFormulaOnFlatCurve()
        {
            double rate = 0.02;
            Function zr = CreateFlatZeroCurve(rate);
            SwaptionHW1 shw1 = new SwaptionHW1(zr);

            double t = 1.0;
            Vector T = new Vector(new double[] { 1.5, 2.0, 2.5 });

            double actual = shw1.ForwardSwapRate(t, T);

            Func<double, double> pzc = time => Math.Exp(-rate * time);
            double den = (T[0] - t) * pzc(T[0]);
            for (int i = 0; i < T.Length - 1; i++)
                den += (T[i + 1] - T[i]) * pzc(T[i + 1]);
            double expected = (pzc(t) - pzc(T[T.Length - 1])) / den;

            Assert.AreEqual(expected, actual, 1e-12);
        }

        [Test]
        public void AMatchesManualClosedFormFormula()
        {
            double rate = 0.02;
            Function zr = CreateFlatZeroCurve(rate);
            double t = 1.0, T = 3.0, alpha = 0.1, sigma = 0.01;

            double actual = SwaptionHW1.A(t, T, alpha, sigma, zr);

            double dT = T - t;
            double firstTerm = sigma * sigma * (alpha * dT - 2.0 * (1.0 - Math.Exp(-alpha * dT))
                + 0.5 * (1.0 - Math.Exp(-2.0 * alpha * dT))) / (2.0 * Math.Pow(alpha, 3.0));
            double alphaIntFirst = zr.Evaluate(T) * T - zr.Evaluate(t) * t;
            double alphaInt = alphaIntFirst + sigma * sigma * (alpha * (T - t) - 2.0 * (Math.Exp(-alpha * t) - Math.Exp(-alpha * T))
                + 0.5 * (Math.Exp(-2.0 * alpha * t) - Math.Exp(-2.0 * alpha * T))) / (2.0 * Math.Pow(alpha, 3.0));
            double expected = firstTerm - alphaInt;

            Assert.AreEqual(expected, actual, 1e-12);
        }

        [Test]
        public void HWSwaptionMatrixCellMatchesDirectHWSwaptionCall()
        {
            Function zr = CreateFlatZeroCurve(0.02);
            SwaptionHW1 shw1 = new SwaptionHW1(zr);

            Vector swaptionMaturity = new Vector(new double[] { 1.0, 2.0 });
            Vector swapDuration = new Vector(new double[] { 2.0 });
            double a = 0.1, sigma = 0.01, deltaK = 0.5;

            Matrix result = shw1.HWSwaptionMatrix(swaptionMaturity, swapDuration, a, sigma, deltaK);

            for (int i = 0; i < swaptionMaturity.Length; i++)
            {
                int npayment = (int)(swapDuration[0] / deltaK);
                Vector swapPayDate = Vector.Linspace(swaptionMaturity[i] + deltaK,
                    swaptionMaturity[i] + swapDuration[0], npayment);
                double fsr = shw1.ForwardSwapRate(swaptionMaturity[i], swapPayDate);
                double expected = shw1.HWSwaption(a, sigma, 1000.0, fsr, swaptionMaturity[i], swapPayDate);

                Assert.AreEqual(expected, result[i, 0], 1e-9);
            }
        }

        [Test]
        public void HWSwaptionWithDeltaKOverloadMatchesExplicitPayDateOverload()
        {
            Function zr = CreateFlatZeroCurve(0.02);
            SwaptionHW1 shw1 = new SwaptionHW1(zr);

            double a = 0.1, sigma = 0.01, l = 1000.0, k = 0.02, maturity = 1.0, duration = 2.0, deltaK = 0.5;

            double actual = shw1.HWSwaption(a, sigma, l, k, maturity, duration, deltaK);

            int npayment = (int)(duration / deltaK);
            Vector swapPayDate = Vector.Linspace(maturity + deltaK, maturity + duration, npayment);
            double expected = shw1.HWSwaption(a, sigma, l, k, maturity, swapPayDate);

            Assert.AreEqual(expected, actual, 1e-12);
        }

        [Test]
        public void HWSwaptionReturnsFiniteNonNegativePrice()
        {
            Function zr = CreateFlatZeroCurve(0.02);
            SwaptionHW1 shw1 = new SwaptionHW1(zr);

            double price = shw1.HWSwaption(0.1, 0.01, 1000.0, 0.02, 1.0, 2.0, 0.5);

            Assert.IsTrue(double.IsFinite(price));
            Assert.GreaterOrEqual(price, 0.0);
        }

        [Test]
        public void ZCBPutMatchesManualFormula()
        {
            Function zr = CreateFlatZeroCurve(0.02);
            SwaptionHW1 shw1 = new SwaptionHW1(zr);

            double a = 0.1, sigma = 0.01, L = 1.0, K = 0.95, T = 1.0, s = 2.0;

            double actual = shw1.ZCBPut(a, sigma, L, K, T, s);

            double sigmap = sigma * (1.0 - Math.Exp(-a * (s - T))) * Math.Sqrt(0.5 * (1.0 - Math.Exp(-2.0 * a * T)) / a) / a;
            Func<double, double> pzc = t => Math.Exp(-zr.Evaluate(t) * t);
            double h = Math.Log(L * pzc(s) / (K * pzc(T))) / sigmap + 0.5 * sigmap;
            double expected = K * pzc(T) * SpecialFunctions.NormCdf(-h + sigmap) - L * pzc(s) * SpecialFunctions.NormCdf(-h);

            Assert.AreEqual(expected, actual, 1e-12);
        }

        [Test]
        public void AlphaTFuncMatchesManualFormulaAfterHWSwaptionSetsState()
        {
            Function zr = CreateFlatZeroCurve(0.02);
            SwaptionHW1 shw1 = new SwaptionHW1(zr);

            double a = 0.1, sigma = 0.01;
            // HWSwaption sets the private a/sigma fields used by alphaTFunc as a side effect.
            shw1.HWSwaption(a, sigma, 1000.0, 0.02, 1.0, 2.0, 0.5);

            double t = 1.5;
            double actual = shw1.alphaTFunc(t);

            double dt = 0.001;
            double f = (zr.Evaluate(t + dt) * (t + dt) - zr.Evaluate(t) * t) / dt;
            double expected = f + sigma * sigma * Math.Pow(1.0 - Math.Exp(-a * t), 2.0) / (2.0 * a * a);

            Assert.AreEqual(expected, actual, 1e-9);
        }
    }
}
