using System;
using DVPLI;
using NUnit.Framework;

namespace CIRProcess
{
    /// <summary>
    /// Tests for <see cref="CIRCap.BondPut"/> and <see cref="CIRCap.CIRCapMatrix"/>.
    /// BondPut's expected value is derived independently from the CIR zero-coupon bond
    /// price formula (P(t,S) and P(t,T)), reusing the exact r0/k/theta/sigma/T/S/X used by
    /// TestBondCall.cs's already-verified BondCall benchmark, via put-call parity:
    /// BondPut = BondCall - P(t,S) + X*P(t,T).
    /// </summary>
    [TestFixture]
    public class TestCIRCap
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void BondPut_MatchesIndependentlyDerivedPutCallParityValue()
        {
            double k = 1.0;
            double theta = 0.02;
            double sigma = 0.08;
            double r0 = 0.01;
            Vector par = new Vector(3);
            par[0] = k;
            par[1] = theta;
            par[2] = sigma;

            double strike = 0.98;
            double noz = 100.0;

            double T = 1.0;
            double S = 2.0;

            double putFairmat = noz * CIRCap.BondPut(r0, 0.0, T, S, strike, par);
            double putBenchmark = 0.069892106145608;

            double maxError = 1e-9;
            Assert.Less(Math.Abs(putFairmat - putBenchmark), maxError);
        }

        [Test]
        public void CIRCapMatrix_FirstCellMatchesDirectBondPutCall()
        {
            double k = 1.0;
            double theta = 0.02;
            double sigma = 0.08;
            double r0 = 0.01;
            Vector par = new Vector(3);
            par[0] = k;
            par[1] = theta;
            par[2] = sigma;

            Vector maturities = (Vector)(new double[] { 2.0, 3.0 });
            Vector strikes = (Vector)(new double[] { 0.02 });
            double tau = 1.0;

            Matrix cap = CIRCap.CIRCapMatrix(maturities, strikes, tau, r0, par);

            double numPrime = 1.0 + strikes[0] * tau;
            double xPrime = 1.0 / numPrime;
            double expectedFirstCell = Math.Max(numPrime * CIRCap.BondPut(r0, 0.0, 1.0, 2.0, xPrime, par), 0.0);

            Assert.That(cap.R, Is.EqualTo(2));
            Assert.That(cap.C, Is.EqualTo(1));
            Assert.That(cap[0, 0], Is.EqualTo(expectedFirstCell).Within(1e-12));
        }

        [Test]
        public void CIRCapMatrix_IsNonDecreasingWithMaturity()
        {
            double k = 1.0;
            double theta = 0.02;
            double sigma = 0.08;
            double r0 = 0.01;
            Vector par = new Vector(3);
            par[0] = k;
            par[1] = theta;
            par[2] = sigma;

            Vector maturities = (Vector)(new double[] { 2.0, 3.0 });
            Vector strikes = (Vector)(new double[] { 0.02 });
            double tau = 1.0;

            Matrix cap = CIRCap.CIRCapMatrix(maturities, strikes, tau, r0, par);

            Assert.That(cap[1, 0], Is.GreaterThanOrEqualTo(cap[0, 0]));
        }
    }
}
