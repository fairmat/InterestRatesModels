using System;
using NUnit.Framework;

namespace HullAndWhiteTwoFactors
{
    /// <summary>
    /// Exposes the protected internal members of <see cref="HW2Context"/> for testing.
    /// </summary>
    internal class TestableHW2Context : HW2Context
    {
        public double Alpha1
        {
            get { return this.alpha1; }
            set { this.alpha1 = value; }
        }

        public double Alpha2
        {
            get { return this.alpha2; }
            set { this.alpha2 = value; }
        }

        public double Sigma1
        {
            get { return this.sigma1; }
            set { this.sigma1 = value; }
        }

        public double Sigma2
        {
            get { return this.sigma2; }
            set { this.sigma2 = value; }
        }

        public double Rho
        {
            get { return this.rho; }
            set { this.rho = value; }
        }

        /// <summary>
        /// Calls the protected internal <see cref="HW2Context.Eta"/> method.
        /// </summary>
        /// <param name="t">The valuation date.</param>
        /// <param name="s">The maturity date.</param>
        /// <returns>The result of the Eta function.</returns>
        public double CallEta(double t, double s)
        {
            return this.Eta(t, s);
        }

        /// <summary>
        /// Calls the protected internal <see cref="HW2Context.BHat"/> method.
        /// </summary>
        /// <param name="t">The valuation date.</param>
        /// <param name="s">The maturity.</param>
        /// <param name="dt">The delta between this t position and the previous one.</param>
        /// <returns>The result of the BHat function.</returns>
        public double CallBHat(double t, double s, double dt)
        {
            return this.BHat(t, s, dt);
        }

        /// <summary>
        /// Calls the protected internal <see cref="HW2Context.CHat"/> method.
        /// </summary>
        /// <param name="t">The valuation date.</param>
        /// <param name="s">The maturity.</param>
        /// <param name="dt">The difference between one approximated time step and the previous one.</param>
        /// <param name="bHat">The pre-calculated value of the BHat method.</param>
        /// <returns>The result of the CHat function.</returns>
        public double CallCHat(double t, double s, double dt, double bHat)
        {
            return this.CHat(t, s, dt, bHat);
        }

        /// <summary>
        /// Calls the protected internal <see cref="HW2Context.Chat"/> method.
        /// </summary>
        /// <param name="t">The valuation date.</param>
        /// <param name="s">The maturity.</param>
        /// <param name="dt">The difference between one approximated time step and the previous one.</param>
        /// <returns>The result of the Chat function.</returns>
        public double CallChat(double t, double s, double dt)
        {
            return this.Chat(t, s, dt);
        }
    }

    [TestFixture]
    public class TestHW2Context
    {
        /// <summary>
        /// Performs the common test suite initialization before each test.
        /// </summary>
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        /// <summary>
        /// Verifies that Eta matches its closed-form value for a realistic set of parameters.
        /// </summary>
        [Test]
        public void Eta_RealisticParameters_MatchesClosedForm()
        {
            TestableHW2Context context = new TestableHW2Context();
            context.Alpha1 = 1.0;
            context.Alpha2 = 0.1;
            context.Sigma1 = 0.01;
            context.Sigma2 = 0.0165;
            context.Rho = 0.6;

            double result = context.CallEta(5.0, 7.0);

            Assert.That(result, Is.EqualTo(0.003282050022501).Within(1e-9));
        }

        /// <summary>
        /// Verifies that Eta matches its closed-form value for a second, distinct set of parameters.
        /// </summary>
        [Test]
        public void Eta_AlternateParameterSet_MatchesClosedForm()
        {
            TestableHW2Context context = new TestableHW2Context();
            context.Alpha1 = 0.1;
            context.Alpha2 = 0.05;
            context.Sigma1 = 0.02;
            context.Sigma2 = 0.01;
            context.Rho = 0.3;

            double result = context.CallEta(2.0, 5.0);

            Assert.That(result, Is.EqualTo(0.008539763709306).Within(1e-9));
        }

        /// <summary>
        /// Verifies that Eta is exactly zero when the valuation date equals t = 0.
        /// </summary>
        [Test]
        public void Eta_AtValuationTime_IsExactlyZero()
        {
            TestableHW2Context context = new TestableHW2Context();
            context.Alpha1 = 0.1;
            context.Alpha2 = 0.05;
            context.Sigma1 = 0.02;
            context.Sigma2 = 0.01;
            context.Rho = 0.3;

            double result = context.CallEta(0.0, 3.0);

            Assert.That(result, Is.EqualTo(0.0).Within(1e-12));
        }

        /// <summary>
        /// Verifies that BHat matches its closed-form value for a realistic set of parameters.
        /// </summary>
        [Test]
        public void BHat_RealisticParameters_MatchesClosedForm()
        {
            TestableHW2Context context = new TestableHW2Context();
            context.Alpha1 = 1.0;
            context.Alpha2 = 0.1;

            double result = context.CallBHat(5.0, 7.0, 5.0 / 512.0);

            Assert.That(result, Is.EqualTo(0.868893584177557).Within(1e-9));
        }

        /// <summary>
        /// Verifies that BHat matches its closed-form value for a second, distinct set of parameters.
        /// </summary>
        [Test]
        public void BHat_AlternateParameterSet_MatchesClosedForm()
        {
            TestableHW2Context context = new TestableHW2Context();
            context.Alpha1 = 0.1;
            context.Alpha2 = 0.05;

            double result = context.CallBHat(2.0, 5.0, 1.0);

            Assert.That(result, Is.EqualTo(2.723568171113940).Within(1e-9));
        }

        /// <summary>
        /// Verifies that CHat matches its closed-form value for a realistic set of parameters.
        /// </summary>
        [Test]
        public void CHat_RealisticParameters_MatchesClosedForm()
        {
            TestableHW2Context context = new TestableHW2Context();
            context.Alpha1 = 1.0;
            context.Alpha2 = 0.1;

            double bHat = context.CallBHat(5.0, 7.0, 5.0 / 512.0);
            double result = context.CallCHat(5.0, 7.0, 5.0 / 512.0, bHat);

            Assert.That(result, Is.EqualTo(1.049136679349510).Within(1e-9));
        }

        /// <summary>
        /// Verifies that CHat matches its closed-form value for a second, distinct set of parameters.
        /// </summary>
        [Test]
        public void CHat_AlternateParameterSet_MatchesClosedForm()
        {
            TestableHW2Context context = new TestableHW2Context();
            context.Alpha1 = 0.1;
            context.Alpha2 = 0.05;

            double bHat = context.CallBHat(2.0, 5.0, 1.0);
            double result = context.CallCHat(2.0, 5.0, 1.0, bHat);

            Assert.That(result, Is.EqualTo(2.584814583271093).Within(1e-9));
        }

        /// <summary>
        /// Verifies that Chat (which computes BHat internally) agrees with CHat given the same
        /// pre-computed BHat value.
        /// </summary>
        [Test]
        public void Chat_MatchesCHatGivenComputedBHat()
        {
            TestableHW2Context context = new TestableHW2Context();
            context.Alpha1 = 1.0;
            context.Alpha2 = 0.1;

            double t = 5.0;
            double s = 7.0;
            double dt = 5.0 / 512.0;
            double bHat = context.CallBHat(t, s, dt);

            double chatResult = context.CallChat(t, s, dt);
            double cHatResult = context.CallCHat(t, s, dt, bHat);

            Assert.That(chatResult, Is.EqualTo(cHatResult).Within(1e-12));
        }

        /// <summary>
        /// Verifies that Eta currently produces a non-finite (NaN or Infinity) result when
        /// alpha1 equals alpha2, since the closed-form expression divides by (alpha1 - alpha2).
        /// </summary>
        [Test]
        public void Eta_WhenAlphasAreEqual_ProducesNonFiniteResult()
        {
            TestableHW2Context context = new TestableHW2Context();
            context.Alpha1 = 0.1;
            context.Alpha2 = 0.1;
            context.Sigma1 = 0.01;
            context.Sigma2 = 0.0165;
            context.Rho = 0.6;

            double result = context.CallEta(5.0, 7.0);

            bool isNonFinite = double.IsNaN(result) || double.IsInfinity(result);
            Assert.That(isNonFinite, Is.True);
        }
    }
}
