/* Copyright (C) 2009-2012 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
 * Author(s): Matteo Tesser (matteo.tesser@fairmat.com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

using System;

namespace HullAndWhiteTwoFactors
{
    /// <summary>
    /// Stores the 5 parameters of the Hull and White two factors model
    /// and implements several utility methods.
    /// </summary>
    /// <remarks>
    /// Fields are being used for performance reason and
    /// can only be used inside this specific assembly.
    /// </remarks>
    [Serializable]
    public class HW2Context
    {
        /// <summary>
        /// Mean-reversion rate for r.
        /// </summary>
        [NonSerialized]
        protected internal double alpha1;

        /// <summary>
        /// Mean-reversion rate for u.
        /// </summary>
        [NonSerialized]
        protected internal double alpha2;

        /// <summary>
        /// Volatility of r.
        /// </summary>
        [NonSerialized]
        protected internal double sigma1;

        /// <summary>
        /// Volatility for u.
        /// </summary>
        [NonSerialized]
        protected internal double sigma2;

        /// <summary>
        /// Correlation between r and u.
        /// </summary>
        [NonSerialized]
        protected internal double rho;

        /// <summary>
        /// Implements the function Eta. See  HW1994.
        /// </summary>
        /// <param name="t">The valuation date.</param>
        /// <param name="s">The maturity date.</param>
        /// <returns>
        /// The value of the function Eta depending on the
        /// current context and the valuation and maturity date.
        /// </returns>
        protected internal double Eta(double t, double s)
        {
            double a_2 = this.alpha1 * this.alpha1;
            double ab = this.alpha1 * this.alpha2;
            double a_p_b = this.alpha1 + this.alpha2;
            double a_a_b = +this.alpha1 * (this.alpha1 - this.alpha2);
            double b_a_b = +this.alpha2 * (this.alpha1 - this.alpha2);
            double Tt = s - t;
            double BtT = (1 - Math.Exp(-this.alpha1 * (Tt))) / this.alpha1;
            double B0T = (1 - Math.Exp(-this.alpha1 * s)) / this.alpha1;
            double B0t = (1 - Math.Exp(-this.alpha1 * t)) / this.alpha1;

            double CtT = Math.Exp(-this.alpha1 * Tt) / a_a_b - Math.Exp(-this.alpha2 * Tt) / b_a_b + 1 / ab;
            double C0T = Math.Exp(-this.alpha1 * s) / a_a_b - Math.Exp(-this.alpha2 * s) / b_a_b + 1 / ab;
            double C0t = Math.Exp(-this.alpha1 * t) / a_a_b - Math.Exp(-this.alpha2 * t) / b_a_b + 1 / ab;

            double g1 = (Math.Exp(-(a_p_b) * s) * (Math.Exp((a_p_b) * t) - 1)) / (/*eps2*/ +(a_p_b) * (this.alpha1 - this.alpha2)) - (Math.Exp(-2 * this.alpha1 * s) * (Math.Exp(2 * this.alpha1 * t) - 1)) / (2 * a_a_b);
            double g2 = (g1 + CtT - C0T + 0.5 * (BtT * BtT) - 0.5 * (B0T * B0T) + t / this.alpha1 - (Math.Exp(-this.alpha1 * (Tt)) - Math.Exp(-this.alpha1 * s)) / a_2) / (ab);
            double g3 = -(Math.Exp(-(a_p_b) * t) - 1) / ((this.alpha1 - this.alpha2) * (a_p_b)) + (Math.Exp(-2 * this.alpha1 * t) - 1) / (2 * a_a_b);
            double g4 = (g3 - C0t - 0.5 * (B0t * B0t) + t / this.alpha1 + (Math.Exp(-this.alpha1 * t) - 1) / a_2) / (ab);
            double g5 = (0.5 * (CtT * CtT) - 0.5 * (C0T * C0T) + g2) / this.alpha2;
            double g6 = (g4 - 0.5 * (C0t * C0t)) / this.alpha2;

            double e = (this.sigma1 * this.sigma1) / (4 * this.alpha1) *
                       (1 - Math.Exp(-2 * this.alpha1 * t)) * (BtT * BtT) -
                       this.rho * this.sigma1 * this.sigma2 * (B0t * C0t * BtT + g4 - g2) -
                       0.5 * (this.sigma2 * this.sigma2) * ((C0t * C0t) * BtT + g6 - g5);

            return e;
        }

        /// <summary>
        /// Bhat function, needed for bond calculation. See HW1994.
        /// </summary>
        /// <param name="t">The valuation date.</param>
        /// <param name="s">The maturity.</param>
        /// <param name="dt">The delta between this t position and the previous one.</param>
        /// <returns>The result of the Bhat function.</returns>
        protected internal double BHat(double t, double s, double dt)
        {
            double BtT = (1 - Math.Exp(-this.alpha1 * (s - t))) / this.alpha1;
            double Btdt = (1 - Math.Exp(-this.alpha1 * (dt))) / this.alpha1;
            return BtT / Btdt * dt;
        }

        /// <summary>
        /// Chat function, needed for bond calculation. See HW1994.
        /// </summary>
        /// <param name="t">The valuation date.</param>
        /// <param name="s">The maturity.</param>
        /// <param name="dt">
        /// The difference between one approximated time steps and the previous one.
        /// </param>
        /// <returns>The result of the Chat function.</returns>
        protected internal double Chat(double t, double s, double dt)
        {
            double a_a_b = this.alpha1 * (this.alpha1 - this.alpha2);

            double CtT = 1 / (a_a_b) * Math.Exp(-this.alpha1 * (s - t)) - 1 / (this.alpha2 * (this.alpha1 - this.alpha2)) * Math.Exp(-this.alpha2 * (s - t)) + 1 / (this.alpha1 * this.alpha2);
            double Ctdt = 1 / (a_a_b) * Math.Exp(-this.alpha1 * dt) - 1 / (this.alpha2 * (this.alpha1 - this.alpha2)) * Math.Exp(-this.alpha2 * dt) + 1 / (this.alpha1 * this.alpha2);

            return CtT - Ctdt * BHat(t, s, dt) / dt;
        }

        /// <summary>
        /// CHat function, needed for bond calculation. See HW1994.
        /// </summary>
        /// <param name="t">The valuation date.</param>
        /// <param name="s">The maturity.</param>
        /// <param name="dt">
        /// The difference between one approximated time steps and the previous one.
        /// </param>
        /// <param name="bHat">The pre-calculated value of method Bhat.</param>
        /// <returns>The result of the CHat function.</returns>
        protected internal double CHat(double t, double s, double dt, double bHat)
        {
            double a_a_b = this.alpha1 * (this.alpha1 - this.alpha2);
            double b_a_b = this.alpha2 * (this.alpha1 - this.alpha2);
            double CtT = 1 / (a_a_b) * Math.Exp(-this.alpha1 * (s - t)) - 1 / (b_a_b) * Math.Exp(-this.alpha2 * (s - t)) + 1 / (this.alpha1 * this.alpha2);
            double Ctdt = 1 / (a_a_b) * Math.Exp(-this.alpha1 * dt) - 1 / (b_a_b) * Math.Exp(-this.alpha2 * dt) + 1 / (this.alpha1 * this.alpha2);

            return CtT - Ctdt * bHat / dt;
        }
    }
}
