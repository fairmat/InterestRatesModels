/* Copyright (C) 2009-2012 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
 * Author(s): Enrico Degiuli (enrico.degiuli@fairmat.com)
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
using DVPLI;
using Fairmat.Statistics;

namespace CIRProcess
{
    /// <summary>
    /// Calculation of Cap price in the CIR model.
    /// </summary>
    public static class CIRCap
    {
        /// <summary>
        /// Evaluates a cap matrix using specified parameters of CIR model.
        /// </summary>
        /// <param name="maturities">Cap Matrix maturities.</param>
        /// <param name="strikes">Cap Matrix strikes.</param>
        /// <param name="tau">Cap tenor (inverse of number of caplets per year.</param>
        /// <param name="r0">ZR evaluated in 0.</param>
        /// <param name="cirParameters">Vector containing k, theta, sigma.</param>
        /// <returns>The cap matrix.</returns>
        public static Matrix CIRCapMatrix(Vector maturities, Vector strikes, double tau, double r0, Vector cirParameters)
        {
            int numCapletYear = (int)Math.Round(1.0 / tau);
            int nStrikes = strikes.Length;
            Vector numPrime = 1 + strikes * tau;
            Vector xPrime = new Vector(numPrime.Length);
            for (int i = 0; i < numPrime.Length; i++)
                xPrime[i] = 1.0 / numPrime[i];
            Vector matAll = Vector.Linspace(tau, maturities[Range.End],
                                            ((int)maturities[Range.End]) * numCapletYear);
            int npay = matAll.Length - 1;
            Matrix cirCapletMatrix = new Matrix(npay, nStrikes);
            for (int i = 0; i < nStrikes; i++)
            {
                for (int j = 0; j < npay; j++)
                {
                    cirCapletMatrix[j, i] = System.Math.Max(numPrime[i] * BondPut(r0, 0.0, matAll[j], matAll[j + 1], xPrime[i], cirParameters), 0.0);
                }
            }

            Matrix cirCapletSum = new Matrix(npay, nStrikes);
            cirCapletSum[0, Range.All] = cirCapletMatrix[0, Range.All];

            for (int i = 0; i < nStrikes; i++)
                cirCapletSum[Range.All, i] = ((Vector)cirCapletMatrix[Range.All, i]).CumSum();
            Matrix cirCap = new Matrix(maturities.Length, strikes.Length);
            int[] index = new int[maturities.Length];
            for (int i = 0; i < maturities.Length; i++)
                index[i] = ((int)Math.Round(maturities[i])) * numCapletYear - 2;
            for (int i = 0; i < maturities.Length; i++)
                cirCap[i, Range.All] = cirCapletSum[index[i], Range.All];
            return cirCap;
        }

        /// <summary>
        /// Calculates a put option on a bond.
        /// </summary>
        /// <param name="r0">The Zr evaluated in zero.</param>
        /// <param name="t">The valuation date (in years fractions).</param>
        /// <param name="T">The option exercise time to buy a bond with maturity S.</param>
        /// <param name="S">The maturity of the bound to buy.</param>
        /// <param name="X">The strike price.</param>
        /// <param name="cirParameters">
        /// A vector containing CIR model parameters k, theta, sigma.
        /// </param>
        /// <returns>The put option price.</returns>
        public static double BondPut(double r0, double t, double T, double S, double X, Vector cirParameters)
        {
            double k = cirParameters[0];
            double theta = cirParameters[1];
            double sigma = cirParameters[2];

            double h = Math.Sqrt(k * k + 2.0 * sigma * sigma);

            double dentS = 2.0 * h + (k + h) * (Math.Exp(h * (S - t)) - 1.0);
            double AtS = Math.Pow(2.0 * h * Math.Exp(0.5 * (k + h) * (S - t)) / dentS, 2.0 * k * theta / (sigma * sigma));
            double BtS = 2.0 * (Math.Exp(h * (S - t)) - 1.0) / dentS;

            double dentT = 2.0 * h + (k + h) * (Math.Exp(h * (T - t)) - 1.0);
            double AtT = Math.Pow(2.0 * h * Math.Exp(0.5 * (k + h) * (T - t)) / dentT, 2.0 * k * theta / (sigma * sigma));
            double BtT = 2.0 * (Math.Exp(h * (T - t)) - 1.0) / dentT;

            double PtS = AtS * Math.Exp(-BtS * r0);
            double PtT = AtT * Math.Exp(-BtT * r0);
            return BondCall(r0, t, T, S, X, cirParameters) - PtS + X * PtT;
        }

        /// <summary>
        /// Calculates a call option on a bond.
        /// </summary>
        /// <param name="r0">The Zr evaluated in zero.</param>
        /// <param name="t">The valuation date (in years fractions).</param>
        /// <param name="T">The option exercise time to buy a bond with maturity S.</param>
        /// <param name="S">The maturity of the bound to buy.</param>
        /// <param name="X">The strike price.</param>
        /// <param name="cirParameters">
        /// A vector containing CIR model parameters k, theta, sigma.
        /// </param>
        /// <returns>The call option price.</returns>
        public static double BondCall(double r0, double t, double T, double S, double X, Vector cirParameters)
        {
            double k = cirParameters[0];
            double theta = cirParameters[1];
            double sigma = cirParameters[2];

            double h = Math.Sqrt(k * k + 2.0 * sigma * sigma);

            double dentS = 2.0 * h + (k + h) * (Math.Exp(h * (S - t)) - 1.0);
            double AtS = Math.Pow(2.0 * h * Math.Exp(0.5 * (k + h) * (S - t)) / dentS, 2.0 * k * theta / (sigma * sigma));
            double BtS = 2.0 * (Math.Exp(h * (S - t)) - 1.0) / dentS;

            double denTS = 2.0 * h + (k + h) * (Math.Exp(h * (S - T)) - 1.0);
            double ATS = Math.Pow(2.0 * h * Math.Exp(0.5 * (k + h) * (S - T)) / denTS, 2.0 * k * theta / (sigma * sigma));
            double BTS = 2.0 * (Math.Exp(h * (S - T)) - 1.0) / denTS;

            double dentT = 2.0 * h + (k + h) * (Math.Exp(h * (T - t)) - 1.0);
            double AtT = Math.Pow(2.0 * h * Math.Exp(0.5 * (k + h) * (T - t)) / dentT, 2.0 * k * theta / (sigma * sigma));
            double BtT = 2.0 * (Math.Exp(h * (T - t)) - 1.0) / dentT;

            double PtS = AtS * Math.Exp(-BtS * r0);
            double PtT = AtT * Math.Exp(-BtT * r0);

            double rstar = Math.Log(ATS / X) / BTS;
            double rho = 2.0 * h / (sigma * sigma * (Math.Exp(h * (T - t)) - 1.0));
            double psi = (k + h) / (sigma * sigma);

            double v1 = 4.0 * k * theta / (sigma * sigma);
            double v2 = 2.0 * rho * rho * r0 * Math.Exp(h * (T - t)) / (rho + psi + BTS);
            double v3 = 2.0 * rho * rho * r0 * Math.Exp(h * (T - t)) / (rho + psi);

            double x1 = 2.0 * rstar * (rho + psi + BTS);
            double x2 = 2.0 * rstar * (rho + psi);

            double c1 = v1 / 0.5 + Math.Floor(v2 / 0.5) + 1;
            double c2 = v1 / 0.5 + Math.Floor(v3 / 0.5) + 1;

            NonCentralChiSquare ncx21 = new NonCentralChiSquare(v1, v2);
            NonCentralChiSquare ncx22 = new NonCentralChiSquare(v1, v3);

            return PtS * ncx21.Cdf(x1) - X * PtT * ncx22.Cdf(x2);
        }
    }
}
