/* Copyright (C) 2009-2012 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
 * Author(s): Michele Furgeri (info@fairmat.com)
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
using System.Collections.Generic;
using System.Text;
using DVPLDOM;
using DVPLI;
using Fairmat.Math;

namespace HullAndWhiteOneFactor
{
   

    /// <summary>
    /// This class implements the HW1 model
    /// for pricing of a Cap (caplets portfolio)
    /// and a Floor (floorlets portfolio).
    /// </summary>
    public class CapHW1
    {
        /// <summary>
        /// A zero rate function used for this object evaluations.
        /// </summary>
        internal Function zeroRateCurve;

        /// <summary>
        /// Initializes a new instance of the <see cref="HullAndWhiteOneFactor.CapHW1"/> class.
        /// </summary>
        /// <param name='zeroratecurve'>
        /// The zero-rate curve reference for the model.
        /// </param>
        public CapHW1(Function zeroratecurve)
        {
            this.zeroRateCurve = zeroratecurve;
        }

        /// <summary>
        /// Calculates a matrix of cap prices within the Hull-White model.
        /// </summary>
        /// <param name='capMaturity'>
        /// Vector of cap maturities.
        /// </param>
        /// <param name='capRate'>
        /// Vector of cap strikes.
        /// </param>
        /// <param name='a'>
        /// Hull-White alpha parameter.
        /// </param>
        /// <param name='sigma'>
        /// Hull-White sigma parameter.
        /// </param>
        /// <param name='deltaK'>
        /// Time period between caplets expressed in year fraction.
        /// </param>
        /// <returns>
        /// The cap prices matrix.
        /// </returns>
        public Matrix HWMatrixCaps(Vector capMaturity, Vector capRate, double a, double sigma, double deltaK)
        {
            Matrix hwCapsMatrix = new Matrix(capMaturity.Length, capRate.Length);
            for (int m = 0; m < capMaturity.Length; m++)
            {
                for (int s = 0; s < capRate.Length; s++)
                {
                    hwCapsMatrix[m, s] = HWCap(a, sigma, capRate[s], deltaK, capMaturity[m]);
                }
            }

            return hwCapsMatrix;
        }

        /// <summary>
        /// Calculates the price of a cap within the Hull-White model.
        /// </summary>
        /// <param name='a'>
        /// Hull-White alpha parameter.
        /// </param>
        /// <param name='sigma'>
        /// Hull-White sigma parameter.
        /// </param>
        /// <param name='K'>
        /// Strike rate.
        /// </param>
        /// <param name='deltaK'>
        /// Time period between caplets expressed in year fraction.
        /// </param>
        /// <param name='T'>
        /// Cap maturity.
        /// </param>
        /// <returns>
        /// The cap price.
        /// </returns>
        public double HWCap(double a, double sigma, double K, double deltaK, double T)
        {
            double KK = 1.0 / (1.0 + K * deltaK);
            double NP = Math.Round(T / deltaK);
            int np = (int)NP;
            double CA = 0;
            double[] s = new double[np];
            s[0] = deltaK;
            double[] caplets = new double[np - 1];
            for (int i = 1; i < s.Length; i++)
            {
                s[i] = s[i - 1] + deltaK;
                caplets[i - 1] = HWCaplet(a, sigma, s[i - 1], s[i], KK);
                CA = CA + caplets[i - 1];
            }

            return CA;
        }

        /// <summary>
        /// Calculates the price of a floor within the Hull-White model.
        /// </summary>
        /// <param name='a'>
        /// Hull-White alpha parameter.
        /// </param>
        /// <param name='sigma'>
        /// Hull-White sigma parameter.
        /// </param>
        /// <param name='K'>
        /// Strike rate.
        /// </param>
        /// <param name='deltaK'>
        /// Time period between floorlets expressed in year fraction.
        /// </param>
        /// <param name='T'>
        /// Floor maturity.
        /// </param>
        /// <returns>
        /// The floor price.
        /// </returns>
        private double HWFloor(double a, double sigma, double K, double deltaK, double T)
        {
            double KK = 1.0 / (1.0 + K * deltaK);
            double NP = Math.Round(T / deltaK);
            int np = (int)NP;
            double FL = 0;
            double[] s = new double[np];
            s[0] = deltaK;
            double[] floorlets = new double[np - 1];
            for (int i = 1; i < s.Length; i++)
            {
                s[i] = s[i - 1] + deltaK;
                floorlets[i - 1] = HWFloorlet(a, sigma, s[i - 1], s[i], KK);
                FL = FL + floorlets[i - 1];
            }

            return FL;
        }

        /// <summary>
        /// Calculates the price of a single caplet within the Hull-White model.
        /// </summary>
        /// <param name='a'>
        /// Hull-White alpha parameter.
        /// </param>
        /// <param name='sigma'>
        /// Hull-White sigma parameter.
        /// </param>
        /// <param name='T'>
        /// Caplet reset date.
        /// </param>
        /// <param name='s'>
        /// Caplet payment date.
        /// </param>
        /// <param name='K'>
        /// Strike rate.
        /// </param>
        /// <returns>
        /// The caplet price.
        /// </returns>
        public double HWCaplet(double a, double sigma, double T, double s, double K)
        {
            double d1 = D1(a, sigma, T, s, K);
            double d2 = D2(a, sigma, T, s, K);
            return (K * PZC(T) * SpecialFunctions.NormCdf(-d2)) - (PZC(s) * SpecialFunctions.NormCdf(-d1));
        }

        /// <summary>
        /// Calculates the price of a single floorlet within the Hull-White model.
        /// </summary>
        /// <param name='a'>
        /// Hull-White alpha parameter.
        /// </param>
        /// <param name='sigma'>
        /// Hull-White sigma parameter.
        /// </param>
        /// <param name='T'>
        /// Floorlet reset date.
        /// </param>
        /// <param name='s'>
        /// Floorlet payment date.
        /// </param>
        /// <param name='K'>
        /// Strike rate.
        /// </param>
        /// <returns>
        /// The floorlet price.
        /// </returns>
        private double HWFloorlet(double a, double sigma, double T, double s, double K)
        {
            double d1 = D1(a, sigma, T, s, K);
            double d2 = D2(a, sigma, T, s, K);
            return (PZC(s) * SpecialFunctions.NormCdf(d1)) - (K * PZC(T) * SpecialFunctions.NormCdf(d2));
        }

        /// <summary>
        /// Calculates d1 factor to be used in pricing formulas.
        /// </summary>
        /// <param name='a'>
        /// Hull-White alpha parameter.
        /// </param>
        /// <param name='sigma'>
        /// Hull-White sigma parameter.
        /// </param>
        /// <param name='T'>
        /// Reset date.
        /// </param>
        /// <param name='s'>
        /// Maturity date.
        /// </param>
        /// <param name='K'>
        /// Strike rate.
        /// </param>
        /// <returns>
        /// The value of the d1 factor.
        /// </returns>
        public double D1(double a, double sigma, double T, double s, double K)
        {
            return (1.0 / SigP(a, sigma, T, s)) * Math.Log(PZC(s) / (PZC(T) * K)) + (SigP(a, sigma, T, s) / 2.0);
        }

        /// <summary>
        /// Calculates d2 factor to be used in pricing formulas.
        /// </summary>
        /// <param name='a'>
        /// Hull-White alpha parameter.
        /// </param>
        /// <param name='sigma'>
        /// Hull-White sigma parameter.
        /// </param>
        /// <param name='T'>
        /// Reset date.
        /// </param>
        /// <param name='s'>
        /// Maturity date.
        /// </param>
        /// <param name='K'>
        /// Strike rate.
        /// </param>
        /// <returns>
        /// The value of the d2 factor.
        /// </returns>
        public double D2(double a, double sigma, double T, double s, double K)
        {
            return D1(a, sigma, T, s, K) - SigP(a, sigma, T, s);
        }

        /// <summary>
        /// Calculate the volatility term to be used in Hull-White pricing formula.
        /// </summary>
        /// <param name='a'>
        /// Hull-White alpha parameter.
        /// </param>
        /// <param name='sigma'>
        /// Hull-White sigma parameter.
        /// </param>
        /// <param name='T'>
        /// Reset date.
        /// </param>
        /// <param name='s'>
        /// Strike rate.
        /// </param>
        /// <returns>
        /// The volatility.
        /// </returns>
        private double SigP(double a, double sigma, double T, double s)
        {
            return (sigma / a) * (1.0 - Math.Exp(-a * (s - T))) * Math.Sqrt((1.0 - Math.Exp(-2.0 * a * T)) / (2.0 * a));
        }

        /// <summary>
        /// Helper function to make functions easier to read.
        /// Just returns the value of the zero rate at position k.
        /// </summary>
        /// <param name="k">The position where to get the value of the zero rate from.</param>
        /// <returns>The value of the zero rate at position k.</returns>
        internal double ZR(double k)
        {
            return this.zeroRateCurve.Evaluate(k);
        }

        /// <summary>
        /// Helper function to make functions easier to read.
        /// Just returns the value of the discount factor at position t.
        /// This is calculated with e^(-ZR(t)*t).
        /// </summary>
        /// <param name="t">The position where to get the value of discount factor from.</param>
        /// <returns>The value of the discount factor at position t.</returns>
        private double PZC(double t)
        {
            if (t == 0)
                return 1.0;
            else
                return Math.Exp(-ZR(t) * t);
        }
    }
}
