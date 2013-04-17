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
using System.Collections.Generic;
using System.Text;
using DVPLDOM;
using DVPLI;
using Fairmat.Optimization;
using Fairmat.Math;

namespace HullAndWhiteOneFactor
{
    /// <summary>
    /// This class provides closed form pricing of for Swaptions (the holder has the right to
    /// pay the fixed rate and receive floating rate) using the HW1 factor model.
    /// </summary>
    public class SwaptionHW1
    {
        // TODO: More expressive names needed.

        /// <summary>
        /// A zero rate function used for this object evaluations.
        /// </summary>
        private Function zeroRateCurve;

        /// <summary>
        /// Hull-White alpha parameter.
        /// </summary>
        private double a;

        /// <summary>
        /// Hull-White sigma parameter.
        /// </summary>
        private double sigma;

        /// <summary>
        /// Time at which calculate bond price in Func() method.
        /// </summary>
        private double t;

        /// <summary>
        /// Interval between bond coupons to be used in bond price calculation in Func() method.
        /// </summary>
        private double dt;

        /// <summary>
        /// Notional to be used in bond price calculation in Func() method.
        /// </summary>
        private double L;

        /// <summary>
        /// Vector of payment times to be used in bond price calculation in Func() method.
        /// </summary>
        private Vector T;

        /// <summary>
        /// Cash flow vector to be used in bond price calculation in Func() method.
        /// </summary>
        private Vector CF;

        /// <summary>
        /// Initializes a new instance of the <see cref="HullAndWhiteOneFactor.SwaptionHW1"/> class.
        /// </summary>
        /// <param name='zeroratecurve'>
        /// The zero-rate curve reference for the model.
        /// </param>
        public SwaptionHW1(Function zeroratecurve)
        {
            this.zeroRateCurve = zeroratecurve;
        }

        /// <summary>
        /// Calculates a matrix of swaption prices within Hull-White model.
        /// </summary>
        /// <param name='swaptionMaturity'>
        /// Vector of swaption maturity.
        /// </param>
        /// <param name='swapDuration'>
        /// Vector of swap duration.
        /// </param>
        /// <param name='a'>
        /// Hull-White alpha parameter.
        /// </param>
        /// <param name='sigma'>
        /// Hull-White sigma parameter.
        /// </param>
        /// <param name='deltaK'>
        /// Time interval between swap coupon expressed in year fraction.
        /// </param>
        /// <returns>
        /// The swaption price matrix.
        /// </returns>
        public Matrix HWSwaptionMatrix(Vector swaptionMaturity, Vector swapDuration, double a, double sigma, double deltaK)
        {
            Matrix result = new Matrix(swaptionMaturity.Length, swapDuration.Length);
            Vector swapPayDate;
            int npayment;
            double FSR;
            for (int i = 0; i < swaptionMaturity.Length; i++)
            {
                for (int j = 0; j < swapDuration.Length; j++)
                {
                    npayment = (int)(swapDuration[j] / deltaK);
                    swapPayDate = Vector.Linspace(swaptionMaturity[i] + deltaK, swaptionMaturity[i] + swapDuration[j], npayment);
                    FSR = this.ForwardSwapRate(swaptionMaturity[i], swapPayDate);
                    result[i, j] = this.HWSwaption(a, sigma, 1000.0, FSR, swaptionMaturity[i], swapPayDate);
                }
            }

            return result;
        }

        /// <summary>
        /// Calculates the forward swap rate.
        /// </summary>
        /// <param name='t'>
        /// The swap starting time.
        /// </param>
        /// <param name='T'>
        /// The vector of swap payment times.
        /// </param>
        /// <returns>
        /// The swap rate.
        /// </returns>
        public double ForwardSwapRate(double t, Vector T)
        {
            double den = (T[0] - t) * this.PZC(T[0]);
            for (int i = 0; i < (T.Length - 1); i++)
                den += (T[i + 1] - T[i]) * this.PZC(T[i + 1]);
            return (this.PZC(t) - this.PZC(T[T.Length - 1])) / den;
        }

        /// <summary>
        /// Price of a swaption in HW1 model: the holder has the right to pay
        /// the fixed rate and receive floating rate.
        /// This is a different version of the main one to do a double check of the swaption.
        /// </summary>
        /// <param name='a'>
        /// Hull-White alpha parameter.
        /// </param>
        /// <param name='sigma'>
        /// Hull-White sigma parameter.
        /// </param>
        /// <param name='l'>
        /// The Notional.
        /// </param>
        /// <param name='k'>
        /// The Strike.
        /// </param>
        /// <param name='T'>
        /// The maturity.
        /// </param>
        /// <param name='s'>
        /// Vector of swaption payment dates.
        /// </param>
        /// <returns>
        /// The swaptions price.
        /// </returns>
        public double HWSwaption(double a, double sigma, double l, double k, double T, Vector s)
        {
            OptimizationSettings options = new OptimizationSettings();
            options.epsilon = 1e-10;

            this.dt = s[0] - T;
            this.CF = k * this.dt + (new Vector(s.Length));
            this.CF[s.Length - 1] = this.CF[s.Length - 1] + 1;
            this.a = a;
            this.sigma = sigma;
            this.t = T;
            this.T = s;
            this.L = l;
            SolutionInfo sol = Fairmat.Optimization.Helper.FSolve(new ObjFunction(this.Func), (Vector)new double[1] { 0.01 },
                (Vector)new double[1] { -1.0 }, (Vector)new double[1] { 1.0 }, options);
            double RK = (double)sol.x;

            Vector X = this.HWBond(a, sigma, RK, this.T, this.t);
            double result = 0;
            for (int i = 0; i < s.Length; i++)
                result += this.CF[i] * this.ZCBPut(a, sigma, 1.0, X[i], T, s[i]);
            result = result * l;
            return result;
        }

        /// <summary>
        /// Function to be used in the FSolve() problem of HWSwaption() method.
        /// This is a version done to double check the HW Swaption, with dt fixed to 0.001.
        /// </summary>
        /// <param name='x'>
        /// Vector of length 1 representing the variable on which the FSolve problem is performed.
        /// </param>
        /// <returns>
        /// The value of the function given the provided vector x.
        /// </returns>
        public double Func(Vector x)
        {
            double result = this.CF.Scalar(this.HWBond(this.a, this.sigma, x[0], this.T, this.t));
            return result - 1.0;
        }

        /// <summary>
        /// Calculates a vector of zero coupon bond prices
        /// (discount factors) within the Hull-White model.
        /// This is a version done to double check the HW swaption.
        /// </summary>
        /// <param name='a'>
        /// Hull-White alpha parameter.
        /// </param>
        /// <param name='sigma'>
        /// Hull-White sigma parameter.
        /// </param>
        /// <param name='r'>
        /// Short rate value.
        /// </param>
        /// <param name='T'>
        /// Vector of bond maturity.
        /// </param>
        /// <param name='t'>
        /// Time at which the valuation is done.
        /// </param>
        /// <returns>
        /// The vector of bond prices.
        /// </returns>
        private Vector HWBond(double a, double sigma, double r, Vector T, double t)
        {
            double y = r - this.alphaTFunc(t);
            Vector result = new Vector(T.Length);
            for (int i = 0; i < T.Length; i++)
            {
                result[i] = Math.Exp(A(t, T[i], a, sigma, this.zeroRateCurve) - y * B(T[i] - t, a));
            }
            return result;
        }

        /// <summary>
        /// Calculates the value of a put option on a zero coupon bond.
        /// </summary>
        /// <param name='a'>
        /// Hull-White alpha parameter.
        /// </param>
        /// <param name='sigma'>
        /// Hull-White sigma parameter.
        /// </param>
        /// <param name='L'>
        /// Zero coupon notional.
        /// </param>
        /// <param name='K'>
        /// Bond rate.
        /// </param>
        /// <param name='T'>
        /// The maturity of the option.
        /// </param>
        /// <param name='s'>
        /// The maturity of the bond.
        /// </param>
        /// <returns>
        /// The put price.
        /// </returns>
        public double ZCBPut(double a, double sigma, double L, double K, double T, double s)
        {
            double h = this.H(a, sigma, L, K, T, s);
            double sigmap = this.SigmaP(a, sigma, T, s);
            return K * this.PZC(T) * SpecialFunctions.NormCdf(-h + sigmap) - L * this.PZC(s) * SpecialFunctions.NormCdf(-h);
        }

        /// <summary>
        /// Calculates H() function to be used in ZCBPut() method.
        /// </summary>
        /// <param name='a'>
        /// Hull-White alpha parameter.
        /// </param>
        /// <param name='sigma'>
        /// Hull-White sigma parameter.
        /// </param>
        /// <param name='L'>
        /// Zero coupon notional.
        /// </param>
        /// <param name='K'>
        /// Bond rate.
        /// </param>
        /// <param name='T'>
        /// The maturity of the option.
        /// </param>
        /// <param name='s'>
        /// The maturity of the bond.
        /// </param>
        /// <returns>
        /// A double with the value of the H() function.
        /// </returns>
        private double H(double a, double sigma, double L, double K, double T, double s)
        {
            double sigmap = this.SigmaP(a, sigma, T, s);
            return Math.Log(L * this.PZC(s) / (K * this.PZC(T))) / sigmap + 0.5 * sigmap;
        }

        /// <summary>
        /// Calculates SigmaP() function to be used in ZCBPut() method.
        /// </summary>
        /// <param name='a'>
        /// Hull-White alpha parameter.
        /// </param>
        /// <param name='sigma'>
        /// Hull-White sigma parameter.
        /// </param>
        /// <param name='T'>
        /// The maturity of the option.
        /// </param>
        /// <param name='s'>
        /// The maturity of the bond.
        /// </param>
        /// <returns>
        /// A double with the value of the SigmaP() function.
        /// </returns>
        private double SigmaP(double a, double sigma, double T, double s)
        {
            return sigma * (1.0 - Math.Exp(-a * (s - T))) * Math.Sqrt(0.5 * (1.0 - Math.Exp(-2.0 * a * T)) / a) / a;
        }

        /// <summary>
        /// Helper function to make functions easier to read.
        /// Just returns the value of the zero rate at position k.
        /// </summary>
        /// <param name="k">The position where to get the value of the zero rate from.</param>
        /// <returns>The value of the zero rate at position k.</returns>
        private double ZR(double k)
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
            return Math.Exp(-this.ZR(t) * t);
        }

        /// <summary>
        /// Numerically calculates the instantaneous forward rate.
        /// </summary>
        /// <param name='t'>
        /// Time at which calculate the forward rate.
        /// </param>
        /// <param name='dt'>
        /// Interval to be used in the numerical derivative.
        /// </param>
        /// <returns>
        /// The value of the instantaneous forward rate.
        /// </returns>
        private double F(double t, double dt)
        {
            if (t == 0)
                return this.ZR(t);
            else
                return (this.ZR(t + dt) * (t + dt) - this.ZR(t) * t) / dt;
        }

        /// <summary>
        /// Calculation of alpha function
        /// </summary>
        /// <param name="t">Time at which calculate the alpha function</param>
        /// <returns>Alpha function value</returns>
        public double alphaTFunc(double t)
        {
            double dt = 0.001;
            return this.F(t, dt) + this.sigma * this.sigma * Math.Pow(1.0 - Math.Exp(-this.a * t), 2.0) / (2.0 * this.a * this.a);
        }

        /// <summary>
        /// Calculates the function A() to be used in the Bond() method.
        /// </summary>
        /// <param name='t'>
        /// The time at which the Bond price will be calculated.
        /// </param>
        /// <param name='T'>
        /// The bond maturity.
        /// </param>
        /// <param name='alpha'>
        /// Hull-White alpha parameter.
        /// </param>
        /// <param name='sigma'>
        /// Hull-White sigma parameter.
        /// </param>
        /// <param name='zeroRateCurve'>
        /// Zero rate curve.
        /// </param>
        /// <returns>
        /// A double with the value of the A() function.
        /// </returns>
        public static double A(double t, double T, double alpha, double sigma, Function zeroRateCurve)
        {
            double dT = T - t;
            double firstTerm = sigma * sigma * (alpha * dT - 2.0 * (1.0 - Math.Exp(-alpha * dT))
                + 0.5 * (1.0 - Math.Exp(-2.0 * alpha * dT))) / (2.0 * Math.Pow(alpha, 3.0));

            return firstTerm - AlphaInt(t, T, alpha, sigma, zeroRateCurve);
        }

        /// <summary>
        /// Calculates the integral of alpha function to be used in the A() method.
        /// </summary>
        /// <param name="t">Lower value defining integration interval.</param>
        /// <param name="T">Upper value defining integration interval.</param>
        /// <param name="alpha">Hull-White alpha parameter.</param>
        /// <param name="sigma">Hull-White sigma parameter.</param>
        /// <param name="zeroRateCurve">Zero rate curve.</param>
        /// <returns>The integral of alpha function between t and T.</returns>
        private static double AlphaInt(double t, double T, double alpha, double sigma, Function zeroRateCurve)
        {
            double firstTerm = zeroRateCurve.Evaluate(T) * T - zeroRateCurve.Evaluate(t) * t;
            return firstTerm + sigma * sigma * (alpha * (T - t) - 2.0 * (Math.Exp(-alpha * t) - Math.Exp(-alpha * T))
                + 0.5 * (Math.Exp(-2.0 * alpha * t) - Math.Exp(-2.0 * alpha * T))) / (2.0 * Math.Pow(alpha, 3.0));
        }

        /// <summary>
        /// Calculates the function B() to be used in the Bond() method.
        /// </summary>
        /// <param name='T'>
        /// The difference between bond maturity time and valuation time.
        /// </param>
        /// <param name='alpha'>
        /// Hull-White alpha parameter.
        /// </param>
        /// <returns>
        /// The value of the B function.
        /// </returns>
        private static double B(double T, double alpha)
        {
            return (1.0 - Math.Exp(-alpha * T)) / alpha;
        }
    }
}
