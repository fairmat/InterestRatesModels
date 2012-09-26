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
                    FSR = ForwardSwapRate(swaptionMaturity[i], swapPayDate);
                    result[i, j] = HWSwaption(a, sigma, 1000.0, FSR, swaptionMaturity[i], swapPayDate);
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
            double den = (T[0] - t) * PZC(T[0]);
            for (int i = 0; i < (T.Length - 1); i++)
                den += (T[i + 1] - T[i]) * PZC(T[i + 1]);
            return (PZC(t) - PZC(T[T.Length - 1])) / den;
        }

        /// <summary>
        /// Price of a swaption in HW1 model: the holder has the right to pay
        /// the fixed rate and receive floating rate.
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
            options.epsilon = 1e-5;

            this.dt = s[0] - T;
            this.CF = k * l * this.dt + (new Vector(s.Length));
            this.CF[s.Length - 1] = this.CF[s.Length - 1] + l;
            this.a = a;
            this.sigma = sigma;
            this.t = T;
            this.T = s;
            this.L = l;
            SolutionInfo sol = Fairmat.Optimization.Helper.FSolve(new ObjFunction(Func), (Vector)new double[1] { 0.01 },
                (Vector)new double[1] { -1.0 }, (Vector)new double[1] { 1.0 }, options);
            double RK = (double)sol.x;
            Vector PK = new Vector(s.Length);
            Vector vec = HWBond(a, sigma, RK, this.T, this.t, this.dt);
            for (int i = 0; i < this.CF.Length; i++)
                PK[i] = this.CF[i] * vec[i];
            double result = 0;
            for (int i = 0; i < s.Length; i++)
                result += ZCBPut(a, sigma, this.CF[i], PK[i], T, s[i]);
            return result;
        }

        public double HWSwaption2(double a, double sigma, double l, double k, double T, Vector s)
        {
            OptimizationSettings options = new OptimizationSettings();
            //options.epsilon = 1e-5;
            options.epsilon = 1e-10;

            this.dt = s[0] - T;
            this.CF = k * this.dt + (new Vector(s.Length));
            this.CF[s.Length - 1] = this.CF[s.Length - 1] + 1;
            this.a = a;
            this.sigma = sigma;
            this.t = T;
            this.T = s;
            this.L = l;
            SolutionInfo sol = Fairmat.Optimization.Helper.FSolve(new ObjFunction(Func2), (Vector)new double[1] { 0.01 },
                (Vector)new double[1] { -1.0 }, (Vector)new double[1] { 1.0 }, options);
            double RK = (double)sol.x;

            Vector X = HWBond2(a, sigma, RK, this.T, this.t, 0.001);
            double result = 0;
            for (int i = 0; i < s.Length; i++)
                result += CF[i]*ZCBPut(a, sigma, 1.0, X[i], T, s[i]);
            result = result * l;
            return result;
        }

        /// <summary>
        /// Function to be used in the FSolve() problem of HWSwaption() method.
        /// </summary>
        /// <param name='x'>
        /// Vector of length 1 representing the variable on which the FSolve problem is performed.
        /// </param>
        /// <returns>
        /// The value of the function given the provided vector x.
        /// </returns>
        public double Func(Vector x)
        {
            double result = this.CF.Scalar(HWBond(this.a, this.sigma, x[0], this.T, this.t, this.dt));
            return result - this.L;
        }

        public double Func2(Vector x)
        {
            double result = this.CF.Scalar(HWBond2(this.a, this.sigma, x[0], this.T, this.t, 0.001));
            return result - 1.0;
        }

        /// <summary>
        /// Calculates a vector of zero coupon bond prices
        /// (discount factors) within the Hull-White model.
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
        /// <param name='dt'>
        /// Time interval between bond coupons.
        /// </param>
        /// <returns>
        /// The vector of bond prices.
        /// </returns>
        private Vector HWBond(double a, double sigma, double r, Vector T, double t, double dt)
        {
            double Pt = PZC(t);
            double Pdt = PZC(t + dt);
            double Bdt = (1.0 - Math.Exp(-a * dt)) / a;

            Vector PT;
            Vector BT;
            Vector Bh;
            Vector Ah;
            Vector result;
            PT = new Vector(T.Length);
            BT = new Vector(T.Length);
            Bh = new Vector(T.Length);
            Ah = new Vector(T.Length);
            result = new Vector(T.Length);
            for (int i = 0; i < T.Length; i++)
            {
                PT[i] = PZC(T[i]);
                BT[i] = (1.0 - Math.Exp(-a * (T[i] - t))) / a;
                Bh[i] = dt * BT[i] / Bdt;
                Ah[i] = Math.Exp(Math.Log(PT[i] / Pt) - Math.Log(Pdt / Pt) * BT[i] / Bdt - sigma * sigma * (1.0 - Math.Exp(-2 * a * t)) * BT[i] * (BT[i] - Bdt) / (4 * a));
                result[i] = Ah[i] * Math.Exp(-Bh[i] * r);
            }

            return result;
        }

        private Vector HWBond2(double a, double sigma, double r, Vector T, double t, double dt)
        {
            double Pt = PZC(t);
            Vector PT;
            Vector BT;
            Vector AT;
            Vector result;
            PT = new Vector(T.Length);
            BT = new Vector(T.Length);
            double f = -( Math.Log(PZC(t + dt)) - Math.Log(PZC(t)) ) / ( dt );
            AT = new Vector(T.Length);
            result = new Vector(T.Length);
            for (int i = 0; i < T.Length; i++)
            {
                PT[i] = PZC(T[i]);
                BT[i] = (1.0 - Math.Exp(-a * (T[i] - t))) / a;
                AT[i] = (PT[i] / Pt)*Math.Exp(BT[i]*f - sigma * sigma * (1.0 - Math.Exp(-2 * a * t)) * BT[i] * BT[i] / (4 * a));
                result[i] = AT[i] * Math.Exp(-BT[i] * r);
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
        private double ZCBPut(double a, double sigma, double L, double K, double T, double s)
        {
            double h = H(a, sigma, L, K, T, s);
            double sigmap = SigmaP(a, sigma, T, s);
            return K * PZC(T) * Fairmat.Statistics.SpecialFunctions.NormCdf(-h + sigmap) - L * PZC(s) * Fairmat.Statistics.SpecialFunctions.NormCdf(-h);
        }

        private double ZCBPut2(double a, double sigma, double L, double K, double T, double s)
        {
            double h = H(a, sigma, L, K, T, s);
            double sigmap = SigmaP(a, sigma, T, s);
            return K * PZC(T) * Fairmat.Statistics.SpecialFunctions.NormCdf(-h + sigmap) - L * PZC(s) * Fairmat.Statistics.SpecialFunctions.NormCdf(-h);
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
            double sigmap = SigmaP(a, sigma, T, s);
            return Math.Log(L * PZC(s) / (K * PZC(T))) / sigmap + 0.5 * sigmap;
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
            return Math.Exp(-ZR(t) * t);
        }
    }
}
