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
    /// Implementation of HW1 Calibration Problem (Swaption matrix based).
    /// </summary>
    public class SwaptionHW1OptimizationProblem : IOptimizationProblem
    {
        /// <summary>
        /// A <see cref="HullAndWhiteOneFactor.SwaptionHW1"/> object used for calibration.
        /// </summary>
        private SwaptionHW1 shw1;

        /// <summary>
        /// A Blacks swaption matrix used for calibration.
        /// </summary>
        private Matrix blackSwaption;

        /// <summary>
        /// Vector of swaption maturities.
        /// </summary>
        private Vector swaptionMaturity;

        /// <summary>
        /// Vector of swaps duration.
        /// </summary>
        private Vector swapDuration;

        /// <summary>
        /// Time interval between swaptions expressed in year fraction.
        /// </summary>
        private double deltaK;

        /// <summary>
        /// Constructor for the HW1 Calibration Problem based on swaption matrices.
        /// </summary>
        /// <param name="shw1">
        /// A <see cref="HullAndWhiteOneFactor.SwaptionHW1"/> object to use.
        /// </param>
        /// <param name="blackSwaption">A Blacks swaption matrix to use.</param>
        /// <param name="swaptionMaturity">Vector of swaption maturities.</param>
        /// <param name="swapDuration">Vector of swaps duration.</param>
        /// <param name="deltaK">Interval between swaptions expressed in year.</param>
        internal SwaptionHW1OptimizationProblem(SwaptionHW1 shw1, Matrix blackSwaption, Vector swaptionMaturity, Vector swapDuration, double deltaK)
        {
            this.shw1 = shw1;
            this.blackSwaption = blackSwaption;
            this.swaptionMaturity = swaptionMaturity;
            this.swapDuration = swapDuration;
            this.deltaK = deltaK;
        }

        #region IOptimizationProblem Members

        /// <summary>
        /// Gets the bounds for the optimization.
        /// </summary>
        public Bounds Bounds
        {
            get
            {
                Bounds b = new Bounds();
                b.Lb = (Vector)new double[] { HW1.alphaLowerBound , 1e-8 };
                //b.Lb = (Vector)new double[] { (UserSettings.GetSettings(typeof(HWPreferences)) as HWPreferences).AlphaLB, 1e-8 };
                b.Ub = (Vector)new double[] { 1 - 1e-8, 0.5 };
                return b;
            }
        }

        /// <summary>
        /// Gets null as we have no linear constrains defined.
        /// </summary>
        public virtual LinearConstraints LinearIneqConstraints
        {
            get
            {
                return null;
            }
        }

        /// <summary>
        /// Gets a value indicating whether there are non linear constraints in this
        /// optimization problem. In this case there are not.
        /// </summary>
        public bool HasNonLinearConstraints
        {
            get
            {
                return false;
            }
        }
        HWCompactSimulator f = new HWCompactSimulator();// { a = x[0], sigma = x[1], zr = hw1Caps.zeroRateCurve };
        /// <summary>
        /// Calibration objective function:
        /// squared difference between black swaptions and HW swaptions.
        /// </summary>
        /// <param name='x'>
        /// The tested [alpha, sigma] vector.
        /// </param>
        /// <returns>
        /// A double containing the squared difference between
        /// black swaptions and the HW swaptions.
        /// </returns>
        public double Obj(Vector x)
        {
            return Obj(x, false);
        }


        internal double Obj(Vector x,bool showDetails)
        {
            var hwSWMatrix = this.shw1.HWSwaptionMatrix(this.swaptionMaturity, this.swapDuration, x[0], x[1], this.deltaK);
            double sum = 0;
            double sump = 0;
            int count = this.swaptionMaturity.Length * this.swapDuration.Length;
            int effective = 0; // number of effective elements
            for (int r = 0; r < this.swaptionMaturity.Length; r++)
                for (int c = 0; c < this.swapDuration.Length; c++)
                    if (this.blackSwaption[r, c] != 0.0)
                    {
                        double deltap = hwSWMatrix[r, c] - this.blackSwaption[r, c];
                        sum += Math.Pow(deltap, 2);
                        sump += Math.Pow(deltap / this.blackSwaption[r, c], 2);
                        effective++;
                    }


            double bias = 0;
            double k = 25000;
            /*
            f.zr = this.shw1.zeroRateCurve;
            f.a = x[0];
            f.sigma = x[1];
            List<double> simDates, fR, avgR;
            f.Simulate(40, out simDates, out fR, out avgR);
            bias += Math.Abs(avgR[avgR.Count - 1] - f.zr.Evaluate(simDates[avgR.Count - 1]));
            Console.WriteLine(f.zr.Evaluate(simDates[avgR.Count - 1]) + "\t" + avgR[avgR.Count - 1] + "\t" + x[0] + "\t" + x[1] + "\t" + Math.Sqrt(sum / count));
            */
            if (showDetails)
            {
                Console.WriteLine("Black Prices");
                Console.WriteLine(blackSwaption);
                Console.WriteLine("HW Prices");
                Console.WriteLine(hwSWMatrix);
                Console.WriteLine("Delta Prices");
                Console.WriteLine(hwSWMatrix - blackSwaption);
                double absErr = Matrix.Abs(hwSWMatrix - blackSwaption).Sum() / effective;
                double avgPrice= blackSwaption.Sum() / effective;
                Console.WriteLine("Abs Error\t"+ absErr+"\t("+ (absErr/ avgPrice).ToString("P") + ")");
                Console.WriteLine("RMSE(r)\t"+ (Math.Sqrt(sump)/effective).ToString("P");
            }
                



            return Math.Sqrt(sum / count )+k*bias;
        }

        /// <summary>
        /// This method is unused but part of the interface.
        /// </summary>
        /// <param name="x">The parameter is not used.</param>
        /// <returns>Nothing the function always throws a NotImplementedException.</returns>
        public DVPLI.Vector Grad(DVPLI.Vector x)
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// This method is unused but part of the interface.
        /// </summary>
        /// <param name="x">The parameter is not used.</param>
        /// <returns>Nothing the function always throws a NotImplementedException.</returns>
        public DVPLI.Vector G(DVPLI.Vector x)
        {
            throw new NotImplementedException();
        }

        #endregion
    }
}
