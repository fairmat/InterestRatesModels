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
        public double Obj(DVPLI.Vector x)
        {
            Matrix hwSWMatrix = this.shw1.HWSwaptionMatrix(this.swaptionMaturity, this.swapDuration, x[0], x[1], this.deltaK);
            double sum = 0;
            int count = this.swaptionMaturity.Length * this.swapDuration.Length;
            for (int r = 0; r < this.swaptionMaturity.Length; r++)
                for (int c = 0; c < this.swapDuration.Length; c++)
                    if (this.blackSwaption[r, c] != 0.0)
                        sum += Math.Pow(hwSWMatrix[r, c] - this.blackSwaption[r, c], 2);
            return Math.Sqrt(sum / count);
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
