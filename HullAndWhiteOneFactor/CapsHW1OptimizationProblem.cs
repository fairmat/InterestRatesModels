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
using Fairmat.Optimization;

namespace HullAndWhiteOneFactor
{
    /*
    [SettingsContainer(@"HW")]
    [Mono.Addins.Extension("/Fairmat/UserSettings")]
    [Serializable]
    public class HWPreferences: ISettings
    {
        [SettingDescription("Alpha Lower Bound")]
        public double AlphaLB = .1;
    }
    */


    /// <summary>
    /// Implementation of HW1 Calibration Problem (Caps matrix based).
    /// </summary>
    public class CapsHW1OptimizationProblem : IOptimizationProblem
    {
        /// <summary>
        /// A <see cref="HullAndWhiteOneFactor.CapHW1"/> object used for calibration.
        /// </summary>
        private CapHW1 hw1Caps;

        /// <summary>
        /// A Blacks caps matrix used for calibration.
        /// </summary>
        private Matrix blackCaps;

        /// <summary>
        /// Vector of cap maturities.
        /// </summary>
        private Vector capMaturity;

        /// <summary>
        /// Vector of cap strikes.
        /// </summary>
        private Vector capRate;

        /// <summary>
        /// Time interval between caplets expressed in year fraction.
        /// </summary>
        private double deltaK;

        /// <summary>
        /// Constructor for the HW1 Calibration Problem based on caps matrices.
        /// </summary>
        /// <param name="hw1Caps">
        /// A <see cref="HullAndWhiteOneFactor.CapHW1"/> object to use.
        /// </param>
        /// <param name="blackCaps">A Blacks caps matrix to use.</param>
        /// <param name="capMaturity">Vector of cap maturities.</param>
        /// <param name="capRate">Vector of cap strikes.</param>
        /// <param name="deltaK">Interval between caplets expressed in year.</param>
        internal CapsHW1OptimizationProblem(CapHW1 hw1Caps, Matrix blackCaps, Vector capMaturity, Vector capRate, double deltaK)
        {
            this.hw1Caps = hw1Caps;
            this.blackCaps = blackCaps;
            this.capMaturity = capMaturity;
            this.capRate = capRate;
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
                //b.Lb = (Vector)new double[] {(UserSettings.GetSettings(typeof(HWPreferences)) as HWPreferences).AlphaLB  , 1e-8 };

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
        /// squared difference between black caps and hw caps.
        /// </summary>
        /// <param name='x'>
        /// The tested [alpha, sigma] vector.
        /// </param>
        /// <returns>
        /// A double containing the squared difference between
        /// black Caps and the HW Caps.
        /// </returns>
        public double Obj(DVPLI.Vector x)
        {
            Matrix hwCapsMatrix = this.hw1Caps.HWMatrixCaps(this.capMaturity, this.capRate, x[0], x[1], this.deltaK);
            double sum = 0;
            for (int r = 0; r < hwCapsMatrix.R; r++)
            {
                for (int c = 0; c < hwCapsMatrix.C; c++)
                {
                    if (this.blackCaps[r, c] != 0.0)
                        sum += Math.Pow(hwCapsMatrix[r, c] - this.blackCaps[r, c], 2);
                }
            }
            double bias = 0;
			/*
            f.zr = hw1Caps.zeroRateCurve;
            f.a = x[0];
            f.sigma = x[1];
            List<double> simDates, fR, avgR;
            f.Simulate(40, out simDates, out fR, out avgR);
            //double bias=f.ExpectedAverageRate(50) - hw1Caps.zeroRateCurve.Evaluate(50);

            for (int z = 0; z < simDates.Count; z++)
                bias += Math.Abs(avgR[z] - hw1Caps.zeroRateCurve.Evaluate(simDates[z]));
			Console.WriteLine(bias);
			*/
            double k = 25000;
            return Math.Sqrt(sum ) + k * bias;
			
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
