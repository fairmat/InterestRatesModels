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
using DVPLDOM;
using DVPLI;
using Fairmat.Finance;
using Fairmat.Optimization;

namespace CIRProcess
{
    /// <summary>
    /// Optimization problem for CIR Cap Calibration.
    /// </summary>
    public class CapCIROptimizationProblem : IOptimizationProblem
    {
        /// <summary>
        /// The starting value of the process usually the zero rate evaluated in 0.
        /// </summary>
        public double r0;

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
        /// Cap tenor.
        /// </summary>
        private double tau;

        /// <summary>
        /// Constructor for the Cox-Ingersoll-Ross Calibration problem based on caps matrices.
        /// </summary>
        /// <param name="blackCaps">A Blacks caps matrix to use.</param>
        /// <param name="capMaturity">Vector of cap maturities.</param>
        /// <param name="capRate">Vector of cap strikes.</param>
        /// <param name="tau">Cap tenor.</param>
        /// <param name="r0">
        /// The starting value of the process usually the zero rate evaluated in 0.
        /// </param>
        public CapCIROptimizationProblem(Matrix blackCaps, Vector capMaturity, Vector capRate, double tau, double r0)
        {
            this.blackCaps = blackCaps;
            this.capMaturity = capMaturity;
            this.capRate = capRate;
            this.tau = tau;
            this.r0 = r0;
        }

        /// <summary>
        /// Constructor for the Cox-Ingersoll-Ross Calibration Problem based on caps matrices,
        /// using an <see cref="InterestRateMarketData"/> to derive the required data.
        /// </summary>
        /// <param name="irmd">
        /// An <see cref="InterestRateMarketData"/> containing the
        /// required information for the optimization problem.
        /// </param>
        public CapCIROptimizationProblem(InterestRateMarketData irmd)
        {
            this.capMaturity = irmd.CapMaturity;
            this.capRate = irmd.CapRate;
            this.tau = irmd.CapTenor;

            PFunction zr = new PFunction(null);
            zr.m_Function.iType = DVPLUtils.EInterpolationType.LINEAR;
            double[,] zrval = (double[,])ArrayHelper.Concat(irmd.ZRMarketDates.ToArray(),
                                                            irmd.ZRMarket.ToArray());
            zr.Expr = zrval;

            this.r0 = zr.Evaluate(0.0);

            BlackModel bm = new BlackModel(zr);
            this.blackCaps = new Matrix(this.capMaturity.Length, this.capRate.Length);
            for (int i = 0; i < this.capMaturity.Length; i++)
            {
                for (int j = 0; j < this.capRate.Length; j++)
                {
                    if (irmd.CapVolatility[i, j] == 0)
                        this.blackCaps[i, j] = 0;
                    else
                        this.blackCaps[i, j] = bm.Cap(this.capRate[j], irmd.CapVolatility[i, j],
                                                      this.tau, this.capMaturity[i]);

                    if (double.IsNaN(this.blackCaps[i, j]))
                        throw new Exception("Error on cap market price calculation");
                }
            }
        }

        #region IOptimizationProblem Members

        /// <summary>
        /// Calibration objective function:
        /// squared difference between black caps and Cox-Ingersoll-Ross caps.
        /// </summary>
        /// <param name='x'>
        /// The tested [kappa, theta, sigma] vector.
        /// </param>
        /// <returns>
        /// The distance (using the L2 norm) between
        /// black Caps and the Cox-Ingersoll-Ross Caps.
        /// </returns>
        public double Obj(DVPLI.Vector x)
        {
            Matrix cirCapMatrix = CIRCap.CIRCapMatrix(this.capMaturity, this.capRate,
                                                      this.tau, this.r0, x);
            double sum = 0;
            for (int r = 0; r < cirCapMatrix.R; r++)
            {
                for (int c = 0; c < cirCapMatrix.C; c++)
                {
                    if (this.blackCaps[r, c] != 0.0)
                        sum += Math.Pow(cirCapMatrix[r, c] - this.blackCaps[r, c], 2);
                }
            }

            return Math.Sqrt(sum);
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
        /// Gets the bounds for the optimization.
        /// </summary>
        public Bounds Bounds
        {
            get
            {
                Bounds b = new Bounds();
                b.Lb = (Vector)new double[] { 1e-5, 1e-5, 1e-5 };
                b.Ub = (Vector)new double[] { 5, 0.5, 0.5 };
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
