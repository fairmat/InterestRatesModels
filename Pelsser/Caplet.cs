/* Copyright (C) 2009-2012 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
 * Author(s): Matteo Tesser (matteo.tesser@fairmat.com)
 *            Michele Furgeri (info@fairmat.com)
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
using System.Threading.Tasks;
using DVPLI;
using DVPLI.TaskScheduler;
using DVPLUtils;
using Fairmat.Statistics;
using Pelsser;
using Fairmat.Math;

namespace Pelsser.Calibration
{
    /// <summary>
    /// Caplet prices using Pelsser Squared Gaussian model.
    /// </summary>
    public class Caplet
    {
        /// <summary>
        /// Keeps the data of a context for this caplet calibrator
        /// and provides helpers to handle these data fields.
        /// </summary>
        private class Context
        {
            /// <summary>
            /// Gets or sets a cache for the pre-calculated values of
            /// s^2*C(0,Ti) for Ti taken from a vector T.
            /// </summary>
            internal Vector Sigma0s { get; set; }

            /// <summary>
            /// Gets or sets a cache for the value of C(deltaK).
            /// </summary>
            internal double CCost { get; set; }

            /// <summary>
            /// Gets or sets the <see cref="SquaredGaussianModel"/> used in this context.
            /// </summary>
            internal SquaredGaussianModel Model { get; set; }

            /// <summary>
            /// Gets or sets a Vector of cap maturities.
            /// </summary>
            internal Vector Mat { get; set; }

            /// <summary>
            /// Gets or sets the forward with the deltaK step.
            /// </summary>
            internal Vector Fwd { get; set; }

            /// <summary>
            /// Gets or sets a vector of strikes.
            /// </summary>
            internal Vector RK { get; set; }

            /// <summary>
            /// Gets or sets the step to use for each iteration.
            /// </summary>
            internal double DeltaK { get; set; }

            /// <summary>
            /// Gets or sets the rows of the volatility surface.
            /// </summary>
            internal Vector TOss { get; set; }

            /// <summary>
            /// Gets or sets the caplets <see cref="Matrix"/>.
            /// </summary>
            internal Matrix Caplets { get; set; }

            /// <summary>
            /// Gets or sets a <see cref="Vector"/> of prices.
            /// </summary>
            internal Vector K { get; set; }

            /// <summary>
            /// Gets or sets a <see cref="Vector"/> with the logarithm of the prices.
            /// </summary>
            internal Vector LogK { get; set; }

            /// <summary>
            /// Gets or sets the starting row, used by the parallel version.
            /// </summary>
            internal int RowStart { get; set; }

            /// <summary>
            /// Gets or sets the ending row, used by the parallel version.
            /// </summary>
            internal int RowEnd { get; set; }

            /// <summary>
            /// Callback for calculating a series of caplets matrix rows.
            /// </summary>
            /// <param name="c">The context where the rows to be calculated are stored.</param>
            internal static void CalculateRowP(object c)
            {
                Context context = c as Context;
                for (int r = context.RowStart; r <= context.RowEnd; r++)
                    context.CalculateRow(r);
            }

            /// <summary>
            /// Calculates an entire caplets matrix row using Pelsser's model.
            /// </summary>
            /// <param name="r">The row index.</param>
            internal void CalculateRow(int r)
            {
                double T = this.Mat[r + 1] + this.DeltaK;
                PelsserCache pc = new PelsserCache(this.Mat[r + 1], T, this.Model);

                double MU0 = this.Model.Mu0(this.Mat[r + 1]);

                int rm = r + 1;
                double nu = MU0 - pc.B * this.Sigma0s[rm];

                double phi = 1 + 2 * this.CCost * this.Sigma0s[rm];
                PelsserCache pc0t = new PelsserCache(0, this.Mat[r + 1], this.Model);
                PelsserCache pc0tau = new PelsserCache(0, T, this.Model);

                // Discount factor
                double p0tau = Math.Exp(pc0tau.A);
                double p0t = Math.Exp(pc0t.A);
                double sqrtSigmaRm = Math.Sqrt(this.Sigma0s[rm]);

                // Execute the subsequent operations on all caplets.
                for (int c = 0; c < this.K.Length; c++)
                {
                    double d = Math.Pow(pc.B, 2) + 4 * this.CCost * (pc.A - this.LogK[c]);

                    // Exclude the negative values of the discriminant.
                    d = Math.Max(0, d);
                    double el = (-pc.B - Math.Sqrt(d)) / (2 * this.CCost);
                    double eic = (-pc.B + Math.Sqrt(d)) / (2 * this.CCost);

                    double nP1 = SpecialFunctions.NormCdf(-(eic * phi - nu) / Math.Sqrt(phi * this.Sigma0s[rm])) + SpecialFunctions.NormCdf((el * phi - nu) / Math.Sqrt(phi * this.Sigma0s[rm]));
                    double nP2 = SpecialFunctions.NormCdf(-(eic - MU0) / sqrtSigmaRm) + SpecialFunctions.NormCdf((el - MU0) / sqrtSigmaRm);

                    double put_zc = -p0tau * nP1 + p0t * this.K[c] * nP2;

                    this.Caplets[r, c] = put_zc;
                }
            }
        }

        /// <summary>
        /// Calculates the integral C(t, T) function with T taking
        /// every element of the parameter mat.
        /// </summary>
        /// <param name="model">The model instance.</param>
        /// <param name="t">The starting point of the integral.</param>
        /// <param name="mat">The vector of ending interval points T.</param>
        /// <returns>The vector of pre-calculated C(t,T).</returns>
        private Vector CtT(SquaredGaussianModel model, double t, Vector mat)
        {
            Vector C = new Vector(mat.Length);
            for (int i = 0; i < mat.Length; i++)
                C[i] = model.C(mat[i] - t);
            return C;
        }

        /// <summary>
        /// Caplet prices calculated as a put on a zero coupon bond.
        /// </summary>
        /// <param name="model">The model to use to execute the calculation.</param>
        /// <param name="mat">
        /// Caplet maturity. This vector starts from zero and increases
        /// of step DeltaK each element till the last one.
        /// </param>
        /// <param name="fwd">Forward with the deltaK step.</param>
        /// <param name="rk">Strike vector (columns).</param>
        /// <param name="deltaK">Amount to use as increase factor.</param>
        /// <param name="tOss">The Maturities.</param>
        /// <returns>A <see cref="Matrix"/> with the caplet prices.</returns>
        public Matrix PGSMCaplets(SquaredGaussianModel model, Vector mat, Vector fwd, Vector rk, double deltaK, Vector tOss)
        {
            double s = model.sigma1.fV();
            int col = rk.Length;

            int NP = (int)(1 + tOss[tOss.Length - 1] * 252);
            double[] dates = new double[NP];
            double step = mat[mat.Length - 1] / (NP - 1);
            for (int z = 0; z < NP; z++)
                dates[z] = step * z;

            DateTime t0 = DateTime.Now;
            model.Setup(dates);
            DateTime t1 = DateTime.Now;

            Vector K = 1.0 / (1 + rk * deltaK);
            double cCost = model.C(deltaK);
            Vector sigma0s = Math.Pow(s, 2) * CtT(model, 0, mat);

            Matrix caplets = new Matrix(mat.Length - 1, rk.Length);
            Matrix caps = new Matrix(tOss.Length, rk.Length);

            // Pre-calculate values.
            Vector logK = Vector.Log(K);

            // Parallel version.
            List<Task> tl = new List<Task>();

            Context context = new Context();
            context.Model = model;
            context.Mat = mat;
            context.Fwd = fwd;
            context.RK = rk;
            context.DeltaK = deltaK;
            context.TOss = tOss;
            context.K = K;
            context.LogK = logK;
            context.Caplets = caplets;
            context.Sigma0s = sigma0s;
            context.CCost = cCost;
            context.RowStart = 0;
            context.RowEnd = (mat.Length - 2) / 2;
            tl.Add(Task.Factory.StartNew(Context.CalculateRowP, context));

            context = new Context();
            context.Model = model;
            context.Mat = mat;
            context.Fwd = fwd;
            context.RK = rk;
            context.DeltaK = deltaK;
            context.TOss = tOss;
            context.K = K;
            context.LogK = logK;
            context.Caplets = caplets;
            context.Sigma0s = sigma0s;
            context.CCost = cCost;
            context.RowStart = (mat.Length - 2) / 2 + 1;
            context.RowEnd = mat.Length - 2 - 1;
            tl.Add(Task.Factory.StartNew(Context.CalculateRowP, context));

            tsScheduler.WaitTaskList(tl);

            // Sequential version.
            /*
            Context Context = new Context();
            Context.Model = Model;
            Context.Mat = Mat;
            Context.Fwd = Fwd;
            Context.RK = RK;
            Context.DeltaK = DeltaK;
            Context.TOss = TOss;
            Context.K = K;
            Context.LogK = LogK;
            Context.Caplets = Caplets;
            Context.Sigma0s = Sigma0s;
            Context.CCost = CCost;

           for (int r = 0; r < Mat.Length - 2; r++)
              Context.CalculateRow(r);
           */

            // Calculate the caps from the caplets.
            for (int r = 0; r < tOss.Length; r++)
            {
                for (int c = 0; c < rk.Length; c++)
                {
                        double current = 0;
                        for (int ci = 0; ci < (int)(tOss[r] / deltaK) - 1; ci++)
                            current += caplets[ci, c];
                        caps[r, c] = current;
                }
            }

            return caps;
        }
    }
}
