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
using DVPLDOM;
using DVPLI;
using Fairmat.Optimization;
using Pelsser;

namespace Pelsser.Calibration
{
    /// <summary>
    /// Describes the optimization problem.
    /// </summary>
    public class PelsserCappletOptimizationProblem : IOptimizationProblem
    {
        /// <summary>
        /// A reference to the project where this optimization problem is solved.
        /// </summary>
        private ProjectROV project;

        /// <summary>
        /// A Pelsser caplet used to resolve the optimization problem.
        /// </summary>
        private Caplet caplet;

        /// <summary>
        /// Vector of cap maturities.
        /// </summary>
        private Vector capMaturity;

        /// <summary>
        /// The forward with the deltaK step.
        /// </summary>
        private Vector fwd;

        /// <summary>
        /// The caplet strike vector (columns).
        /// </summary>
        private Vector capK;

        /// <summary>
        /// The interval between caplets expressed in year.
        /// </summary>
        private double deltaK;

        /// <summary>
        /// The caplet maturities.
        /// </summary>
        private Vector capMat;

        /// <summary>
        /// A Blacks caps matrix used for calibration.
        /// </summary>
        private Matrix blackCaps;

        /// <summary>
        /// Constructor for the Pelsser Calibration Problem based on caps matrices.
        /// </summary>
        /// <param name="project">The project where this calibration is being done.</param>
        /// <param name="caplet">The Pelsser Caplet object to use.</param>
        /// <param name="capMaturity">Vector of cap maturities.</param>
        /// <param name="fwd">Forward with the deltaK step.</param>
        /// <param name="capK">The caplet strike vector (columns).</param>
        /// <param name="deltaK">Interval between caplets expressed in year.</param>
        /// <param name="capMat">The caplet maturities.</param>
        /// <param name="blackCaps">A Blacks caps matrix to use.</param>
        internal PelsserCappletOptimizationProblem(ProjectROV project, Caplet caplet, Vector capMaturity, Vector fwd, Vector capK, double deltaK, Vector capMat, Matrix blackCaps)
        {
            this.project = project;
            this.caplet = caplet;
            this.capMaturity = capMaturity;
            this.fwd = fwd;
            this.capK = capK;
            this.deltaK = deltaK;
            this.capMat = capMat;
            this.blackCaps = blackCaps;
            Console.WriteLine("Pelsser Caps problem on " + this.capMaturity.Length + " x " + this.capK.Length + " elements");
        }

        #region IOptimizationProblem Members

        /// <summary>
        /// Gets the definitions of the bounds where to seek for alpha and sigma.
        /// </summary>
        public Bounds Bounds
        {
            get
            {
                // Alpha must be positive, given the term alphaT[i] = F(current, dt) + 2 * Math.Exp(-alpha1Temp * current) * INT;
                Bounds bounds = new Bounds();
                bounds.Lb = new DVPLI.Vector() { 1e-8, 1e-6 };
                bounds.Ub = new DVPLI.Vector() { 1, .1 };
            
                return bounds;
            }
        }

        /// <summary>
        /// Gets a value indicating whether there are non linear constraints in this
        /// optimization problem. In this case there are.
        /// </summary>
        public bool HasNonLinearConstraints
        {
            get
            {
                return true;
            }
        }

        /// <summary>
        /// Gets null as we have no linear constrains defined.
        /// </summary>
        public LinearConstraints LinearIneqConstraints
        {
            get
            {
                return null;
            }
        }

        /// <summary>
        /// Creates a representation of the Pelsser constraints.
        /// </summary>
        /// <remarks>
        /// The method is static in order to be used by other classes.
        /// </remarks>
        /// <param name="project">The project where to evaluate the constraints.</param>
        /// <param name="x">The vector of x values.</param>
        /// <param name="maturity">The cap maturity to use to evaluate the constraints.</param>
        /// <returns>
        /// A vector with the representation of the Pelsser constraints.
        /// </returns>
        public static DVPLI.Vector PelsserConstraint(ProjectROV project, DVPLI.Vector x, double maturity)
        {
            // This is modeled as a one-dimensional bound.
            SquaredGaussianModel model = Assign(project, x);
            Vector res = new Vector(1);

            double dt = 1.0 / (2 * 252);
            for (double t = 0; t <= maturity; t += dt)
                res[0] += Math.Max(0, -model.F2(t, dt)); // The square of F.

            
            dt = .1317;
            for (double t = 0; t <= 2; t += dt)
                res[0] +=Math.Max(0,-model.F2(t, dt)); // The square of F.
          

            dt = .25;
            for (double t = 0; t <= 2; t += dt)
                res[0] +=Math.Max(0,-model.F2(t, dt)); // The square of F.
            
            
            return res;
        }

        /// <summary>
        /// Implements the G() function of the interface which describes
        /// a non linear inequality constraint.
        /// </summary>
        /// <param name="x">The vector of x values.</param>
        /// <returns>
        /// A vector with the representation of the Pelsser constraints.
        /// </returns>
        public DVPLI.Vector G(DVPLI.Vector x)
        {
            return PelsserConstraint(this.project, x, 50);
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
        /// Assigns new trial values for x in order to use it
        /// for the evaluation of the constraints and the objective function.
        /// </summary>
        /// <param name="project">The project where to evaluate.</param>
        /// <param name="x">The vector of x values.</param>
        /// <returns>A <see cref="SquaredGaussianModel"/> set with the values.</returns>
        internal static SquaredGaussianModel Assign(ProjectROV project, Vector x)
        {
            Engine.Parser.NewContext();

            SquaredGaussianModel model = project.Processes[0].Plugin as SquaredGaussianModel;

            // Assign the new values.
            ModelParameter a1 = model.a1 as ModelParameter;
            ModelParameter s1 = model.sigma1 as ModelParameter;
            a1.m_Value = (RightValue)x[0];
            s1.m_Value = (RightValue)x[1];

            bool errors=model.Parse(null);
            if (errors) 
                 throw new Exception("Cannot Assing model");
           
            return model;
        }

        /// <summary>
        /// Calibration objective function:
        /// squared difference between black caps and Pelsser caps.
        /// </summary>
        /// <param name='x'>
        /// The tested [alpha, sigma] vector.
        /// </param>
        /// <returns>
        /// A scalar containing the sum of squared differences between
        /// black Caps and the Pelsser Caps.
        /// </returns>
        public double Obj(DVPLI.Vector x)
        {
            SquaredGaussianModel model = Assign(this.project, x);

            try
            {
                Matrix caps = this.caplet.PGSMCaplets(model, this.capMaturity, this.fwd, this.capK, this.deltaK, this.capMat);
                double sum = 0;
                for (int r = 0; r < caps.R; r++)
                {
                    for (int c = 0; c < caps.C; c++)
                    {
                        if (this.blackCaps[r, c] != 0.0)
                            sum += Math.Pow(caps[r, c] - this.blackCaps[r, c], 2);
                    }
                }

                return Math.Sqrt(sum/(caps.R*caps.C));
            }
            catch (Exception)
            {
                return 1000 * x.Norm(NormType.L2);
            }
        }

        #endregion
    }
}
