/* Copyright (C) 2009-2012 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
 * Author(s): Matteo Tesser (matteo.tesser@fairmat.com)
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
using System.Runtime.Serialization;
using DVPLDOM;
using DVPLI;
using Mono.Addins;

namespace CIRProcess
{
    /// <summary>
    /// Implementation of CIR model for pricing Credits Default Swaps.
    /// </summary>
    [Serializable]
    public unsafe class CIR2 : IExtensibleProcess, IPopulable, IMarkovSimulator, IParsable,
                               IPostSimulationTransformation, IDefaultProbabiliy
    {
        #region SerializedFields

        /// <summary>
        /// The mean reversion speed (first component).
        /// </summary>
        public IModelParameter k1;

        /// <summary>
        /// The mean reversion speed (second component).
        /// </summary>
        public IModelParameter k2;

        /// <summary>
        /// The long term mean (first component).
        /// </summary>
        public IModelParameter theta1;

        /// <summary>
        /// The long term mean (second component).
        /// </summary>
        public IModelParameter theta2;

        /// <summary>
        /// The volatility (first component).
        /// </summary>
        public IModelParameter sigma1;

        /// <summary>
        /// The volatility (second component).
        /// </summary>
        public IModelParameter sigma2;

        /// <summary>
        /// The starting value for the process (first component).
        /// </summary>
        public IModelParameter startingValue1;

        /// <summary>
        /// The starting value for the process (second component).
        /// </summary>
        public IModelParameter startingValue2;

        /// <summary>
        /// Reference to the zero rate.
        /// </summary>
        [ExternalSymbolReference("ZR", typeof(PFunction))]
        public IModelParameter ZRReference;

        #endregion SerializedFields

        /// <summary>
        /// Keeps a cached evaluated copy of the zero rate.
        /// </summary>
        [NonSerialized]
        private Function zr;

        /// <summary>
        /// The squared value of the parameter k (first component), used for caching reasons.
        /// </summary>
        [NonSerialized]
        private double k12;

        /// <summary>
        /// The squared value of the parameter k (second component), used for caching reasons.
        /// </summary>
        [NonSerialized]
        private double k22;

        /// <summary>
        /// The squared value of the parameter sigma (first component), used for caching reasons.
        /// </summary>
        [NonSerialized]
        private double sigma12;

        /// <summary>
        /// The squared value of the parameter sigma (second component), used for caching reasons.
        /// </summary>
        [NonSerialized]
        private double sigma22;

        /// <summary>
        /// Temporary value used to speedup computation. See  <see cref="Parse"/> method.
        /// </summary>
        [NonSerialized]
        private double gamma1;

        /// <summary>
        /// Temporary value used to speedup computation. See  <see cref="Parse"/> method.
        /// </summary>
        [NonSerialized]
        private double gamma2;

        /// <summary>
        /// Gets the process type name.
        /// </summary>
        public static string ProcessType
        {
            get
            {
                return "Interest Rate Models/CIR (Two Factors)";
            }
        }

        /// <summary>
        /// Ensure the parameters are correct.
        /// </summary>
        /// <param name='context'>
        /// The underlying project.
        /// </param>
        /// <returns>
        /// False if there were no parse errors.
        /// </returns>
        public bool Parse(IProject context)
        {
            bool errors = false;
            errors = BoolHelper.AddBool(errors, this.startingValue1.Parse(context));
            errors = BoolHelper.AddBool(errors, this.startingValue2.Parse(context));
            errors = BoolHelper.AddBool(errors, this.k1.Parse(context));
            errors = BoolHelper.AddBool(errors, this.k2.Parse(context));
            errors = BoolHelper.AddBool(errors, this.theta1.Parse(context));
            errors = BoolHelper.AddBool(errors, this.theta2.Parse(context));
            errors = BoolHelper.AddBool(errors, this.sigma1.Parse(context));
            errors = BoolHelper.AddBool(errors, this.sigma2.Parse(context));

            errors = BoolHelper.AddBool(errors, this.ZRReference.Parse(context));

            // Stores some temporary parameters derived from the main model parameters.
            if (!errors)
            {
                this.zr = (Function)this.ZRReference.fVRef();

                // In this way they are not stochastic.
                this.k12 = Math.Pow(this.k1.fV(), 2);
                this.k22 = Math.Pow(this.k2.fV(), 2);
                this.sigma12 = Math.Pow(this.sigma1.fV(), 2);
                this.sigma22 = Math.Pow(this.sigma2.fV(), 2);
                this.gamma1 = Math.Sqrt(this.k12 + 2 * this.sigma12);
                this.gamma2 = Math.Sqrt(this.k22 + 2 * this.sigma22);
            }

            return errors;
        }

        /// <summary>
        /// Default constructor which prepares an instance
        /// of the two factors Cox-Ingersoll-Ross model with some sample values.
        /// </summary>
        public CIR2()
        {
            DefaultValues();
        }

        /// <summary>
        /// Checks if the value is negative and floors it to zero, instead of attempting
        /// a square root of negative numbers.
        /// </summary>
        /// <param name="x">The value to do the square root of.</param>
        /// <returns>0 if x is negative, otherwise the square root of x.</returns>
        private double AdjSqrt(double x)
        {
            if (x < 0) return 0;
            else
                return Math.Sqrt(x);
        }

        #region IMarkovSimulator Members

        /// <summary>
        /// This function defines the drift in the Cox-Ingersoll-Ross Two factors Markov process.
        /// The formula to calculate the A component is
        /// A[0] = k1 * (theta1 - previous state)
        /// A[1] = k2 * (theta2 - previous state).
        /// </summary>
        /// <param name="i">The parameter is not used.</param>
        /// <param name="x">The state vector at the previous state.</param>
        /// <param name="a">The output of the function.</param>
        public void a(int i, double* x, double* a)
        {
            a[0] = this.k1.fV() * (this.theta1.fV() - x[0]);
            a[1] = this.k2.fV() * (this.theta2.fV() - x[1]);
        }

        /// <summary>
        /// This function defines the volatility in the
        /// Cox-Ingersoll-Ross Two factors Markov process.
        /// The formula to calculate the B component is
        /// B[0] = sqrt(previous state) * sigma1
        /// B[1] = sqrt(previous state) * sigma2.
        /// </summary>
        /// <param name="i">The parameter is not used.</param>
        /// <param name="x">The state vector at the previous state.</param>
        /// <param name="b">The output of the function.</param>
        public void b(int i, double* x, double* b)
        {
            b[0] = AdjSqrt(x[0]) * this.sigma1.fV();
            b[1] = AdjSqrt(x[1]) * this.sigma2.fV();
        }

        /// <summary>
        /// This function calculated drift and volatility in the
        /// Cox-Ingersoll-Ross Two factors Markov process.
        /// The formula to calculate the B component is
        /// B[0] = sqrt(previous state) * sigma1
        /// B[1] = sqrt(previous state) * sigma2.
        /// </summary>
        /// <param name="i">The parameter is not used.</param>
        /// <param name="x">The state vector at the previous state.</param>
        /// <param name="b">The output drift.</param>
        /// <param name="b">The output volatility.</param>
        public void ab(int i, double* x, double* a,double *b)
        {
            a[0] = this.k1.fV() * (this.theta1.fV() - x[0]);
            a[1] = this.k2.fV() * (this.theta2.fV() - x[1]);
            
            b[0] = AdjSqrt(x[0]) * this.sigma1.fV();
            b[1] = AdjSqrt(x[1]) * this.sigma2.fV();
        }



        /// <summary>
        /// Gets the starting point for the process, which is
        /// defined in <see cref="startingValue1"/> and <see cref="startingValue2"/>.
        /// </summary>
        public double[] x0
        {
            get
            {
                double[] r0 = new double[2];
                r0[0] = this.startingValue1.fV();
                r0[1] = this.startingValue2.fV();
                return r0;
            }
        }

        /// <summary>
        /// Sets the passed array with a Boolean stating if the process
        /// must be simulated as a log-normal process.
        /// </summary>
        /// <param name="isLog">
        /// A reference to the array to be set with the required information.
        /// The array contains a two members set to false as this simulation
        /// is not log-normal.
        /// </param>
        public void isLog(ref bool[] isLog)
        {
            isLog[0] = false;
            isLog[1] = false;
        }

        /// <summary>
        /// Gets details about the structure of the functions A and B of the Markov
        /// process.
        /// In this case drift and volatility are state dependant and not time dependant.
        /// </summary>
        public DynamicInfo DynamicInfo
        {
            get
            {
                return new DynamicInfo(false, true, false, true);
            }
        }

        #endregion IMarkovSimulator Members

        #region IExtensibleProcess Members

        /// <summary>
        /// Sets the descriptions for the components.
        /// </summary>
        private void SetDescription()
        {
            this.startingValue1.Description = "starting value factor 1";
            this.startingValue2.Description = "starting value factor 2";

            this.k1.Description = "k1";
            this.k2.Description = "k2";
            this.theta1.Description = "Theta1";
            this.theta2.Description = "Theta2";
            this.sigma1.Description = "Sigma1";
            this.sigma2.Description = "Sigma2";

            this.ZRReference.Description = "zero rate function reference";
        }

        /// <summary>
        /// Sets the default values for the components.
        /// </summary>
        public void DefaultValues()
        {
            this.k1 = new ModelParameter(0.026758131);
            this.k2 = new ModelParameter(0.226406137);
            this.theta1 = new ModelParameter(0.023147821);
            this.theta2 = new ModelParameter(0.01);
            this.sigma1 = new ModelParameter(0.1);
            this.sigma2 = new ModelParameter(0.0001);
            this.startingValue1 = new ModelParameter(0.00001);
            this.startingValue2 = new ModelParameter(0.00001);
            this.ZRReference = (ModelParameter)"@ZR";
            SetDescription();
        }

        /// <summary>
        /// Gets a value indicating whether FullSimulation is implemented, in this
        /// case it doesn't so it always returns false.
        /// </summary>
        public bool ImplementsFullSimulation
        {
            get
            {
                return false;
            }
        }

        /// <summary>
        /// Gets a value indicating whether a Markov based simulation is implemented, in this
        /// case it does so it always returns true.
        /// </summary>
        public bool ImplementsMarkovBasedSimulation
        {
            get
            {
                return true;
            }
        }

        /// <summary>
        /// Gets the ProcessInfo for this plugin, in this case
        /// Interest Rate Models/CIR (Two Factors).
        /// </summary>
        public ProcessInfo ProcessInfo
        {
            get
            {
                return new ProcessInfo(CIR2.ProcessType);
            }
        }

        /// <summary>
        /// Called by Simulator after parse.
        /// In this case does nothing.
        /// </summary>
        /// <param name='simulationDates'>The parameter is not used.</param>
        public void Setup(double[] simulationDates)
        {
            return;
        }

        /// <summary>
        /// Gets the information required in order to allow the simulation to run.
        /// CIR has:
        /// * 0 latent components in the state components.
        /// * 2 component of noise.
        /// * 2 state component.
        /// * The components names are s1 and s2.
        /// </summary>
        public SimulationInfo SimulationInfo
        {
            get
            {
                SimulationInfo simulationInfo = new SimulationInfo();
                simulationInfo.LatentSize = 0;
                simulationInfo.NoiseSize = 2;
                simulationInfo.StateSize = 2;
                simulationInfo.StateDescription = new string[] { "s1", "s2" };
                return simulationInfo;
            }
        }

        #endregion

        #region IExportableContainer Members

        /// <summary>
        /// Creates a list of all the sub-objects that can be edited.
        /// </summary>
        /// <param name="recursive">
        /// The parameter is not used.
        /// </param>
        /// <returns>
        /// The created list with all the sub objects that can be edited.
        /// </returns>
        public List<IExportable> ExportObjects(bool recursive)
        {
            List<IExportable> l = new List<IExportable>();
            l.Add(this.startingValue1);
            l.Add(this.startingValue2);

            l.Add(this.k1);
            l.Add(this.k2);
            l.Add(this.theta1);
            l.Add(this.theta2);
            l.Add(this.sigma1);
            l.Add(this.sigma2);

            l.Add(this.ZRReference);

            return l;
        }

        #endregion

        #region IPopulable Members

        /// <summary>
        /// This is part of the IPopulable interface but it's not implemented.
        /// </summary>
        /// <param name="names">The parameter is not used.</param>
        /// <param name="values">The parameter is not used.</param>
        public void Populate(string[] names, double[] values)
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Helper function to make functions easier to read.
        /// Just returns the value of the zero rate at position t.
        /// </summary>
        /// <param name="t">The position where to get the value of the zero rate from.</param>
        /// <returns>The value of the zero rate at position t.</returns>
        private double GetZR(double t)
        {
            return this.zr.Evaluate(t);
        }

        #endregion

        #region IPostSimulationTransformation Members
        /// <summary>
        /// Calculates the probabilities for one scenario.
        /// </summary>
        /// <param name="dates">Dates (in Year fractions).</param>
        /// <param name="outDynamic">Input and output matrix.</param>
        public void Transform(double[] dates, IMatrixSlice outDynamic)
        {
            for (int j = 0; j < dates.Length; j++)
            {
                // Truncate to zero the first two components.
                outDynamic[j, 0] = Math.Max(0.0000000001, outDynamic[j, 0]);
                outDynamic[j, 1] = Math.Max(0.0000000001, outDynamic[j, 1]);
            }
        }

        #endregion

        #region IDefaultProbabiliy Members

        /// <summary>
        /// Calculates the probability of a default at date t (conditional
        /// of being in a given state (realization) of  process s1, s2),
        /// being on period s = dates[j].
        /// </summary>
        /// <param name="dynamic">The underlying's realizations.</param>
        /// <param name="dates">The corresponding underlying's dates.</param>
        /// <param name="j">The wanted date index.</param>
        /// <param name="t">The date at which we want to estimate the default prob.</param>
        /// <returns>The default probability estimation.</returns>
        public double CDSDefaultP(IReadOnlyMatrixSlice dynamic, double[] dates, int j, double t)
        {
            double tau = t - dates[j];

            double RRate = 0.4;
            double s1 = dynamic[j, 0];
            double s2 = dynamic[j, 1];

            double final = t;
            double ZR = GetZR(final);

            // PTau.
            double discf = Math.Exp(-ZR * tau);

            double num_A1 = 2 * this.gamma1 * Math.Exp((this.k1.fV() + this.gamma1) * tau / 2.0);
            double num_A2 = 2 * this.gamma2 * Math.Exp((this.k2.fV() + this.gamma2) * tau / 2.0);

            double eg1 = Math.Exp(this.gamma1 * tau) - 1;
            double eg2 = Math.Exp(this.gamma2 * tau) - 1;

            double den_A1 = (this.k1.fV() + this.gamma1) * eg1 + 2 * this.gamma1;
            double den_A2 = (this.k2.fV() + this.gamma2) * eg2 + 2 * this.gamma2;
            double A1 = Math.Pow(num_A1 / den_A1, (2 * this.k1.fV() * this.theta1.fV()) / this.sigma12);
            double A2 = Math.Pow(num_A2 / den_A2, (2 * this.k2.fV() * this.theta2.fV()) / this.sigma22);

            double B1 = 2.0 * eg1 / den_A1;
            double B2 = 2.0 * eg2 / den_A2;

            double gg = discf * A1 * Math.Exp(-B1 * s1) * A2 * Math.Exp(-B2 * s2);
            double gd = RRate * discf + (1 - RRate) * gg;

            double y = -(1.0 / tau) * Math.Log(gd);
            double sd = y - ZR;

            double pdcum = (1 - Math.Exp(-sd * tau)) / (1.0 - RRate);

            return pdcum;
        }

        #endregion
    }
}
