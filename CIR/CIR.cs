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
using DVPLDOM;
using DVPLI;
using Mono.Addins;

namespace CIRProcess
{
    /// <summary>
    /// Implementation of one Cox-Ingersoll-Ross model.
    /// </summary>
    [Serializable]
    public class CIR : IExtensibleProcessIR, IMarkovSimulator, IParsable,
                       IEstimationResultPopulable, IGreeksDerivativesInfo
    {
        #region SerializedFields

        /// <summary>
        /// The k parameter of the model. It's the mean reversion speed.
        /// </summary>
        public IModelParameter k;

        /// <summary>
        /// The theta parameter of the model. It's the long term mean.
        /// </summary>
        public IModelParameter theta;

        /// <summary>
        /// The sigma parameter of the model. It's the volatility.
        /// </summary>
        public IModelParameter sigma;

        /// <summary>
        /// The starting point for the model.
        /// </summary>
        public IModelParameter r0;
        #endregion

        /// <summary>
        /// Keeps the name of the implementation.
        /// </summary>
        internal static string CIRLabel = "CIR";

        /// <summary>
        /// Temporary evaluated version of <see cref="k"/>, used to speed up calculation.
        /// </summary>
        [NonSerialized]
        private double alphaTemp;

        /// <summary>
        /// Temporary evaluated version of <see cref="theta"/>, used to speed up calculation.
        /// </summary>
        [NonSerialized]
        private double gammaTemp;

        /// <summary>
        /// Temporary evaluated version of <see cref="sigma"/>, used to speed up calculation.
        /// </summary>
        [NonSerialized]
        private double sigmaTemp;

        /// <summary>
        /// Temporary value used in many calculations <see cref="Setup"/> method.
        /// </summary>
        [NonSerialized]
        private double d;

        /// <summary>
        /// Temporary value used in many calculations <see cref="Setup"/> method.
        /// </summary>
        [NonSerialized]
        private double nu;

        /// <summary>
        /// A list of the parameter names to be shown to the user.
        /// </summary>
        internal static string[] parameterNames = { "k", "Theta", "Sigma", "r0" };

        /// <summary>
        /// Parameter values to be stored to load them in the model.
        /// </summary>
        [NonSerialized]
        private double[] parameterValues = new double[4];

        /// <summary>
        /// Default constructor, sets all the values to a default of 0.001.
        /// </summary>
        public CIR()
        {
            this.parameterValues = new double[] { 0.001, 0.001, 0.001, 0.001 };
            SetParametersValue();
        }

        /// <summary>
        /// Prepares the serialized model parameters with the provided values.
        /// </summary>
        private void SetParametersValue()
        {
            this.k = new ModelParameter(this.parameterValues[0], parameterNames[0]);
            this.theta = new ModelParameter(this.parameterValues[1], parameterNames[1]);
            this.sigma = new ModelParameter(this.parameterValues[2], parameterNames[2]);
            this.r0 = new ModelParameter(this.parameterValues[3], parameterNames[3]);
        }

        #region IParsable Members

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
            BoolHelper.AddBool(errors, this.k.Parse(context));
            BoolHelper.AddBool(errors, this.theta.Parse(context));
            BoolHelper.AddBool(errors, this.sigma.Parse(context));
            BoolHelper.AddBool(errors, this.r0.Parse(context));

            if (!errors)
            {
                this.alphaTemp = this.k.fV();
                this.gammaTemp = this.theta.fV();
                this.sigmaTemp = this.sigma.fV();
            }

            return errors;
        }

        #endregion

        #region IExtensibleProcessIR Members

        /// <summary>
        /// Calculates the value of a Bond under the Cox-Ingersoll-Ross model.
        /// </summary>
        /// <param name='dynamic'>
        /// The simulated process.
        /// </param>
        /// <param name='dates'>
        /// The vector of reference dates.
        /// </param>
        /// <param name='i'>
        /// The index at which the state variable must be sampled.
        /// </param>
        /// <param name='t'>
        /// The date in years/fractions at at which the state variable must be sampled.
        /// </param>
        /// <param name='s'>
        /// The maturity of the bond.
        /// </param>
        /// <returns>
        /// The value of the bound at index i using the Cox-Ingersoll-Ross model.
        /// </returns>
        public double Bond(IReadOnlyMatrixSlice dynamic, double[] dates, int i, double t, double s)
        {
            double r = dynamic[i, 0];
            double T = s - t;
            double den = (this.alphaTemp + this.d) * (Math.Exp(this.d * T) - 1.0) + 2.0 * this.d;
            double A = Math.Pow(2.0 * this.d * Math.Exp(0.5 * (this.alphaTemp + this.d) * T) / den, this.nu);
            double B = 2.0 * (Math.Exp(this.d * T) - 1.0) / den;
            return A * Math.Exp(-r * B);
        }

        #endregion

        #region IExtensibleProcess Members

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
        /// Gets the ProcessInfo for this plugin, in this case CIR.
        /// </summary>
        public ProcessInfo ProcessInfo
        {
            get
            {
                return new ProcessInfo(CIRLabel);
            }
        }

        /// <summary>
        /// Called by Simulator after parse.
        /// Calculates the variables d and nu for faster execution.
        /// </summary>
        /// <param name='dates'>The parameter is not used.</param>
        public void Setup(double[] dates)
        {
            this.d = Math.Sqrt(this.alphaTemp * this.alphaTemp + 2.0 * this.sigmaTemp * this.sigmaTemp);
            this.nu = 2.0 * this.alphaTemp * this.gammaTemp / (this.sigmaTemp * this.sigmaTemp);
        }

        /// <summary>
        /// Gets the information required in order to allow the simulation to run.
        /// CIR has:
        /// * 0 latent components in the state components.
        /// * 1 component of noise.
        /// * 1 state component.
        /// * The unique component is named short rate.
        /// * No privileged component.
        /// </summary>
        public SimulationInfo SimulationInfo
        {
            get
            {
                SimulationInfo s = new SimulationInfo();
                s.LatentSize = 0;
                s.NoiseSize = 1;
                s.StateDescription = new string[] { "short rate" };
                s.StateSize = 1;
                s.DefaultComponent = -1;
                return s;
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
            List<IExportable> parameters = new List<IExportable>();
            parameters.Add(this.k);
            parameters.Add(this.theta);
            parameters.Add(this.sigma);
            parameters.Add(this.r0);
            return parameters;
        }

        #endregion

        #region IMarkovSimulator Members

        /// <summary>
        /// Gets details about the structure of the functions A and B of the Markov
        /// process.
        /// In this case drift is state dependant and not time dependant and
        /// volatility is neither state or time dependent.
        /// </summary>
        public DynamicInfo DynamicInfo
        {
            get
            {
                return new DynamicInfo(false, true, false, false);
            }
        }

        /// <summary>
        /// This function defines the drift in the Cox-Ingersoll-Ross Markov process.
        /// The formula to calculate the A component is
        /// A = k * (theta - max(previous state, 0)).
        /// </summary>
        /// <param name="i">The parameter is not used.</param>
        /// <param name="x">The state vector at the previous state.</param>
        /// <param name="a">The output of the function.</param>
        public unsafe void a(int i, double* x, double* a)
        {
            a[0] = this.alphaTemp * (this.gammaTemp - Math.Max(x[0], 0));
        }

        /// <summary>
        /// This function defines the volatility in the Cox-Ingersoll-Ross Markov process.
        /// The formula to calculate the B component is
        /// B = sigma * sqrt(max(previous state, 0)).
        /// </summary>
        /// <param name="i">The parameter is not used.</param>
        /// <param name="x">The state vector at the previous state.</param>
        /// <param name="b">The output of the function.</param>
        public unsafe void b(int i, double* x, double* b)
        {
            b[0] = this.sigmaTemp * Math.Sqrt(Math.Max(x[0], 0));
        }

        /// <summary>
        /// This function calculates the drift and the volatility in the Cox-Ingersoll-Ross Markov process.
        /// as 
        /// A = k * (theta - max(previous state, 0)).
        /// B = sigma * sqrt(max(previous state, 0)).
        /// </summary>
        /// <param name="i">The parameter is not used.</param>
        /// <param name="x">The state vector at the previous state.</param>
        /// <param name="b">The output of the function.</param>
        public unsafe void ab(int i, double* x, double* a,double* b)
        {
            a[0] = this.alphaTemp * (this.gammaTemp - Math.Max(x[0], 0));
            b[0] = this.sigmaTemp * Math.Sqrt(Math.Max(x[0], 0));
        }



        /// <summary>
        /// Sets the passed array with a Boolean stating if the process
        /// must be simulated as a log-normal process.
        /// </summary>
        /// <param name="isLog">
        /// A reference to the array to be set with the required information.
        /// The array contains a single member set to false as this simulation
        /// is not log-normal.
        /// </param>
        public void isLog(ref bool[] isLog)
        {
            isLog[0] = false;
        }

        /// <summary>
        /// Gets the starting point for the process, which is defined in <see cref="r0"/>.
        /// </summary>
        public double[] x0
        {
            get
            {
                return new double[] { this.r0.fV() };
            }
        }

        #endregion

        /// <summary>
        /// Populate editable fields from name and value vectors
        /// specific to the CIR process.
        /// </summary>
        /// <param name="container">
        /// The stochastic process which is being referenced to.
        /// </param>
        /// <param name="estimate">
        /// The estimation result which contains values and names of parameters.
        /// </param>
        public void Populate(IStochasticProcess container, EstimationResult estimate)
        {
            bool found;
            this.parameterValues = new double[4];
            for (int i = 0; i < parameterNames.Length; i++)
            {
                this.parameterValues[i] = PopulateHelper.GetValue(parameterNames[i],
                                                                  estimate.Names, estimate.Values,
                                                                  out found);
                Console.WriteLine("ParameterValues[{0}] = {1}\t ParameterNames[{0}] = {2}",
                                  i, this.parameterValues[i], parameterNames[i]);
            }

            SetParametersValue();
        }

        #region IGreeksDerivativesInfo implementation

        /// <summary>
        /// Gets the factors for Delta Greek derivative.
        /// </summary>
        /// <returns>
        /// Null as the functionality is not implemented.
        /// </returns>
        public IModelParameter[] GetDeltaFactors()
        {
            // TODO: fixme how to handle short rate processes?
            return null;
        }

        /// <summary>
        /// Gets the factors for Vega Greek derivative.
        /// </summary>
        /// <returns>
        /// A model parameter containing the sigma value of HW.
        /// </returns>
        public IModelParameter[] GetVegaFactors()
        {
            return new IModelParameter[] { this.sigma };
        }

        #endregion IGreeksDerivativesInfo implementation
    }
}
