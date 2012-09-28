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
using DVPLDOM;
using DVPLI;
using ParallelVectors;

namespace Pelsser
{
    /// <summary>
    /// Implementation of the Pelsser's Squared Gaussian model for interest rates.
    /// </summary>
    [Serializable]
    public class SquaredGaussianModel : IExtensibleProcessIR, IZeroRateReference, IMarkovSimulator,
                                        IParsable, IPostSimulationTransformation, IPopulable,
                                        IVectorialMarkovSimulator, IGreeksDerivativesInfo,
                                        IOpenCLCode
    {
        #region SerializedFields

        /// <summary>
        /// Reference to the zero rate.
        /// </summary>
        public IModelParameter zr;

        /// <summary>
        /// The mean reversion rate.
        /// </summary>
        public IModelParameter a1;

        /// <summary>
        /// The diffusion parameter.
        /// </summary>
        public IModelParameter sigma1;

        /// <summary>
        /// Unused but required for compatibility.
        /// </summary>
        private IModelParameter lamda0;

        #endregion

        /// <summary>
        /// Temporary cache for calculated values.
        /// </summary>
        [NonSerialized]
        private Dictionary<PelsserKey, PelsserCache> cache;

        /// <summary>
        /// Temporary zero rate function, used to optimize the simulation.
        /// </summary>
        [NonSerialized]
        private Function zeroRateCurve;

        /// <summary>
        /// DO NOT TOUCH AND DO NOT USE THIS member variable, it is kept only for compatibility.
        /// </summary>
        [Obsolete]
        private double[] mDates = null;

        /// <summary>
        /// Discretization of the alpha with the simulator step.
        /// </summary>
        [NonSerialized]
        private double[] alphaT0;

        /// <summary>
        /// Temporary variable for the transformation (mDailyDates).
        /// </summary>
        [NonSerialized]
        private double[] alphaT;

        /// <summary>
        /// Temporary variable to keep the precalculated values of C(t, T),
        /// for the various deltaT (from 0 to T).
        /// </summary>
        [NonSerialized]
        private double[] cDeltaT;

        /// <summary>
        /// Temporary variable to keep the precalculated values of D(t, T),
        /// for the various deltaT (from 0 to T).
        /// </summary>
        [NonSerialized]
        private double[] dDeltaT;

        /// <summary>
        /// Temporarily stores alpha1 to avoid recalculating it each time.
        /// </summary>
        [NonSerialized]
        private double alpha1Temp;

        /// <summary>
        /// Temporarily stores sigma1 to avoid recalculating it each time.
        /// </summary>
        [NonSerialized]
        private double sigma1Temp;

        /// <summary>
        /// Temporarily stores the squared version of sigma1 to avoid recalculating it each time.
        /// </summary>
        [NonSerialized]
        private double sigma1SquaredTemp;

        /// <summary>
        /// The gamma variable used in the Pelsser model to simplify formulas.
        /// </summary>
        [NonSerialized]
        private double gamma;

        /// <summary>
        /// Keeps a cached copy of the dates the simulation will be run with.
        /// </summary>
        [NonSerialized]
        private double[] cacheDates;

        /// <summary>
        /// Keeps the readable description of the zero rate.
        /// </summary>
        private static string zrDescription = "Zero Rate";

        /// <summary>
        /// Keeps the readable description of the a1 (alpha) model variable.
        /// </summary>
        private static string a1Description = "Alpha";

        /// <summary>
        /// Keeps the readable description of the sigma1 model variable.
        /// </summary>
        private static string sigma1Description = "Sigma";

        /// <summary>
        /// Gets or sets the cacheDates, which are the cached values of dates
        /// passed during the Setup.
        /// </summary>
        public double[] CacheDates
        {
            get
            {
                return this.cacheDates;
            }

            set
            {
                this.cacheDates = value;
            }
        }

        /// <summary>
        /// Default constructor. Builds a new SquaredGaussianModel with these default values:
        /// * sigma = 0.06.
        /// * lamba = -0.9.
        /// * alpha = 0.09.
        /// * no zero rate reference.
        /// </summary>
        public SquaredGaussianModel()
            : this(0.06, -0.9, 0.09, string.Empty)
        {
        }

        /// <summary>
        /// Constructor which builds a new SquareGaussianModel with provided values.
        /// </summary>
        /// <param name="sigma">The sigma factor of the model.</param>
        /// <param name="lambda">The parameter is not used.</param>
        /// <param name="a1">The alpha factor of the model.</param>
        /// <param name="zeroRateReference">A reference to a zero rate.</param>
        public SquaredGaussianModel(double sigma, double lambda, double a1, string zeroRateReference)
        {
            this.sigma1 = new ModelParameter(sigma, sigma1Description);
            this.a1 = new ModelParameter(a1, a1Description);
            this.zr = new ModelParameter(zeroRateReference, zrDescription);
        }

        /// <summary>
        /// Gets <see cref="alphaT"/> at the position i.
        /// </summary>
        /// <param name="i">The position where to get the alphaT value.</param>
        /// <returns>The requested value.</returns>
        public double AlphaT(int i)
        {
            return this.alphaT[i];
        }

        #region IParsable Members

        /// <summary>
        /// Ensure the parameters are correct.
        /// </summary>
        /// <param name='p_Context'>
        /// The underlying project.
        /// </param>
        /// <returns>
        /// False if there were no parse errors.
        /// </returns>
        public bool Parse(IProject p_Context)
        {
            bool errors = false;
            BoolHelper.AddBool(errors, this.a1.Parse(p_Context));
            BoolHelper.AddBool(errors, this.sigma1.Parse(p_Context));

            if (this.zr.Expression.IndexOf("@") == -1)
            {
                p_Context.AddError(this.zr.Expression + " is not a reference to a zero rate curve");
            }

            object zrReference = Engine.Parser.EvaluateAsReference(this.zr.Expression);
            if (!Engine.Parser.GetParserError())
            {
                this.zeroRateCurve = zrReference as Function;
                if (this.zeroRateCurve == null)
                {
                    errors = true;

                    p_Context.AddError("Cannot find the Zero Rate Curve! " + this.zr.Expression);
                }
            }
            else
                errors = true;

            if (!errors)
            {
                this.alpha1Temp = this.a1.fV();
                this.sigma1Temp = this.sigma1.fV();
                this.sigma1SquaredTemp = Math.Pow(this.sigma1Temp, 2);
                CalculateGamma();
            }

            return errors;
        }

        #endregion

        #region IExtensibleProcessIR Members

        /// <summary>
        /// Calculates the value of a Bond under the Pelsser model.
        /// </summary>
        /// <param name='dynamic'>
        /// The simulated process.
        /// </param>
        /// <param name='dates'>
        /// The vector of reference dates.
        /// </param>
        /// <param name='i'>
        /// The index at which the state variables must be sampled.
        /// </param>
        /// <param name='t'>
        /// The date in years/fractions at at which the state variables must be sampled.
        /// </param>
        /// <param name='s'>
        /// The maturity of the bond.
        /// </param>
        /// <returns>The value of the bound at index i using the Pelsser model.</returns>
        public double Bond(IReadOnlyMatrixSlice dynamic, double[] dates, int i, double t, double s)
        {
            // Get the value of the short rate.
            double y = Math.Sqrt(dynamic[i, 0]) - this.alphaT0[i];
            PelsserKey k = new PelsserKey(t, s);
            PelsserCache cachedValue = null;
            lock (this.cache)
            {
                if (this.cache.ContainsKey(k))
                    cachedValue = this.cache[k];
                else
                {
                    cachedValue = new PelsserCache(t, s, this);

                    // Insert the value in the cache.
                    this.cache.Add(k, cachedValue);
                }
            }

            double v = Math.Exp(cachedValue.A - y * cachedValue.B - (y * y) * cachedValue.CtT0);
            return v;
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
        /// Gets the ProcessInfo for this plug-in, in this case Pelsser Squared Gaussian Model.
        /// </summary>
        public ProcessInfo ProcessInfo
        {
            get
            {
                return new ProcessInfo("Pelsser Squared Gaussian Model");
            }
        }

        /// <summary>
        /// This will simply calculate the gamma factor used in the Pelsser
        /// formulation and set the relative temporary variable.
        /// </summary>
        /// <remarks>
        /// Gamma = sqrt(a^2 + 2 * sigma^2).
        /// </remarks>
        public void CalculateGamma()
        {
            this.gamma = Math.Sqrt(this.alpha1Temp * this.alpha1Temp + 2 * this.sigma1Temp * this.sigma1Temp);
        }

        /// <summary>
        /// Precalculates several values and functions of the model in order
        /// to cache them for later use.
        /// </summary>
        /// <param name="dates">
        /// The vector of dates which will be used to simulate the model.
        /// </param>
        internal void CalculateValueForCache(double[] dates)
        {
            this.CacheDates = new double[dates.Length];
            dates.CopyTo(this.CacheDates, 0);

            double dt = dates[1] - dates[0];
            double INT = Int(0, dates[0]);
            this.alphaT = new double[dates.Length];
            this.cDeltaT = new double[dates.Length];
            this.dDeltaT = new double[dates.Length];

            this.dDeltaT[0] = 1;
            this.cDeltaT[0] = 0;
            this.alphaT[0] = F(dates[0], dt) + 2 * Math.Exp(-this.alpha1Temp * dates[0]) * INT;

            DateTime t0 = DateTime.Now;

            for (int i = 1; i < dates.Length; i++)
            {
                double current = dates[i];

                // If it's not like this keep the last calculated value.
                if (i < dates.Length - 1)
                    dt = dates[i + 1] - current;
                INT += Int(dates[i - 1], current);

                this.alphaT[i] = F(current, dt) + 2 * Math.Exp(-this.alpha1Temp * current) * INT;
                this.cDeltaT[i] = C(current);
                this.dDeltaT[i] = D(current);
            }

            DateTime t1 = DateTime.Now;

            if (double.IsNaN(INT))
            {
                throw new Exception(string.Format("Pelsser model: error while initializing: a={0} s={1}, M={2}",
                                                  this.alpha1Temp, this.sigma1Temp, dates[dates.Length - 1]));
            }

            DateTime t2 = DateTime.Now;
        }

        /// <summary>
        /// Called by Simulator after parse.
        /// Initializes here time-dependant but not state dependent variables.
        /// </summary>
        /// <param name='dates'>
        /// The dates at which the process realizations will be requested.
        /// </param>
        public void Setup(double[] dates)
        {
            this.cache = new Dictionary<PelsserKey, PelsserCache>();

            this.CacheDates = null;

            // Do not touch, it's here for compatibility reasons.
            this.mDates = null;

            CalculateValueForCache(dates);
            this.alphaT0 = this.alphaT;

            return;
        }

        /// <summary>
        /// Gets the information required in order to allow the simulation to run.
        /// Pelsser Squared Gaussian has:
        /// * 0 latent components in the state components.
        /// * 1 component of noise.
        /// * 1 state component.
        /// * The unique component is named short rate.
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
            parameters.Add(this.a1);
            parameters.Add(this.sigma1);
            parameters.Add(this.zr);
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
        /// This function defines the drift in the Pelsser Markov process.
        /// The formula to calculate the A component is
        /// A = - alpha * previous State.
        /// </summary>
        /// <param name="i">The parameter is not used.</param>
        /// <param name="x">The state vector at the previous state.</param>
        /// <param name="a">The output of the function.</param>
        public unsafe void a(int i, double* x, double* a)
        {
            a[0] = -this.alpha1Temp * x[0];
        }

        /// <summary>
        /// This function defines the volatility in the Pelsser Markov process.
        /// The formula to calculate the B component is
        /// B = sigma.
        /// </summary>
        /// <param name="i">The parameter is not used.</param>
        /// <param name="x">The parameter is not used.</param>
        /// <param name="b">The output of the function.</param>
        public unsafe void b(int i, double* x, double* b)
        {
            b[0] = this.sigma1Temp;
        }

        /// <summary>
        /// Sets the passed array with a Boolean stating if the process
        /// must be simulated as a log-normal process.
        /// In this case it should not.
        /// </summary>
        /// <param name="isLog">
        /// A reference to the array to be set with the required information.
        /// </param>
        public void isLog(ref bool[] isLog)
        {
            isLog[0] = false;
        }

        /// <summary>
        /// Gets the starting point for the process y.
        /// In this case it starts from zero as the <see cref="Transform"/>
        /// will handle the conversion from y to the rate r.
        /// </summary>
        public double[] x0
        {
            get
            {
                return new double[] { 0 };
            }
        }

        #endregion

        #region IPostSimulationTransformation Implementation

        /// <summary>
        /// Handles the conversion, after the simulation, from y to the rate r.
        /// </summary>
        /// <param name="dates">The parameter is not used.</param>
        /// <param name="outDynamic">The input and output components of the transformation.</param>
        public void Transform(double[] dates, IMatrixSlice outDynamic)
        {
            for (int j = 0; j < dates.Length; j++)
            {
                outDynamic[j, 0] = Math.Pow(outDynamic[j, 0] + this.alphaT[j], 2);
            }
        }

        #endregion IPostSimulationTransformation Implementation

        #region IOpenCLCode implementation

        /// <summary>
        /// Gets the arguments needed for an OpenCL simulation.
        /// alpha1, sigma1 and alphaT are used in this context.
        /// </summary>
        public List<Tuple<string, object>> Arguments
        {
            get
            {
                List<Tuple<string, object>> args = new List<Tuple<string, object>>();
                args.Add(new Tuple<string, object>("alpha1Temp", this.alpha1Temp));
                args.Add(new Tuple<string, object>("sigma1Temp", this.sigma1Temp));
                args.Add(new Tuple<string, object>("alphaT", this.alphaT));
                return args;
            }
        }

        /// <summary>
        /// Gets the OpenCL code used to calculate A, B and the POSTRANFORM.
        /// </summary>
        public Dictionary<string, string> Code
        {
            get
            {
                Dictionary<string, string> sources = new Dictionary<string, string>();
                sources.Add("B", "*b = sigma1Temp;");
                sources.Add("A", "*a = -alpha1Temp*x[0];");
                sources.Add("POSTRANSFORM", "for(int j = 0; j < stepN; ++j) { x[j] = pown(x[j]+alphaT[j], 2); }");
                return sources;
            }
        }

        /// <summary>
        /// Gets a value indicating whether the plugin OpenCL implementation is usable.
        /// This plugin can always run through the OpenCL simulator so it always returns true.
        /// </summary>
        public bool IsOpenCLUsable
        {
            get
            {
                return true;
            }
        }

        #endregion

        /// <summary>
        /// Helper function to make functions easier to read.
        /// Just returns the value of the zero rate at position t.
        /// </summary>
        /// <param name="t">The position where to get the value of the zero rate from.</param>
        /// <returns>The value of the zero rate at position t.</returns>
        private double ZR(double t)
        {
            return this.zeroRateCurve.Evaluate(t);
        }

        /// <summary>
        /// Calculates the instantaneous forward.
        /// </summary>
        /// <param name="t">The time where to calculate the instantaneous forward.</param>
        /// <param name="dt">
        /// The delta of time to consider to calculate the instantaneous forward.
        /// </param>
        /// <returns>The instantaneous forward for the provided parameters.</returns>
        public double f(double t, double dt)
        {
            if (t == 0)
                return ZR(t);
            else
                return (ZR(t) * t - ZR(t - dt) * (t - dt)) / dt;
        }

        /// <summary>
        /// Calculates the C(t, T) function of the model.
        /// </summary>
        /// <param name="deltaT">The delta between T and t.</param>
        /// <returns>The result of the C function of the model.</returns>
        public double C(double deltaT)
        {
            return (Math.Exp(2 * this.gamma * (deltaT)) - 1) / ((this.alpha1Temp + this.gamma) * Math.Exp(2 * this.gamma * (deltaT)) + (this.gamma - this.alpha1Temp));
        }

        /// <summary>
        /// Calculates the D(t, T) function of the model.
        /// </summary>
        /// <param name="deltaT">The delta between T and t.</param>
        /// <returns>The result of the D function of the model.</returns>
        public double D(double deltaT)
        {
            return (2 * this.gamma * Math.Exp(this.gamma * (deltaT))) / ((this.alpha1Temp + this.gamma) * Math.Exp(2 * this.gamma * (deltaT)) + (this.gamma - this.alpha1Temp));
        }

        /// <summary>
        /// Calculates the sigma function(t, T) of the model.
        /// </summary>
        /// <param name="deltaT">The delta between T and t.</param>
        /// <returns>The result of the sigma function of the model.</returns>
        protected double SIG(double deltaT)
        {
            return this.sigma1Temp * this.sigma1Temp * C(deltaT);
        }

        /// <summary>
        /// Calculates F(T, dt) of the model, but doesn't apply the square root.
        /// This value is, as such equal to F-squared.
        /// </summary>
        /// <param name="T">The time where to calculate the function.</param>
        /// <param name="dt">
        /// The delta of time to consider to calculate the function.
        /// </param>
        /// <returns>The result of the F(T, dt) function.</returns>
        public double F2(double T, double dt)
        {
            return f(T, dt) - SIG(T);
        }

        /// <summary>
        /// Calculates F(T, dt) of the model.
        /// </summary>
        /// <param name="T">The time where to calculate the function.</param>
        /// <param name="dt">
        /// The delta of time to consider to calculate the function.
        /// </param>
        /// <returns>The result of the F(T, dt) function.</returns>
        public double F(double T, double dt)
        {
            return Math.Sqrt(F2(T, dt));
        }

        /// <summary>
        /// Calculates the integral inside alpha(T).
        /// </summary>
        /// <param name="t1">The starting point of the integration.</param>
        /// <param name="t2">The ending point of the integration.</param>
        /// <returns>The value of the integral.</returns>
        public double Int(double t1, double t2)
        {
            if (t1 == t2) return 0;
            double deltat = t2 - t1;
            double i = (Math.Exp(this.alpha1Temp * t1) * SIG(t1) * F(t1, deltat) + Math.Exp(this.alpha1Temp * t2) * SIG(t2) * F(t2, deltat)) / 2.0;
            return i * deltat;
        }

        /// <summary>
        /// Calculates the mu(t, T, y) used by the Caplets.
        /// </summary>
        /// <param name="T">The span of time to calculate the function on.</param>
        /// <returns>The result of the calculation.</returns>
        public double Mu0(double T)
        {
            if (T == 0) return 0;

            // Initialization of the integral.
            double integralMu = 0;
            int I = DVPLDOM.AdaptiveTimeDiscretization.DiscreteTime(T, this.CacheDates);
            double delta = this.CacheDates[1] - this.CacheDates[0];
            double[] b = B(0, I, delta);

            for (int i = 0; i < I; i++)
            {
                integralMu += this.dDeltaT[I - i] * b[i];
            }

            integralMu *= this.sigma1SquaredTemp;

            return -integralMu * delta;
        }

        /// <summary>
        /// Calculates the function A(t, T) of the model, which is an integral.
        /// </summary>
        /// <param name="ti">The starting Index from where to calculate.</param>
        /// <param name="si">The ending Index of the integration.</param>
        /// <param name="delta">The delta of time where to execute the calculation.</param>
        /// <param name="btT">
        /// A populated array with the results of the function B(s, T) of the model.
        /// </param>
        /// <returns>The result of the function A(t, T).</returns>
        public double A(int ti, int si, double delta, double[] btT)
        {
            // Initialization of the integral.
            double integralA = 0;
            int N = si - ti;

            for (int i = 0; i < N; i++)
            {
                integralA += 0.5 * this.sigma1SquaredTemp * Math.Pow(btT[i], 2) - this.sigma1SquaredTemp * this.cDeltaT[N - i] - Math.Pow(this.alphaT[ti + i], 2);
            }

            return integralA * delta;
        }

        /// <summary>
        /// Calculates the function B(t, T) of the model,
        /// which is an integral, for a single element.
        /// </summary>
        /// <param name="ti">The starting Index from where to calculate.</param>
        /// <param name="si">The ending Index of the integration.</param>
        /// <param name="delta">The delta of time where to execute the calculation.</param>
        /// <returns>The result of the function B(t, T).</returns>
        public double BSingle(int ti, int si, double delta)
        {
            // Initializes B integral.
            double integralB = 0;

            for (int i = ti; i < si; i++)
                integralB += this.alphaT[i] / this.dDeltaT[si - i];
            return 2 * this.dDeltaT[si - ti] * integralB * delta;
        }

        /// <summary>
        /// Calculates the function B(t, T) of the model,
        /// which is an integral, for all elements.
        /// </summary>
        /// <param name="ti">The starting Index from where to calculate.</param>
        /// <param name="si">The ending Index of the integration.</param>
        /// <param name="delta">The delta of time where to execute the calculation.</param>
        /// <returns>The result of the function B(t, T).</returns>
        public double[] B(int ti, int si, double delta)
        {
            double[] intermediates = new double[si - ti + 1];

            for (int i = ti; i < si; i++)
                intermediates[i - ti] = this.alphaT[i] / this.dDeltaT[si - i];

            double[] result = new double[si - ti];
            double IB = 0;
            for (int i = si - 1; i >= ti; i--)
            {
                IB += intermediates[i - ti];
                result[i - ti] = 2 * this.dDeltaT[si - i] * IB * delta;
            }

            return result;
        }

        /// <summary>
        /// Associate the process to a zero rate defined in the Fairmat model
        /// (e.g. @zr1).
        /// </summary>
        /// <param name='zr'>
        /// The zero rate reference.
        /// </param>
        public void SetZeroRateReference(string zr)
        {
            this.zr = new ModelParameter(zr, zrDescription);
        }

        /// <summary>
        /// Gets the zero rate reference.
        /// </summary>
        /// <returns>
        /// The zero rate reference.
        /// </returns>
        public string GetZeroRateReference()
        {
            return this.zr.Expression;
        }

        #region IPopulable Members

        /// <summary>
        /// Populate editable fields from name and value vectors
        /// specific to Pelsser.
        /// </summary>
        /// <param name="names">
        /// An array with the names of the variable,
        /// will search for alpha1 (or a), sigma1 (or sigma) and lamba0.
        /// </param>
        /// <param name="values">The values associated to the parameters in names.</param>
        public void Populate(string[] names, double[] values)
        {
            bool found;
            this.a1 = new ModelParameter(PopulateHelper.GetValue("alpha1", "a", names, values, out found), a1Description);
            this.sigma1 = new ModelParameter(PopulateHelper.GetValue("sigma1", "sigma", names, values, out found), sigma1Description);
        }

        #endregion

        #region IVectorialMarkovSimulator Members

        /// <summary>
        /// This function defines the drift in the Pelsser Markov process.
        /// The formula to calculate the A component is
        /// A = - alpha * previous State.
        /// This is the version which handles a vectorial execution.
        /// </summary>
        /// <param name="i">The parameter is not used.</param>
        /// <param name="x">The state Matrix at the previous state.</param>
        /// <param name="a">The output of the function.</param>
        public void va(int i, Matrix x, Matrix a)
        {
            VectorNP x1 = new VectorNP(x.GetRowReference(0));
            VectorNP tmp = -this.alpha1Temp * x1;
            tmp.CopyTo(a);
        }

        /// <summary>
        /// This function defines the volatility in the Pelsser markov process.
        /// The formula to calculate the B component is
        /// B = sigma.
        /// This is the version which handles a vectorial execution.
        /// </summary>
        /// <param name="i">The parameter is not used.</param>
        /// <param name="x">The parameter is not used.</param>
        /// <param name="b">The output of the function.</param>
        public void vb(int i, Matrix x, Matrix b)
        {
            VectorNP tmp = new VectorNP(b);
            b.Fill(this.sigma1Temp);
        }

        #endregion

        /// <summary>
        /// Returns the factors for Delta/Gamma Greek derivatives.
        /// In this case none so null is returned.
        /// </summary>
        /// <returns>Always null.</returns>
        public IModelParameter[] GetDeltaFactors()
        {
            return null;
        }

        /// <summary>
        /// Returns the factors for Vega Greek derivative.
        /// In this case just <see cref="sigma1"/>.
        /// </summary>
        /// <returns>
        /// An array of <see cref="IModelParameter"/> containing only <see cref="sigma1"/>
        /// </returns>
        public IModelParameter[] GetVegaFactors()
        {
            return new IModelParameter[] { this.sigma1 };
        }
    }
}
