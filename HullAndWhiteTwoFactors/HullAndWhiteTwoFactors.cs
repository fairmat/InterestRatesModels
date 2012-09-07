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
using ParallelVectors;

namespace HullAndWhiteTwoFactors
{
    /// <summary>
    /// Implements the Hull And White Two factors model simulation. For a general reference see
    /// "Numerical procedures for implementing term structure models II: Two-Factor Models",
    /// Journal of Derivatives, 2, 1 (Winter 1994b) 37- 48, referred in the comments as HW1994.
    /// </summary>
    [Serializable]
    public class HW2 : HW2Context, IExtensibleProcessIR, IZeroRateReference, IMarkovSimulator,
                       IParsable, IPopulable, IGreeksDerivativesInfo
    {
        #region Serialized Parameters

        /// <summary>
        /// Reference to a zero rate curve.
        /// </summary>
        public IModelParameter _zr;

        /// <summary>
        /// Mean-reversion rate for r (short rate).
        /// </summary>
        public IModelParameter _a1;

        /// <summary>
        /// Volatility of r  (short rate).
        /// </summary>
        public IModelParameter _a2;

        /// <summary>
        /// Mean-reversion rate for u (latent component).
        /// </summary>
        public IModelParameter _s1;

        /// <summary>
        /// Mean-reversion rate for r (short rate).
        /// </summary>
        public IModelParameter _s2;

        /// <summary>
        /// Mean-reversion rate for u (latent component).
        /// </summary>
        public IModelParameter _rho;

        /// <summary>
        /// Drift adjustment: can be used for adding risk premium or quanto adjustments.
        /// </summary>
        [OptionalField(VersionAdded = 2)]
        private IModelParameter driftAdjustment;

        #endregion Serialized Parameters

        /// <summary>
        /// A Data structure where pre-calculated values are stored.
        /// </summary>
        [NonSerialized]
        private SPCache cache;

        /// <summary>
        /// HW2 project context (available after calling Parse method).
        /// </summary>
        [NonSerialized]
        private Project context;

        /// <summary>
        /// Reference to the parsed zero rate function.
        /// </summary>
        [NonSerialized]
        private Function zeroRateCurve;

        /// <summary>
        /// Describe how the simulator should treat the HW2 process. Kept for compatibility.
        /// </summary>
        private static bool transformedSimulation = false;

        /// <summary>
        /// Cached value of Pow(sigma1, 2).
        /// </summary>
        [NonSerialized]
        private double sigma1Pow2;

        /// <summary>
        /// Cached value of alpha2 - alpha1.
        /// </summary>
        [NonSerialized]
        private double alpha2minus1;

        /// <summary>
        /// Pre-calculated theta. Available after Setup.
        /// </summary>
        [NonSerialized]
        private double[] theta;

        /// <summary>
        /// Simulated horizon, caches Project.GetTotalTime().
        /// </summary>
        [NonSerialized]
        private double horizon;

        /// <summary>
        /// Simulation dates. Available after setup.
        /// </summary>
        [NonSerialized]
        private double[] mDates = null;

        /// <summary>
        /// Keeps the readable description of the a1 model variable.
        /// </summary>
        private static string a1Description = "Alpha 1";

        /// <summary>
        /// Keeps the readable description of the a2 model variable.
        /// </summary>
        private static string a2Description = "Alpha 2";

        /// <summary>
        /// Keeps the readable description of the s1 model variable.
        /// </summary>
        private static string s1Description = "Sigma ";

        /// <summary>
        /// Keeps the readable description of the s2 model variable.
        /// </summary>
        private static string s2Description = "Sigma 2";

        /// <summary>
        /// Keeps the readable description of the zeroRate model variable.
        /// </summary>
        private static string zeroRateDescription = "Zero Rate";

        /// <summary>
        /// Keeps the readable description of the rho model variable.
        /// </summary>
        private static string rhoDescription = "rho";

        /// <summary>
        /// Keeps the readable description for drift adjustment.
        /// </summary>
        private const string driftAdjustmentDescription = "Drift Adjustment";

        /// <summary>
        /// Initializes a new instance of the <see cref="HullAndWhiteTwoFactors.HW2"/> class.
        /// This is the default constructor and sets the
        /// default values of several components of the model:
        /// * alpha 1 = 0.1
        /// * sigma 1 = 0.001
        /// * alpha 2 = 0.01
        /// * sigma 2 = 0.001
        /// * rho = 0.001
        /// * No zero rate.
        /// * drift correction = 0.
        /// </summary>
        public HW2()
        {
            this._zr = new ModelParameter(string.Empty, zeroRateDescription);
            this._a1 = new ModelParameter(0.1, a1Description);
            this._s1 = new ModelParameter(0.001, s1Description);
            this._a2 = new ModelParameter(0.01, a2Description);
            this._s2 = new ModelParameter(0.001, s2Description);
            this._rho = new ModelParameter(0.001, rhoDescription);
            this.driftAdjustment = new ModelParameter(0, driftAdjustmentDescription);
        }

        /// <summary>
        /// Initializes optional fields after deserialization.
        /// </summary>
        /// <param name='context'>
        /// The parameter is not used.
        /// </param>
        [OnDeserialized]
        public void OnDeserialized(StreamingContext context)
        {
            if (this.driftAdjustment == null)
            {
                this.driftAdjustment = new ModelParameter(0, driftAdjustmentDescription);
            }
        }

        #region IPopulable Members

        /// <summary>
        /// Populate editable fields from name and value vectors
        /// specific to HW2.
        /// </summary>
        /// <param name="names">
        /// An array with the names of the variable,
        /// will search for alpha1 (or a1), sigma1 (or sigma), alpha2 (or a2), sigma2 and rho.
        /// </param>
        /// <param name="values">The values associated to the parameters in names.</param>
        public void Populate(string[] names, double[] values)
        {
            bool found = false;
            this._a1 = new ModelParameter(PopulateHelper.GetValue("alpha1", "a1", names, values, out found), a1Description);
            this._s1 = new ModelParameter(PopulateHelper.GetValue("sigma1", "sigma", names, values, out found), s1Description);
            this._a2 = new ModelParameter(PopulateHelper.GetValue("alpha2", "a2", names, values, out found), a2Description);
            this._s2 = new ModelParameter(PopulateHelper.GetValue("sigma2", names, values, out found), s2Description);
            this._rho = new ModelParameter(PopulateHelper.GetValue("rho", names, values, out found), rhoDescription);
        }

        #endregion

        /// <summary>
        /// Evaluates the referenced zr.
        /// </summary>
        /// <param name="t">The evaluation date.</param>
        /// <returns>The evaluated zr.</returns>
        private double Zr(double t)
        {
            return this.zeroRateCurve.Evaluate(t);
        }

        #region IZeroRateReference members

        /// <summary>
        /// Associate the process to a zero rate defined in the Fairmat model
        /// (e.g. @zr1).
        /// </summary>
        /// <param name='zr'>
        /// The zero rate reference.
        /// </param>
        public void SetZeroRateReference(string zr)
        {
            this._zr = new ModelParameter(zr, "Zero Rate");
        }

        /// <summary>
        /// Gets the zero rate reference.
        /// </summary>
        /// <returns>
        /// The zero rate reference.
        /// </returns>
        public string GetZeroRateReference()
        {
            return this._zr.Expression;
        }

        #endregion

        /// <summary>
        /// Called by Simulator after parse.
        /// Initializes here time-dependant but not state dependent variables.
        /// </summary>
        /// <param name='dates'>
        /// The dates at which the process realizations will be requested.
        /// </param>
        public void Setup(double[] dates)
        {
            Preprocessing();

            this.mDates = dates;
            this.cache = new SPCache();
            int length = dates.Length;

            this.theta = new double[length];
            double dt = dates[1] - dates[0];

            for (int i = 0; i < length; i++)
            {
                this.theta[i] = Theta(dates[i], dt);
                if (i < length - 1)
                    dt = dates[i + 1] - dates[i];
            }
        }

        /// <summary>
        /// Parses the data and ensures the parameters are correct.
        /// </summary>
        /// <param name='p_Context'>
        /// The underlying project.
        /// </param>
        /// <returns>
        /// True if the parameter are correct.
        /// </returns>
        public bool Parse(IProject p_Context)
        {
            this.context = p_Context as Project;
            bool errors = false;

            BoolHelper.AddBool(errors, this._a1.Parse(p_Context));
            BoolHelper.AddBool(errors, this._a2.Parse(p_Context));

            BoolHelper.AddBool(errors, this._s1.Parse(p_Context));
            BoolHelper.AddBool(errors, this._s2.Parse(p_Context));
            BoolHelper.AddBool(errors, this._rho.Parse(p_Context));

            BoolHelper.AddBool(errors, this.driftAdjustment.Parse(p_Context));

            if (this._zr.Expression.IndexOf("@") == -1)
            {
                p_Context.AddError(this._zr.Expression +
                                   " is not a reference to a zero rate curve");
            }

            // Checks for the model constraints: alpha1 != alhpa2
            if (Math.Abs(this._a1.fV() - this._a2.fV()) < 10e-5)
            {
                p_Context.AddError("H&W2:  alpha1 and alpha2 must be different");
            }

            object zr_reference = Engine.Parser.EvaluateAsReference(this._zr.Expression);
            if (!Engine.Parser.GetParserError())
            {
                this.zeroRateCurve = zr_reference as Function;
                if (this.zeroRateCurve == null)
                {
                    errors = true;
                    p_Context.AddError("Cannot find the Zero Rate Curve! " + this._zr.Expression);
                }
            }
            else
                errors = true;

            if (!errors)
            {
                base.alpha1 = this._a1.fV();
                base.sigma1 = this._s1.fV();
                this.sigma1Pow2 = System.Math.Pow(this._s1.fV(), 2);
            }

            CorrelationMatrix R = (p_Context as ProjectProcess).Processes.r;
            int index = (p_Context as ProjectProcess).Processes.GetProcessCorrelationIndex(this);

            // Index is -1 is when the process is not still in the process list.
            if (index != -1)
            {
                // Updates the correlation in the global correlation matrix.
                R.Set(index, index + 1, this._rho);
            }

            return errors;
        }

        /// <summary>
        /// Gets the ProcessInfo for this process. In this case H&amp;W2.
        /// </summary>
        public ProcessInfo ProcessInfo
        {
            get
            {
                ProcessInfo pi = new ProcessInfo("H&W2");
                return pi;
            }
        }

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
            List<IExportable> list = new List<IExportable>();
            list.Add(this._zr);
            list.Add(this._a1);
            list.Add(this._s1);
            list.Add(this._a2);
            list.Add(this._s2);
            list.Add(this._rho);
            list.Add(this.driftAdjustment);
            return list;
        }

        /// <summary>
        /// Calculates and caches the model parameters for use during the simulation.
        /// </summary>
        private void CalculateInnerModelParameters()
        {
            double k12 = 0; // Proxy for rho
            // Get the rho from the correlation matrix.
            if (!GenericSolverParameters.m_ConsiderCorrMatrixAsCholesky)
                base.rho = this._rho.fV();
            else
                k12 = this._rho.fV();

            base.alpha1 = this._a1.fV();
            base.sigma1 = this._s1.fV();
            base.alpha2 = this._a2.fV();
            if (!GenericSolverParameters.m_ConsiderCorrMatrixAsCholesky)
            {
                base.sigma2 = this._s2.fV();
            }
            else
            {
                double k2 = this._s2.fV();
                base.sigma2 = Math.Sqrt(k12 * k12 + k2 * k2);
                base.rho = k12 / base.sigma2;
            }

            base.alpha2 = this._a2.fV();
            base.sigma2 = this._s2.fV();
        }

        /// <summary>
        /// Does a preliminary processing of all the parameters
        /// and required data which will be used during the simulation.
        /// </summary>
        private void Preprocessing()
        {
            CalculateInnerModelParameters();

            // Pre-calculate the parameters.
            this.horizon = this.context.GetTotalTime();

            // Pre-calculate the theta.
            int length = this.context.timeDiscretization.Periods;
            ITimeDiscretization td = this.context.timeDiscretization;
            this.theta = new double[length];
            for (int i = 0; i < length; i++)
                this.theta[i] = Theta(td.ContinuousTime(i), td.dt(i));
        }

        /// <summary>
        /// Calculates the standard normal cumulative distribution function.
        /// </summary>
        /// <param name="t">The date.</param>
        /// <returns>The value of the phi in the model.</returns>
        private double Phi(double t)
        {
            double phi = (1 - Math.Exp(-base.alpha1 * t)) / base.alpha1;
            double gamma = Math.Exp(-base.alpha1 * t) / (base.alpha1 * (-this.alpha2minus1)) -
                           Math.Exp(-base.alpha2 * t) / (base.alpha2 * (-this.alpha2minus1)) +
                           1.0 / (base.alpha1 * base.alpha2);

            return 0.5 * (Math.Pow(base.sigma1 * phi, 2.0) + Math.Pow(base.sigma2 * gamma, 2.0)) +
                   base.rho * base.sigma1 * base.sigma2 * phi * gamma;
        }

        /// <summary>
        /// Calculates either the right numerical derivative or
        /// the central numerical derivative of Phi.
        /// </summary>
        /// <param name="t">The date.</param>
        /// <param name="dt">The numerical derivative step.</param>
        /// <returns>The calculated value.</returns>
        private double DPhi(double t, double dt)
        {
            if (t == 0)
                return (Phi(dt) - Phi(0)) / dt;

            if (t == this.horizon)
                return (Phi(this.horizon) - Phi(this.horizon - dt)) / dt;

            return (Phi(t + dt) - Phi(t - dt)) / (2 * dt);
        }

        /// <summary>
        /// Calculates the theta element of the Hull And White formula which will be
        /// stored in the theta in order to be used during simulation.
        /// </summary>
        /// <param name="t">The position in which this value will be calculated.</param>
        /// <param name="dt">The delta between this t position and the previous one.</param>
        /// <returns>The requested value of theta, determined by the input values.</returns>
        protected double Theta(double t, double dt)
        {
            return DF(t, dt) + DPhi(t, dt) + base.alpha1 * (F(t, dt) + Phi(t));
        }

        /// <summary>
        /// Numerically calculates the instantaneous forward rate.
        /// </summary>
        /// <param name='t'>
        /// Time at which calculate the forward rate.
        /// </param>
        /// <param name='dt'>
        /// Interval to be used in the numerical derivative.
        /// </param>
        /// <returns>
        /// The value of the instantaneous forward rate.
        /// </returns>
        protected double F(double t, double dt)
        {
            if (t == 0)
                return Zr(t);
            else
                return (Zr(t) * t - Zr(t - dt) * (t - dt)) / dt;
        }

        /// <summary>
        /// Numerically calculates the derivative of function F().
        /// </summary>
        /// <param name='t'>
        /// Time at which calculate the derivative.
        /// </param>
        /// <param name='dt'>
        /// Semi-interval to be used in the numerical derivative.
        /// </param>
        /// <returns>
        /// The derivative of the function F().
        /// </returns>
        protected double DF(double t, double dt)
        {
            if (t == 0)
                return (F(dt, dt) - F(0, dt)) / dt;

            if (t == this.horizon)
                return (F(this.horizon, dt) - F(this.horizon - dt, dt)) / dt;

            return (F(t + dt, dt) - F(t - dt, dt)) / (2 * dt);
        }

        /// <summary>
        /// Vectorial version of method a.
        /// </summary>
        /// <param name="i">The discrete time-step.</param>
        /// <param name="x">The actual state matrix.</param>
        /// <param name="a">The output drift matrix.</param>
        public unsafe void va(int i, Matrix x, Matrix a)
        {
            // Gets the reference to the state and drift components.
            VectorNP x1 = new VectorNP(x.GetRowReference(0));
            VectorNP x2 = new VectorNP(x.GetRowReference(1));

            VectorNP a1 = new VectorNP(a.GetRowReference(0));
            VectorNP a2 = new VectorNP(a.GetRowReference(1));

            if (transformedSimulation)
            {
                VectorNP delta_r = this.theta[i] - alpha1 * x1;
                if (DVPLI.SolverAssumptions.UseRiskNeutralMeasure)
                    delta_r += this.driftAdjustment.vfV();

                delta_r.CopyTo(a1);
            }
            else
            {
                // Does a straight simulation.
                VectorNP delta_r = this.theta[i] + x2 - alpha1 * x1;
                delta_r.CopyTo(a2);
            }

            VectorNP delta_r2 = -alpha2 * x1;
            delta_r2.CopyTo(a2);
        }

        #region IBond Members

        /// <summary>
        /// Calculates the value of a Bond under the Hull and White Two factors model.
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
        /// <returns>The value of the bond at index i using the HW2 model.</returns>
        public double Bond(IReadOnlyMatrixSlice dynamic, double[] dates, int i, double t, double s)
        {
            double R = dynamic[i, 0];
            double U = dynamic[i, 1];

            double dt = this.context.timeDiscretization.dt(i);
            double Ps = ZCB(s);
            double Pt = ZCB(t);
            double Pdt = ZCB(t + dt);

            if (Engine.Log)
                if (SolverContext.SimulationRealization == 0)
                    Log.Write(string.Format("BOND i={0} Ps={1} Pt={2}", i, Ps, Pt));

            return BondHW2(this.context, this.cache, this, R, U, Ps, Pt, Pdt, s, t, dt);
        }

        /// <summary>
        /// Calculates the value of a Zero Coupon Bond.
        /// </summary>
        /// <param name="t">The bond maturity.</param>
        /// <returns>The bond value.</returns>
        private double ZCB(double t)
        {
            return Math.Exp(-this.zeroRateCurve.Evaluate(t) * t);
        }

        #endregion

        /// <summary>
        /// Price of a zcb at t > 0 with maturity TT > t
        /// as a function of the rate R, of the factor u at t
        /// and of the current zr curve (ZR)
        /// using the Hull &amp; White two-factor model
        /// <para>
        /// dr = (theta(t) + u - a * r) * dt + sigma1 * dz1
        /// du = -b * du * dt + sigma2 * dz2
        /// with dz1 * dz2 = rho * dt
        /// </para>
        /// (see Hull-White (1994) Journal of Derivatives).
        /// </summary>
        /// <param name="context">The underlying project context.</param>
        /// <param name="cache">Data structure where pre-calculated values are stored.</param>
        /// <param name="hw2">An HW2Context object to use for the evaluation.</param>
        /// <param name="R">Rate of return from time t to t+dt.</param>
        /// <param name="u">Factor u at t.</param>
        /// <param name="PT">Current price of a zcb maturing at TT.</param>
        /// <param name="P_t">Current price of a zcb maturing at t.</param>
        /// <param name="Pdt">Current price of a zcb maturing at t + dt.</param>
        /// <param name="TT">The bond maturity.</param>
        /// <param name="t">The initial date of the zcb.</param>
        /// <param name="dt">The length of time steps.</param>
        /// <returns>The value of the bond using the HW2 model.</returns>
        private static double BondHW2(Project context, SPCache cache, HW2Context hw2, double R, double u, double PT, double P_t, double Pdt, double TT, double t, double dt)
        {
            double a_a_b = +hw2.alpha1 * (hw2.alpha1 - hw2.alpha2);
            double b_a_b = +hw2.alpha2 * (hw2.alpha1 - hw2.alpha2);
            double ab = hw2.alpha1 * hw2.alpha2;

            double TT_t = TT - t;

            int t_idx = context.DiscreteTime(t);
            int TT_idx = context.DiscreteTime(TT);

            double aHat = 0;
            double bHat = 0;
            double cHat = 0;
#if _USE_CACHE_

#if !_SWAPTIONW_
            lock (cache)
#endif
            {
                if (!cache.IsCached_Eta[t_idx])
                {
                    cache.Cached_Eta[t_idx] = Eta(t, t + dt, a, b, sigma1, sigma2, rho);
                    cache.IsCached_Eta[t_idx] = true;
                }
            }

            int TT_idx_rel= TT_idx-t_idx;
            if(cache.cache.J>TT_idx_rel)
            {
                HWCache2 cached_element = null;
#if !_SWAPTIONW_
                lock (cache.cache)
#endif
                {
                    if (!cache.cache.IsCached(t_idx, TT_idx_rel))
                    {

                        cached_element = cache.cache.InsertToCache(t_idx, TT_idx_rel);
                        cached_element.Eta2 = Eta(t, TT, a, b, sigma1, sigma2, rho);
                        cached_element.Bhat = BHat(t, TT, dt, a);
                        cached_element.Chat = CHat(t, TT, dt, a, b, cached_element.Bhat);

                        double BtT = (1 - Math.Exp(-a * (TT_t))) / a;
                        double Btdt = (1 - Math.Exp(-a * (dt))) / a;
                        cached_element.Ahat = PT / P_t * Math.Exp(-BtT / Btdt * Math.Log(Pdt / P_t)
                                              - cached_element.Eta2
                                              + BtT / Btdt * cache.Cached_Eta[t_idx]);

                }
            }

            cached_element=cache.cache.Get(t_idx,TT_idx_rel);
            bHat=cached_element.Bhat;
            cHat=cached_element.Chat;
            aHat=cached_element.Ahat;
        }
        else
        {
            double BtT = (1 - Math.Exp(-a*(TT_t)))/a;
            double Btdt = (1 - Math.Exp(-a*(dt)))/a;

            aHat = PT/P_t *Math.Exp( -BtT/Btdt * Math.Log( Pdt/P_t )
                   - Eta(t,TT,a,b,sigma1,sigma2,rho)
                   + BtT/Btdt * cache.Cached_Eta[t_idx]);

            bHat=BHat(t,TT,dt,a);
            cHat=CHat(t,TT,dt,a,b,bHat);
        }
#else
            double BtT = (1 - Math.Exp(-hw2.alpha1 * (TT_t))) / hw2.alpha1;
            double Btdt = (1 - Math.Exp(-hw2.alpha1 * (dt))) / hw2.alpha1;

            double l_Pdt_Pt = Math.Log(Pdt / P_t);
            double eta_TT = hw2.Eta(t, TT);
            double eta_dt = hw2.Eta(t, t + dt);

            double arg = -BtT / Btdt * (l_Pdt_Pt - eta_dt) - eta_TT;

            if (Engine.Log)
                if (SolverContext.SimulationRealization == 0)
                    Log.Write(string.Format("arg {0} eta_TT {1} eta_dt {2}", arg, eta_TT, eta_dt));
            aHat = (PT / P_t) * Math.Exp(arg);

            bHat = hw2.BHat(t, TT, dt);
            cHat = hw2.CHat(t, TT, dt, bHat);
#endif
            double bond = aHat * Math.Exp(-bHat * R - cHat * u);

#if     _CHECK_OWERFLOWS_
        if(Double.IsInfinity(bond))
            bond= Double.MaxValue/1000000;
#endif
            return bond;
        }

        #region IMarkovSimulator Members

        /// <summary>
        /// Gets details about the structure of the functions A and B of the Markov
        /// process.
        /// In this case drift and volatility are only state dependant and not time
        /// dependant.
        /// </summary>
        public DynamicInfo DynamicInfo
        {
            get
            {
                DynamicInfo di = new DynamicInfo(false, true, false, false);
                return di;
            }
        }

        /// <summary>
        /// Sets the passed array with a Booleans stating if the process
        /// components must be simulated as a log-normal process.
        /// </summary>
        /// <param name="isLog">
        /// A reference to the array to be set with the required information.
        /// </param>
        public void isLog(ref bool[] isLog)
        {
            isLog[0] = false;
            isLog[1] = false;
        }

        /// <summary>
        /// This function calculates the drift in the HW2 Markov process.
        /// </summary>
        /// <param name="i">The time step of the simulation.</param>
        /// <param name="x">The state vector at the previous state.</param>
        /// <param name="a">The drift,output of the function.</param>
        public unsafe void a(int i, double* x, double* a)
        {
            if (transformedSimulation)
            {
                // Does a transformed simulation.
                double delta_r = this.theta[i] - base.alpha1 * x[0] + this.driftAdjustment.fV();
                a[0] = delta_r;
            }
            else
            {
                // Does a straight simulation.
                double delta_r = this.theta[i] + x[1] - base.alpha1 * x[0] + this.driftAdjustment.fV();
                a[0] = delta_r;
            }

            // Unobservable factor.
            a[1] = -base.alpha2 * x[1];
        }

        /// <summary>
        /// This function updates the volatility parameters (vector b) of the HW2 process.
        /// </summary>
        /// <param name="i">The parameter is not used.</param>
        /// <param name="x">The parameter is not used.</param>
        /// <param name="b">The output of the function.</param>
        public unsafe void b(int i, double* x, double* b)
        {
            b[0] = base.sigma1;
            b[1] = base.sigma2;
        }

        /// <summary>
        /// Gets the two-components starting points.
        /// </summary>
        public double[] x0
        {
            get
            {
                double[] s = new double[2];
                double dt = this.mDates[1];
                s[0] = Zr(dt);
                return s;
            }
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
        /// Gets the information required for the simulator in order to control
        /// the plug-in through the <see cref="IMarkovSimulator"/> interface.
        /// Hull And White two factors has:
        /// * 1 latent components in the state components.
        /// * 2 components of noise.
        /// * 2 state components.
        /// </summary>
        public SimulationInfo SimulationInfo
        {
            get
            {
                SimulationInfo s = new SimulationInfo();
                s.LatentSize = 1;
                s.NoiseSize = 2;
                s.StateDescription = new string[] { "short rate", "latent component" };
                s.StateSize = 2;
                return s;
            }
        }

        #endregion

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
        /// A vector containing the sigma of the observable and non observable processes.
        /// </returns>
        public IModelParameter[] GetVegaFactors()
        {
            return new IModelParameter[] { this._s1, this._s2 };
        }
    }
}
