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
using System.IO;
using System.Text;
using DVPLDOM;
using DVPLI;
using Fairmat.Finance;
using Fairmat.Optimization;
using Fairmat.MarketData;
using Mono.Addins;


namespace Pelsser.Calibration
{
    /// <summary>
    /// Implementation of Pelsser Calibration against caps matrices.
    /// </summary>
    [Extension("/Fairmat/Estimator")]
    public class CapletEstimator : IEstimatorEx, IDescription, IMenuItemDescription
    {
        #region IEstimator Members

        /// <summary>
        /// Gets the value requested by the interface ProvidesTo,
        /// returning SquaredGaussianModel as the type.
        /// </summary>
        public Type ProvidesTo
        {
            get
            {
                return typeof(Pelsser.SquaredGaussianModel);
            }
        }

        /// <summary>
        /// Gets the types required by the estimator in order to work:
        /// InterestRateMarketData is the only required type for this estimator.
        /// </summary>
        /// <param name="settings">The parameter is not used.</param>
        /// <param name="multivariateRequest">The parameter is not used.</param>
        /// <returns>An array containing the type InterestRateMarketData.</returns>
        public virtual EstimateRequirement[] GetRequirements(IEstimationSettings settings, EstimateQuery query)
        {
            return new EstimateRequirement[] { new EstimateRequirement(typeof(InterestRateMarketData)) };
        }

        /// <summary>
        /// Attempts a calibration through <see cref="PelsserCappletOptimizationProblem"/>
        /// using caps matrices.
        /// </summary>
        /// <param name="data">The data to be used in order to perform the calibration.</param>
        /// <param name="settings">The parameter is not used.</param>
        /// <param name="controller">The controller which may be used to cancel the process.</param>
        /// <returns>The results of the calibration.</returns>
        public EstimationResult Estimate(List<object> data, IEstimationSettings settings = null, IController controller = null, Dictionary<string, object> properties = null)
        {
            InterestRateMarketData dataset = data[0] as InterestRateMarketData;
            MatrixMarketData normalVol = null;
            if (data.Count > 1)
                normalVol = (MatrixMarketData)data[1];

            // initialize the result
            EstimationResult result;

            // names of the parameters
            string[] names = new string[] { "alpha1", "sigma1" };

            // we just need ZR Market data for dummy calibration 
            if ( (dataset.ZRMarket != null) && false)
            {
                var dummyValues = new double[]{ 0.014, 0.001 };
                result = new EstimationResult(names, dummyValues);

                result.ZRX = (double[])dataset.ZRMarketDates.ToArray();
                result.ZRY = (double[])dataset.ZRMarket.ToArray();
                result.Objects = new object[1];
                result.Objects[0] = 0.0;

                return result;
            }


            if ((dataset.ZRMarket == null) || (dataset.CapVolatility == null))
            {
                result = new EstimationResult();
                result.ErrorMessage = "Not enough data to calibrate.\n" +
                    "The estimator needs a ZRMarket and a CapVolatility " +
                    "defined inside InterestRateMarketData";
                return result;
            }

            // Backup the dates
            DateTime effectiveDate = DateTime.Now.Date;
            DateTime valuationDate = DateTime.Now.Date;
            if (Document.ActiveDocument != null)
            {
                effectiveDate = Document.ActiveDocument.ContractDate;
                valuationDate = Document.ActiveDocument.SimulationStartDate;
            }

            // Creates the Context.
            Document doc = new Document();
            ProjectROV prj = new ProjectROV(doc);
            doc.Part.Add(prj);

            Function zr = new PFunction(null);
            zr.VarName = "zr";
            // Load the zr.
            double[,] zrvalue = (double[,])ArrayHelper.Concat(dataset.ZRMarketDates.ToArray(), dataset.ZRMarket.ToArray());
            zr.Expr = zrvalue;
            
            prj.Symbols.Add(zr);


            var bm = BlackModelFactory(zr);
            if (bm is BachelierNormalModel)
            {
                (bm as BachelierNormalModel).Tenor = dataset.CapTenor;
            }

            double deltak = dataset.CapTenor;


            Matrix capVol = normalVol != null ? normalVol.Values:dataset.CapVolatility;
            Vector capMat = normalVol != null ? normalVol.RowValues: dataset.CapMaturity;
            Vector capK = normalVol != null ? GetStrikes(normalVol,zr): dataset.CapRate;

            var preferences = settings as Fairmat.Calibration.CapVolatilityFiltering;

            // Matrix calculated with black.
            var blackCaps = new Matrix(capMat.Length, capK.Length);
            for (int m = 0; m < capMat.Length; m++)
            {
                for (int s = 0; s < capK.Length; s++)
                {
                    bool skip = false;
                    if (preferences != null)
                    {
                        if (capK[s] < preferences.MinCapRate || capK[s] > preferences.MaxCapRate ||
                           capMat[m] < preferences.MinCapMaturity || capMat[m] > preferences.MaxCapMaturity)
                                {skip = true; }
                    }

                    var adjS = s;
                    if (IsAtmMatrix(normalVol)) { adjS = 0; }
                    if (capVol[m, adjS] == 0 || skip)
                        blackCaps[m, s] = 0;
                    else
                        blackCaps[m, s] = bm.Cap(capK[s], capVol[m, adjS], deltak, capMat[m]);
                }
            }

            if (blackCaps.IsNaN())
            {
                Console.WriteLine("Black caps matrix has non real values:");
                Console.WriteLine(blackCaps);
                throw new Exception("Cannot calculate Black caps");
            }

            // Maturity goes from 0 to the last item with step deltaK.
            Vector maturity = new Vector((int)(1.0 + capMat[capMat.Length - 1] / deltak));
            for (int l = 0; l < maturity.Length; l++)
                maturity[l] = deltak * l;

            Vector fwd = new Vector(maturity.Length - 1);
            for (int i = 0; i < fwd.Length; i++)
            {
                fwd[i] = bm.Fk(maturity[i + 1], deltak);
            }

            // Creates a default Pelsser model.
            Pelsser.SquaredGaussianModel model = new Pelsser.SquaredGaussianModel();
            model.a1 = (ModelParameter)0.014;
            model.sigma1 = (ModelParameter)0.001;
            model.zr = (ModelParameter)"@zr";
            StochasticProcessExtendible iex = new StochasticProcessExtendible(prj, model);
            prj.Processes.AddProcess(iex);

            prj.Parse();

            DateTime t0 = DateTime.Now;
            Caplet cp = new Caplet();

            PelsserCappletOptimizationProblem problem = new PelsserCappletOptimizationProblem(prj, cp, maturity, fwd, capK, deltak, capMat, blackCaps);

            IOptimizationAlgorithm solver = new QADE();
            IOptimizationAlgorithm solver2 = new SteepestDescent();

            DESettings o = new DESettings();
            o.NP = 35;
            o.TargetCost = 0.0025;
            o.MaxIter = 10;
            o.Verbosity = Math.Max(1, Engine.Verbose);
            o.controller = controller;
            // Parallel evaluation is not supported for this calibration.
            o.Parallel = false;
            o.Debug = true;
            SolutionInfo solution = null;

            Vector x0 = (Vector)new double[] { 0.1, 0.1 };

            solution = solver.Minimize(problem, o, x0);
            if (solution.errors)
                return new EstimationResult(solution.message);

            o.epsilon = 10e-7;
            o.h = 10e-7;
            o.MaxIter = 1000;
            o.Debug = true;
            o.Verbosity = Math.Max(1, Engine.Verbose);

            if (solution != null)
                solution = solver2.Minimize(problem, o, solution.x);
            else
                solution = solver2.Minimize(problem, o, x0);
            
            if (solution.errors)
                return new EstimationResult(solution.message);

            Console.WriteLine(solution);

            result = new EstimationResult(names, solution.x);

            result.ZRX = (double[])dataset.ZRMarketDates.ToArray();
            result.ZRY = (double[])dataset.ZRMarket.ToArray();
            result.Objects = new object[1];
            result.Objects[0] = solution.obj;
            //result.Fit = solution.obj;//Uncomment in 1.6
            // Restore the dates
            if (Document.ActiveDocument != null)
            {
                Document.ActiveDocument.ContractDate = effectiveDate;
                Document.ActiveDocument.SimulationStartDate = valuationDate;
            }

            return result;
        }

        /// <summary>
        /// Returns the Strike vector embedded in the volatility matrix. 
        /// The method manage the convention that if rates/strikes [defined in .ColumnValues fields] are equal to -1,
        /// the strike vector is ATM and must be derived from the term structure. 
        /// </summary>
        /// <param name="normalVol">The matrix</param>
        /// <param name="zr">The term structure.</param>
        /// <returns>The strike vector.</returns>
        static Vector GetStrikes(MatrixMarketData normalVol, IFunction zr)
        {
            if (IsAtmMatrix(normalVol))
            {
                // builds the vector of ATM strikes
                var s = new Vector(normalVol.RowValues.Length);
                for (int j = 0; j < s.Length; j++)
                    s[j] = zr.Evaluate(normalVol.RowValues[j]);
                return s;
            }
            else
                return normalVol.ColumnValues;
        }

        private static bool IsAtmMatrix(MatrixMarketData normalVol)
        {
            return normalVol != null &&
                   normalVol.ColumnValues.Length == 1 &&
                   normalVol.ColumnValues[0] == -1;
        }

        #endregion

        public virtual string Description
        {
            get { return "Calibrate against Caplet prices"; }
        }

        public IEstimationSettings DefaultSettings
        {
            get { return new Fairmat.Calibration.CapVolatilityFiltering(); }
        }

        protected virtual BlackModel BlackModelFactory(Function zr)
        {
            return new BlackModel(zr);
        }

        public virtual string ToolTipText
        {
            get { return "Caplet Estimator"; }
        }
    }
}
