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
using Fairmat.Finance;
using Fairmat.Optimization;
using Mono.Addins;

namespace HullAndWhiteOneFactor
{
    /// <summary>
    /// Implementation of HW1 Calibration against caps matrices.
    /// </summary>
    [Extension("/Fairmat/Estimator")]
    public class CapsHWEstimator : IEstimator,IEstimatorEx, IMenuItemDescription
    {
        /// <summary>
        /// Gets the tooltip for the implemented calibration function.
        /// </summary>
        public string ToolTipText
        {
            get
            {
                return "Calibrate Using Caps";
            }
        }

        /// <summary>
        /// Gets the description of the implemented calibration function.
        /// </summary>
        public string Description
        {
            get
            {
                return "Calibrate HW1 using Caps";
            }
        }

        #region IEstimator Members

        /// <summary>
        /// Gets the value requested by the interface ProvidesTo,
        /// returning HW1 as the type.
        /// </summary>
        public virtual Type ProvidesTo
        {
            get
            {
                return typeof(HW1);
            }
        }

        /// <summary>
        /// Gets the types required by the estimator in order to work:
        /// InterestRateMarketData is the only required type for this estimator.
        /// </summary>
        /// <param name="settings">The parameter is not used.</param>
        /// <param name="multivariateRequest">The parameter is not used.</param>
        /// <returns>An array containing the type InterestRateMarketData.</returns>
        public EstimateRequirement[] GetRequirements(IEstimationSettings settings, EstimateQuery query)
        {
            return new EstimateRequirement[] { new EstimateRequirement(typeof(InterestRateMarketData)) };
        }

        /*
        * IEstimatorEx.
        */
        public IEstimationSettings DefaultSettings
        {
            get { return new Fairmat.Calibration.CapVolatilityFiltering(); }
        }
       


        /// <summary>
        /// Attempts a calibration through <see cref="CapsHW1OptimizationProblem"/>
        /// using caps matrices.
        /// </summary>
        /// <param name="data">The data to be used in order to perform the calibration.</param>
        /// <param name="settings">The parameter is not used.</param>
        /// <param name="controller">The controller which may be used to cancel the process.</param>
        /// <returns>The results of the calibration.</returns>
        public EstimationResult Estimate(List<object> data, IEstimationSettings settings = null, IController controller = null, Dictionary<string, object> properties = null)
        {
            InterestRateMarketData dataset = data[0] as InterestRateMarketData;

            // names of the parameters
            string[] names = { "Alpha", "Sigma" };

            // initialize the result 
            EstimationResult result;

            


            PFunction zr = new PFunction(null);
            zr.VarName = "zr";

            var preferences = settings as Fairmat.Calibration.CapVolatilityFiltering;

            // Loads ZR
            double[,] zrvalue = (double[,])ArrayHelper.Concat(dataset.ZRMarketDates.ToArray(), dataset.ZRMarket.ToArray());
            zr.Expr = zrvalue;

            BlackModel bm = new BlackModel(zr);

            double deltak = dataset.CapTenor;

            if (dataset.CapVolatility == null)
                return new EstimationResult("Cap not available at requested date");


            Matrix capVolatility = dataset.CapVolatility;
            Vector capMaturity = dataset.CapMaturity;
            Vector capRate = dataset.CapRate;
            double a = 0.1;
            double sigma = 0.1;

            // Matrix calculated with Black.
            Matrix blackCaps = new Matrix(capMaturity.Length, capRate.Length);
            Matrix logic = new Matrix(capMaturity.Length, capRate.Length);
            
            for (int m = 0; m < capMaturity.Length; m++)
            {
                for (int s = 0; s < capRate.Length; s++)
                {
                    blackCaps[m, s] = bm.Cap(capRate[s], capVolatility[m, s], deltak, capMaturity[m]);
                    if (double.IsNaN(blackCaps[m, s]))
                    {
                        bm.Cap(capRate[s], capVolatility[m, s], deltak, capMaturity[m]);
                        throw new Exception("Malformed black caps");
                    }

                    if (blackCaps[m, s] == 0.0)
                    {
                        logic[m, s] = 0.0;
                    }
                    else
                    {
                        logic[m, s] = 1.0;
                    }

                    //filter
                    if (preferences != null)
                    {
                        if (capRate[s] < preferences.MinCapRate || capRate[s] > preferences.MaxCapRate ||
                            capMaturity[m]<preferences.MinCapMaturity|| capMaturity[m]>preferences.MaxCapMaturity)
                                {logic[m, s] = 0; blackCaps[m, s] = 0;}
                    }


                }
            }

            DateTime t0 = DateTime.Now;
            CapHW1 hw1Caps = new CapHW1(zr);
            Matrix caps = hw1Caps.HWMatrixCaps(capMaturity, capRate, a, sigma, deltak);

            for (int m = 0; m < capMaturity.Length; m++)
            {
                for (int s = 0; s < capRate.Length; s++)
                {
                    caps[m, s] = logic[m, s] * caps[m, s];
                }
            }

            CapsHW1OptimizationProblem problem = new CapsHW1OptimizationProblem(hw1Caps, blackCaps, capMaturity, capRate, deltak);
            Vector provaparam = new Vector(2);

            var solver = new QADE();

            IOptimizationAlgorithm solver2 = new SteepestDescent();

            DESettings o = new DESettings();
            o.NP = 20;
            o.MaxIter = 10;
            o.Verbosity = 1;
            o.Parallel = false;
            SolutionInfo solution = null;
            Vector x0 = new Vector(new double[] { 0.05, 0.01 });
            o.controller = controller;


            var seed = AttributesUtility.RetrieveAttributeOrDefaultValue(properties, "Seed", -1);
            if (seed != -1)
            {
                o.Repeatable = true;
                o.RandomSeed = seed;
            }

            solution = solver.Minimize(problem, o, x0);

            o.epsilon = 10e-8;
            o.h = 10e-8;
           
         
            o.MaxIter = 100;
            solution = solver2.Minimize(problem, o, solution.x);
            if (solution.errors)
                return new EstimationResult(solution.message);
            Console.WriteLine("Solution:");
            Console.WriteLine(solution);

            //solution.x[0] *= 3;

            result = new EstimationResult(names, solution.x);

            result.ZRX = (double[])dataset.ZRMarketDates.ToArray();
            result.ZRY = (double[])dataset.ZRMarket.ToArray();

            return result;
        }

        #endregion


       


    }
}
