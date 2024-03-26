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
using Fairmat.Optimization;
using Mono.Addins;

namespace CIRProcess
{
    /// <summary>
    /// Implements the estimator for CIR Cap calibration.
    /// </summary>
    [Extension("/Fairmat/Estimator")]
    public class CapCIREstimator : IEstimator, IMenuItemDescription
    {
        /// <summary>
        /// Gets the tooltip for the implemented calibration function.
        /// </summary>
        public string ToolTipText
        {
            get
            {
                return "Calibrate using Caps";
            }
        }

        /// <summary>
        /// Gets the description of the implemented calibration function.
        /// </summary>
        public string Description
        {
            get
            {
                return "Calibrate from cap matrix";
            }
        }

        #region IEstimator Members

        /// <summary>
        /// Gets the value requested by the interface ProvidesTo,
        /// returning <see cref="CIR"/> as the type.
        /// </summary>
        public virtual Type ProvidesTo
        {
            get
            {
                return typeof(CIR);
            }
        }

        /// <summary>
        /// Gets the types required by the estimator in order to work:
        /// InterestRateMarketData is the only required type for this estimator.
        /// </summary>
        /// <param name="settings">The parameter is not used.</param>
        /// <param name="query">The parameter is not used.</param>
        /// <returns>An array containing the type InterestRateMarketData.</returns>
        public EstimateRequirement[] GetRequirements(IEstimationSettings settings, EstimateQuery query)
        {
            return new EstimateRequirement[] { new EstimateRequirement(typeof(InterestRateMarketData)) };
        }

        /// <summary>
        /// Attempts a calibration through <see cref="CapsCIROptimizationProblem"/>
        /// using caps matrices.
        /// </summary>
        /// <param name="data">The data to be used in order to perform the calibration.</param>
        /// <param name="settings">The parameter is not used.</param>
        /// <param name="controller">A controller used for the optimization process.</param>
        /// <returns>The results of the calibration.</returns>
        public EstimationResult Estimate(List<object> data, IEstimationSettings settings = null, IController controller = null, Dictionary<string, object> properties = null)
        {
            InterestRateMarketData dataset = data[0] as InterestRateMarketData;

            // intialize the result
            EstimationResult result; 

            

            // Creates the context.
            Document doc = new Document();
            ProjectROV prj = new ProjectROV(doc);
            doc.Part.Add(prj);

            CapCIROptimizationProblem problem = new CapCIROptimizationProblem(dataset);

            // initial guess
            Vector x0 = new Vector(new double[] { 1, 0.01, 0.05 });

            if (settings != null && settings.DummyCalibration)
            {
                Vector valuesDummy = new Vector(4);
                valuesDummy[Range.New(0, 2)] = x0;
                valuesDummy[3] = problem.r0;

                result = new EstimationResult(CIR.parameterNames, valuesDummy);

                return result;
            }


            IOptimizationAlgorithm solver = new QADE();
            IOptimizationAlgorithm solver2 = new SteepestDescent();

            DESettings o = new DESettings();
            o.NP = 50;
            o.MaxIter = 50;
            o.Verbosity = 1;
            o.Parallel = false;
            o.controller = controller;

            int seed = AttributesUtility.RetrieveAttributeOrDefaultValue(properties, "Seed", -1);
            if (seed != -1)
            {
                o.Repeatable = true;
                o.RandomSeed = seed; 
            }
               

            SolutionInfo solution = null;

            solution = solver.Minimize(problem, o, x0);
            if (solution.errors)
                return new EstimationResult(solution.message);

            o.epsilon = 10e-10;
            o.h = 10e-10;
            o.MaxIter = 1000;

            if (solution != null)
                solution = solver2.Minimize(problem, o, solution.x);
            else
                solution = solver2.Minimize(problem, o, x0);

            if (solution.errors)
                return new EstimationResult(solution.message);

            Console.WriteLine("Solution:");
            Console.WriteLine(solution);
            string[] names = CIR.parameterNames;
            Vector values = new Vector(4);
            values[Range.New(0, 2)] = solution.x;
            values[3] = problem.r0;

            result = new EstimationResult(names, values);

            return result;
        }

        #endregion
    }
}
