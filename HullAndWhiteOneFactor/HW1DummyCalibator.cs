/* Copyright (C) 2017 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
 * Author(s):  Matteo Tesser (matteo.tesser@fairmat.com)
 *            
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
using System.Linq;
using DVPLDOM;
using DVPLI;
using Fairmat.Finance;
using Fairmat.Optimization;
using Fairmat.MarketData;
using Mono.Addins;

namespace HullAndWhiteOneFactor
{
    /// <summary>
    /// Implementation of HW1 Calibration using current ZR but dummy parameters for alpha and sigma
    /// </summary>
    [Extension("/Fairmat/Estimator")]
    public class HW1DummyCalibator : IEstimatorEx, IMenuItemDescription
    {
        public IEstimationSettings DefaultSettings
        {
            get { return null; }
        }


        public Type ProvidesTo
        {
            get { return typeof(HW1);}
        }
    

        public string ToolTipText
        {
            get { return "Provides default values of HW1 process"; }
        }
    

        public string Description
        {
            get { return "Provides default values of HW1 process"; }
        }

        public EstimationResult Estimate(List<object> data, IEstimationSettings settings = null, IController controller = null, Dictionary<string, object> properties = null)
        {
            string[] names = new string[] { "Alpha", "Sigma" };

            var result = new EstimationResult(names, new double[] { 0.1, 0.05 });
            var forwardingCurve = data[0] as DiscountingCurveMarketData;
            result.ZRX = (double[])forwardingCurve.Durations;
            result.ZRY = (double[])forwardingCurve.Values;
            return result;
        }

        public EstimateRequirement[] GetRequirements(IEstimationSettings settings, EstimateQuery query)
        {
            return new EstimateRequirement[] { 
                new EstimateRequirement(typeof(DiscountingCurveMarketData))};
        }

       


    }
}
