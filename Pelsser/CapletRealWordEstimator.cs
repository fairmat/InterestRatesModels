using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Mono.Addins;
using DVPLI;
using Fairmat.MarketData;
using Fairmat.Finance;
using Pelsser;

namespace Pelsser.Calibration
{
    /// <summary>
    /// Implementation of HW1 real world calibration (swaption matrix and historical zr based).
    /// </summary>
    [Extension("/Fairmat/Estimator")]
    class CapletRealWordEstimator : IEstimatorEx, IMenuItemDescription
    {
        /// <summary>
        /// Gets the tooltip for the implemented calibration function.
        /// </summary>
        public string ToolTipText
        {
            get
            {
                return "Real world calibration using caps and zero rate curve historical series";
            }
        }

        /// <summary>
        /// Gets the description of the implemented calibration function.
        /// </summary>
        public string Description
        {
            get
            {
                return "Real world calibration using caps and zero rate curve historical series";
            }
        }


        /// <summary>
        /// Gets the value requested by the interface ProvidesTo,
        /// returning HW1 as the type.
        /// </summary>
        public virtual Type ProvidesTo
        {
            get
            {
                return typeof(SquaredGaussianModel);
            }
        }

        /// <summary>
        /// Gets the types required by the estimator in order to work:
        /// InterestRateMarketData and zero rate curve historical serie are required for this estimator.
        /// </summary>
        /// <param name="settings">The parameter is not used.</param>
        /// <param name="multivariateRequest">The parameter is not used.</param>
        /// <returns>An array containing the type InterestRateMarketData.</returns>
        public EstimateRequirement[] GetRequirements(IEstimationSettings settings, EstimateQuery query)
        {
            return new EstimateRequirement[] { new EstimateRequirement(typeof(InterestRateMarketData)),
                new EstimateRequirement(typeof(DiscountingCurveMarketData[]))};
        }

        /// <summary>
        /// Attempts a calibration through <see cref="SwaptionHW1OptimizationProblem"/>
        /// using swaption matrices.
        /// </summary>
        /// <param name="data">The data to be used in order to perform the calibration.</param>
        /// <param name="settings">The parameter is not used.</param>
        /// <param name="controller">The controller which may be used to cancel the process.</param>
        /// <returns>The results of the calibration.</returns>
        public EstimationResult Estimate(List<object> data, IEstimationSettings settings = null, IController controller = null, Dictionary<string, object> properties = null)
        {
            LambdaCalibrationSettings lsettings = (LambdaCalibrationSettings)settings;
            if (lsettings.Years > lsettings.BondMaturity)
                throw new Exception("Bond maturity has to be greater of the historical series time span.");

            InterestRateMarketData irmd = data[0] as InterestRateMarketData;
            List<object> IrmdData = new List<object>();
            IrmdData.Add(irmd);
            CapletEstimator CapEst = new CapletEstimator();
            EstimationResult er1 = CapEst.Estimate(IrmdData);

            DiscountingCurveMarketData[] dcmd = Array.ConvertAll<IMarketData, DiscountingCurveMarketData>
                ((IMarketData[])data[1], el => (DiscountingCurveMarketData)el);
            MarketPriceOfRiskCalculator mporc = new MarketPriceOfRiskCalculator();
            double lambda = mporc.Calculate(dcmd, lsettings);

            string[] names = new string[er1.Names.Length + 1];
            for (int i = 0; i < er1.Names.Length; i++)
                names[i] = er1.Names[i];
            names[er1.Names.Length] = "lambda0";
            Vector values = new Vector(er1.Values.Length + 1);
            values[new Range(0, values.Length - 2)] = (Vector)er1.Values;
            values[Range.End] = lambda;
            EstimationResult result = new EstimationResult(names, values);
            return result;
        }

        public IEstimationSettings DefaultSettings
        {
            get
            {
                return new LambdaCalibrationSettings();
            }
        }
    }
}
