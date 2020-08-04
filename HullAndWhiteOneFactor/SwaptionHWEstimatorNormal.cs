using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using DVPLDOM;
using DVPLI;
using Fairmat.Finance;
using Fairmat.Optimization;
using Fairmat.Calibration;
using Mono.Addins;
using Fairmat.MarketData;

namespace HullAndWhiteOneFactor
{

    /// <summary>
    /// Extend swaption estimator and uses normal instead of log normal vol
    /// </summary>
    [Extension("/Fairmat/Estimator")]
    public class SwaptionHWEstimatorNormal: SwaptionHWEstimator
    {

        public override string Description
        {
            get
            {
                return "Calibrate HW1 using Swaption [Normal]";
            }
        }

        /// <summary>
        /// List required market data inlcuding normal vol
        /// </summary>
        /// <param name="settings"></param>
        /// <param name="query"></param>
        /// <returns></returns>
        public override EstimateRequirement[] GetRequirements(IEstimationSettings settings, EstimateQuery query)
        {
            return new EstimateRequirement[] {
            new EstimateRequirement(typeof(InterestRateMarketData)),
            new EstimateRequirement() { Field = "Vol-ATM-N", MarketDataType = typeof(MatrixMarketData), TickerReplacement = GetSwaptionVolatilityTicker(query.Market)} };
        }

        protected override BlackModel BlackModelFactory(Function zr)
        {
            return new BachelierNormalModel(zr);
        }

        static string GetSwaptionVolatilityTicker(string marketExpression)
        {
            string market = MarketDataRequestManager.GetMarket(marketExpression);
            string ticker = $"SwaptionVol{CurrencyTickerBuilder.MarketToCurrency(market)}";
            return ticker;
        }
    }


    /// <summary>
    /// Wrapper which offers normal vol functionality to legacy version of HW
    /// </summary>
    [Extension("/Fairmat/Estimator")]
    public class SwaptionHWEstimatorNormalLegacy : SwaptionHWEstimatorNormal
    {
        public override Type ProvidesTo
        {
            get
            {
                return typeof(StocasticProcessHW);
            }
        }
    }
}
