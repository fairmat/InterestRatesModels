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

namespace Pelsser.Calibration
{

    /// <summary>
    /// Extend caplet estimator using normal instead of log normal vol
    /// </summary>
    [Extension("/Fairmat/Estimator")]
    public class CappletEstimatorNormal : CapletEstimator
    {

        public override string Description
        {
            get { return "Calibrate against Caplet prices (Normal Volatilities)"; }
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
            new EstimateRequirement() { Field = "Vol-ATM-N", MarketDataType = typeof(MatrixMarketData), TickerReplacement = GetCapletVolatilityTicker(query.Market)} };
        }

        protected override BlackModel BlackModelFactory(Function zr)
        {
            return new BachelierNormalModel(zr);
        }

        static string GetCapletVolatilityTicker(string marketExpression)
        {
            string market = MarketDataRequestManager.GetMarket(marketExpression);
            string ticker = $"CapVol{CurrencyTickerBuilder.MarketToCurrency(market)}";
            return ticker;
        }

        public override string ToolTipText
        {
            get { return "Caplet Estimator (Normal Volatilities)"; }
        }
    }

}
