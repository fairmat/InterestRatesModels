using DVPLI;

namespace TestCommon
{
    /// <summary>
    /// Shared factory methods for building <see cref="InterestRateMarketData"/> instances
    /// used across multiple test fixtures.
    /// </summary>
    public static class TestMarketDataFactory
    {
        public static InterestRateMarketData CreateCapMarketData(Vector capMaturity, Matrix capVolatility)
        {
            return new InterestRateMarketData
            {
                ZRMarketDates = new Vector(new double[] { 0, 1, 2, 5, 10 }),
                ZRMarket = new Vector(new double[] { 0.01, 0.015, 0.017, 0.02, 0.025 }),
                CapMaturity = capMaturity,
                CapRate = new Vector(new double[] { 0.01, 0.02 }),
                CapTenor = 0.5,
                CapVolatility = capVolatility
            };
        }
    }
}
