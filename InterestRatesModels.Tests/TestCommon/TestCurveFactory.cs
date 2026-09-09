using DVPLDOM;
using DVPLUtils;

namespace TestCommon
{
    /// <summary>
    /// Shared factory methods for building <see cref="Function"/> zero rate curves
    /// used across multiple test fixtures.
    /// </summary>
    public static class TestCurveFactory
    {
        public static Function CreateFlatZeroCurve(double rate)
        {
            Function zeroratecurve = new PFunction(null);
            zeroratecurve.Expr = new double[,] { { 0, rate }, { 50, rate } };
            (zeroratecurve as PFunction).m_Function.iType = EInterpolationType.LINEAR;
            return zeroratecurve;
        }
    }
}
