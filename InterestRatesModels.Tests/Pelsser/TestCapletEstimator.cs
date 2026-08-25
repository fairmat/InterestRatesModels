using System;
using System.Collections.Generic;
using System.Reflection;
using DVPLDOM;
using DVPLI;
using Fairmat.MarketData;
using NUnit.Framework;

namespace Pelsser.Calibration
{
    [TestFixture]
    public class TestCapletEstimator
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        private static object InvokePrivateStatic(string methodName, params object[] args)
        {
            MethodInfo method = typeof(CapletEstimator).GetMethod(methodName,
                BindingFlags.NonPublic | BindingFlags.Static);
            return method.Invoke(null, args);
        }

        private static InterestRateMarketData CreateMarketData()
        {
            return new InterestRateMarketData
            {
                ZRMarketDates = new Vector(new double[] { 0, 1, 2, 5, 10 }),
                ZRMarket = new Vector(new double[] { 0.01, 0.015, 0.017, 0.02, 0.025 }),
                CapMaturity = new Vector(new double[] { 1, 2 }),
                CapRate = new Vector(new double[] { 0.01, 0.02 }),
                CapTenor = 0.5,
                CapVolatility = new Matrix(new double[,] { { 0.20, 0.22 }, { 0.21, 0.23 } })
            };
        }

        [Test]
        public void IsAtmMatrixIsTrueOnlyForSingleColumnValueMinusOne()
        {
            MatrixMarketData atm = new MatrixMarketData { ColumnValues = new Vector(new double[] { -1 }) };
            MatrixMarketData notAtm = new MatrixMarketData { ColumnValues = new Vector(new double[] { 0.01, 0.02 }) };
            MatrixMarketData wrongValue = new MatrixMarketData { ColumnValues = new Vector(new double[] { 0.01 }) };

            Assert.IsTrue((bool)InvokePrivateStatic("IsAtmMatrix", atm));
            Assert.IsFalse((bool)InvokePrivateStatic("IsAtmMatrix", notAtm));
            Assert.IsFalse((bool)InvokePrivateStatic("IsAtmMatrix", wrongValue));
            Assert.IsFalse((bool)InvokePrivateStatic("IsAtmMatrix", (MatrixMarketData)null));
        }

        [Test]
        public void GetStrikesForAtmMatrixEvaluatesZeroCurveAtEachRowDate()
        {
            MatrixMarketData atm = new MatrixMarketData
            {
                ColumnValues = new Vector(new double[] { -1 }),
                RowValues = new Vector(new double[] { 1, 2, 5 })
            };
            Function zr = new PFunction(null);
            zr.Expr = new double[,] { { 0, 0.01 }, { 10, 0.03 } };

            Vector strikes = (Vector)InvokePrivateStatic("GetStrikes", atm, zr);

            Assert.AreEqual(3, strikes.Length);
            for (int i = 0; i < strikes.Length; i++)
                Assert.AreEqual(zr.Evaluate(atm.RowValues[i]), strikes[i], 1e-12);
        }

        [Test]
        public void GetStrikesForNonAtmMatrixReturnsColumnValues()
        {
            MatrixMarketData nonAtm = new MatrixMarketData
            {
                ColumnValues = new Vector(new double[] { 0.01, 0.02 }),
                RowValues = new Vector(new double[] { 1, 2 })
            };
            Function zr = new PFunction(null);
            zr.Expr = new double[,] { { 0, 0.01 }, { 10, 0.03 } };

            Vector strikes = (Vector)InvokePrivateStatic("GetStrikes", nonAtm, zr);

            Assert.That(strikes, Is.SameAs(nonAtm.ColumnValues));
        }

        [Test]
        public void EstimateReturnsErrorWhenZRMarketIsMissing()
        {
            CapletEstimator estimator = new CapletEstimator();
            InterestRateMarketData dataset = CreateMarketData();
            dataset.ZRMarket = null;

            EstimationResult result = estimator.Estimate(new List<object> { dataset });

            Assert.AreEqual(
                "Not enough data to calibrate.\n" +
                "The estimator needs a ZRMarket and a CapVolatility " +
                "defined inside InterestRateMarketData",
                result.ErrorMessage);
        }

        [Test]
        public void EstimateReturnsErrorWhenCapVolatilityIsMissing()
        {
            CapletEstimator estimator = new CapletEstimator();
            InterestRateMarketData dataset = CreateMarketData();
            dataset.CapVolatility = null;

            EstimationResult result = estimator.Estimate(new List<object> { dataset });

            Assert.AreEqual(
                "Not enough data to calibrate.\n" +
                "The estimator needs a ZRMarket and a CapVolatility " +
                "defined inside InterestRateMarketData",
                result.ErrorMessage);
        }

        [Test]
        public void EstimateDummyCalibrationReturnsFixedGuess()
        {
            CapletEstimator estimator = new CapletEstimator();
            InterestRateMarketData dataset = CreateMarketData();
            Fairmat.Calibration.CapVolatilityFiltering settings = new Fairmat.Calibration.CapVolatilityFiltering
            {
                DummyCalibration = true
            };

            EstimationResult result = estimator.Estimate(new List<object> { dataset }, settings);

            CollectionAssert.AreEqual(new string[] { "alpha1", "sigma1" }, result.Names);
            CollectionAssert.AreEqual(new double[] { 0.014, 0.001 }, result.Values);
            CollectionAssert.AreEqual((double[])dataset.ZRMarketDates.ToArray(), result.ZRX);
            CollectionAssert.AreEqual((double[])dataset.ZRMarket.ToArray(), result.ZRY);
            Assert.AreEqual(1, result.Objects.Length);
            Assert.AreEqual(0.0, result.Objects[0]);
        }

        [Test]
        public void EstimateDummyCalibrationIsSkippedWhenZRMarketIsMissing()
        {
            CapletEstimator estimator = new CapletEstimator();
            InterestRateMarketData dataset = CreateMarketData();
            dataset.ZRMarket = null;
            Fairmat.Calibration.CapVolatilityFiltering settings = new Fairmat.Calibration.CapVolatilityFiltering
            {
                DummyCalibration = true
            };

            // Falls through to the validation branch instead of taking the dummy shortcut.
            EstimationResult result = estimator.Estimate(new List<object> { dataset }, settings);

            Assert.IsFalse(string.IsNullOrEmpty(result.ErrorMessage));
        }

        [Test]
        public void EstimateThrowsWhenSecondDataItemIsNotAMatrixMarketData()
        {
            CapletEstimator estimator = new CapletEstimator();
            InterestRateMarketData dataset = CreateMarketData();

            Assert.Throws<InvalidCastException>(() =>
                estimator.Estimate(new List<object> { dataset, "not a matrix market data" }));
        }
    }
}
