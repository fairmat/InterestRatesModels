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
using DVPLDOM;
using DVPLI;
using DVPLSolver;
using NUnit.Framework;

namespace Pelsser
{
    [TestFixture]
    public class TestPelsser
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void Test()
        {
            Engine.MultiThread = true;

            Document doc = new Document();
            ProjectROV rov = new ProjectROV(doc);
            doc.Part.Add(rov);

            AFunction zerorate = new AFunction(rov);
            zerorate.VarName = "zr";
            zerorate.m_IndependentVariables = 1;
            zerorate.m_Value = (RightValue)0.05;

            rov.Symbols.Add(zerorate);

            int n_sim = 5000;
            int n_steps = 900;
            SquaredGaussianModel process = new SquaredGaussianModel();
            process.a1 = (ModelParameter)0.1;
            process.sigma1 = (ModelParameter)0.01;
            process.zr = (ModelParameter)"@zr";

            StochasticProcessExtendible s = new StochasticProcessExtendible(rov, process);
            rov.Processes.AddProcess(s);

            // Set the discounting.
            RiskFreeInfo rfi = rov.GetDiscountingModel() as RiskFreeInfo;
            rfi.ActualizationType = EActualizationType.RiskFree;
            rfi.m_deterministicRF = 0.0;

            OptionTree op = new OptionTree(rov);
            op.PayoffInfo.PayoffExpression = "bond(t;10;@v1)";

            // Set the simulation maturity.
            op.PayoffInfo.Timing.EndingTime.m_Value = (RightValue)2.0;
            op.PayoffInfo.European = true;
            rov.Map.Root = op;

            rov.NMethods.Technology = ETechType.T_SIMULATION;
            rov.NMethods.PathsNumber = n_sim;
            rov.NMethods.SimulationSteps = n_steps;

            ROVSolver solver = new ROVSolver();
            solver.BindToProject(rov);
            solver.DoValuation(-1);

            if (rov.HasErrors)
            {
                Console.WriteLine(rov.m_RuntimeErrorList[0]);
            }

            Assert.IsFalse(rov.HasErrors);

            ResultItem price = rov.m_ResultList[0] as ResultItem;
            Console.WriteLine("Bond Test Value = " + price.value.ToString());

            Assert.LessOrEqual(Math.Abs(0.6702 - price.value), .01);

            // Try to do some simple tests and check the results.
            double b0_10 = process.Bond(DynamicParam(0, process), process.CacheDates, 0, 0, 10);
            Console.WriteLine("Bond(0,10) = " + b0_10);

            Assert.LessOrEqual(Math.Abs(b0_10 - 0.606513), .001);

            double b7_10 = process.Bond(DynamicParam(0.00427631, process), process.CacheDates, 0, 7, 10);
            Console.WriteLine("Bond(7,10) = " + b7_10);

            Assert.LessOrEqual(Math.Abs(b7_10 - 0.856374), .001);

            double b7_30 = process.Bond(DynamicParam(0.00427631, process), process.CacheDates, 0, 7, 30);
        }

        /// <summary>
        /// Format the y parameter so it can be made compatible with <see cref="Pelsser.Bond"/>.
        /// </summary>
        /// <param name="y">The y parameter wanted inside Bond.</param>
        /// <param name="process">The Pelsser process which will be used.</param>
        /// <returns>The Matrix to pass to Bond.</returns>
        private Matrix DynamicParam(double y, SquaredGaussianModel process)
        {
            double alphaT = process.F(process.CacheDates[0], process.CacheDates[1] - process.CacheDates[0]) +
                                2 * Math.Exp(-process.a1.V() * process.CacheDates[0]) * process.Int(0, process.CacheDates[0]);
            return new Matrix(new double[] { Math.Pow(y + alphaT, 2) });
        }
    }

    [TestFixture]
    public class TestPelsser2
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void Test()
        {
            // Comparison with the benchmark of the expected value of P(1,2)
            double mcPrice;
            double mcDevST;
            PelsserBond.Calculate(1.0, 2.0, 2.0, 100000, 200, out mcPrice, out mcDevST);

            // The benchmark value with 35000 paths
            // double theoreticalValue = 0.95121;
            double RightValue = 0.951408282184948;

            Console.WriteLine("Benchmark  Bond = " + RightValue);
            Console.WriteLine("Fairmat Bond = " + mcPrice);
            Console.WriteLine("Standard Deviation = " + mcDevST);
            double tol = 5.0 * mcDevST;
            Assert.Less(Math.Abs(RightValue - mcPrice), tol);
        }
    }

    [TestFixture]
    public class TestPelsser3
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void Test()
        {
            // Comparison with the benchmark of the expected value of  P(5,6)
            double mcPrice;
            double mcDevST;
            PelsserBond.Calculate(5.0, 6.0, 6.0, 100000, 200, out mcPrice, out mcDevST);

            // The benchmark value with 35000 paths
            // double theoreticalValue = 0.95109;
            double RightValue = 0.951276263137499;

            Console.WriteLine("Benchmark  Bond = " + RightValue);
            Console.WriteLine("Fairmat Bond = " + mcPrice);
            Console.WriteLine("Standard Deviation = " + mcDevST);
            double tol = 5.0 * mcDevST;
            Assert.Less(Math.Abs(RightValue - mcPrice), tol);
        }
    }

    [TestFixture]
    public class TestPelsser4
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test, Category("BigTest")]
        public void Test()
        {
            // Try several different simulation dates with the same bound.
            double mcDevST;
            int numTries = 21;
            double min, max, val;
            min = 0.0;
            max = 2.0;
            Vector timeVector = Vector.Linspace(min, max, numTries);
            Vector priceVector = new Vector(numTries);
            for (int j = 0; j < numTries; j++)
            {
                PelsserBond.Calculate(0.0, 1.0, timeVector[j], 10000, 100, out val, out mcDevST);
                priceVector[j] = val;
            }

            // The benchmark value with 35000 paths
            double theoreticalValue = 0.95121;
            Console.WriteLine("Benchmark  Bond = " + theoreticalValue);
            Console.WriteLine("Fairmat Bond = ");
            for (int j = 0; j < numTries; j++)
                Console.WriteLine("Simulation Time = " + timeVector[j] + "\t Bond Price = " + priceVector[j]);
        }
    }

    /// <summary>
    /// This test simulates a bound at several dates, the time horizon is the one bond expiration.
    /// </summary>
    [TestFixture]
    public class TestPelsser5
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }
        
        [Test]
        [Ignore("disabled unsupported test")]
        public virtual void Test()
        {
            Assert.IsTrue(Run(1.0));
        }

        protected bool Run(double h)
        {
            // Try several different simulation dates with the same bond.
            double mcDevST;
            double val;
            PelsserBond.Calculate(0.0, 1.0, h, 10000, 100, out val, out mcDevST);

            double theoreticalValue = 0.951228651454765;

            double diff = theoreticalValue - val;
            Console.WriteLine("Simulation horizon\t" + h);
            Console.WriteLine("Diff:\t" + (100 * diff / theoreticalValue) + "%");
            if (Math.Abs(diff) > 10e-6)
                return false;
            return true;
        }
    }

    /// <summary>
    /// Simulation expiration less than 1.
    /// </summary>
    [TestFixture]
    public class TestPelsser6 : TestPelsser5
    {
        [Test]
        [Ignore("disabled unsupported test")]
        public override void Test()
        {
            Assert.IsTrue(Run(.2));
        }
    }

    /// <summary>
    /// Simulation expiration greater than 1.
    /// </summary>
    [TestFixture]
    public class TestPelsser7 : TestPelsser5
    {
        [Test]
        [Ignore("disabled unsupported test")]
        public override void Test()
        {
            Assert.IsTrue(Run(7.5));
        }
    }

    public class PelsserBond
    {
        public static void Calculate(double t1, double t2, double simEnd, int numSim, int numSteps, out double val, out double stDev)
        {
            Engine.MultiThread = true;

            Document doc = new Document();
            ProjectROV rov = new ProjectROV(doc);
            doc.Part.Add(rov);

            AFunction zerorate = new AFunction(rov);
            zerorate.VarName = "zr";
            zerorate.m_IndependentVariables = 1;
            zerorate.m_Value = (RightValue)0.05;

            rov.Symbols.Add(zerorate);

            // To be changed to 350000.
            int n_sim = numSim;
            int n_steps = numSteps;
            SquaredGaussianModel process = new SquaredGaussianModel();
            process.a1 = (ModelParameter)0.1;
            process.sigma1 = (ModelParameter)0.01;
            process.zr = (ModelParameter)"@zr";

            StochasticProcessExtendible s = new StochasticProcessExtendible(rov, process);
            rov.Processes.AddProcess(s);

            // Set the discounting.
            RiskFreeInfo rfi = rov.GetDiscountingModel() as RiskFreeInfo;
            rfi.ActualizationType = EActualizationType.RiskFree;
            rfi.m_deterministicRF = 0.0;

            OptionTree op = new OptionTree(rov);
            op.PayoffInfo.PayoffExpression = "bond(" + t1.ToString() + ";" + t2.ToString() + ";@v1)";

            // Here we put the simulation maturity.
            op.PayoffInfo.Timing.EndingTime.m_Value = (RightValue)simEnd;
            op.PayoffInfo.European = true;
            rov.Map.Root = op;

            rov.NMethods.Technology = ETechType.T_SIMULATION;
            rov.NMethods.PathsNumber = n_sim;
            rov.NMethods.SimulationSteps = n_steps;

            ROVSolver solver = new ROVSolver();
            solver.BindToProject(rov);
            solver.DoValuation(-1);

            if (rov.HasErrors)
            {
                Console.WriteLine(rov.m_RuntimeErrorList[0]);
            }

            ResultItem price = rov.m_ResultList[0];
            val = price.value;
            stDev = price.stdDev / Math.Sqrt(numSim);
        }
    }
}
