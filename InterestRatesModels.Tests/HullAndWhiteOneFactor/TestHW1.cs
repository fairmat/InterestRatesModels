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
using DVPLSolver;
using Fairmat.Statistics;
using Mono.Addins;
using NUnit.Framework;

namespace HullAndWhiteOneFactor
{
    /// <summary>
    /// Implementation of HW1 Calibration (swaption matrix based) Test.
    /// </summary>
    [TestFixture]
    public class TestHW1
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void Test()
        {
            ProjectROV r = HullAndWhite1("bond(t;2;@V1)", 1, 1, .05);
            r.Container.NMethods.Technology = ETechType.T_SIMULATION;
            r.Container.NMethods.m_UseRepeatableSequence = true;
            bool parse = r.Parse();
            Assert.IsTrue(parse);

            AnalysisValuation valuator = new AnalysisValuation();
            valuator.BindToProject(r);
            valuator.RunAnalysis(-1);

            if (r.HasErrors)
            {
                Console.WriteLine("Errors:");
                foreach (Exception ex in r.m_RuntimeErrorList)
                    Console.WriteLine(ex.Message);
            }

            Assert.IsFalse(r.HasErrors);

            double v = r.m_ResultList[0].m_Value;
            Console.WriteLine("v = " + v.ToString());

            Assert.Less(Math.Abs(v - 0.9136), 0.0001);
        }

        private static ProjectROV HullAndWhite1(string payoff, double maturity, double a1, double sigma1)
        {
            Document doc = new Document();
            ProjectROV rov1 = new ProjectROV(doc);
            doc.Part.Add(rov1);

            // Create the zero rate curve
            double a = 0.08;
            double b = 0.05;
            double c = -0.18;
            string zr = string.Format("{0} - {1}*exp({2}*x1)", a, b, c);

            AFunction zero_rate = new AFunction(rov1);
            zero_rate.VarName = "ZeroRate";
            zero_rate.m_IndependentVariables = 1;
            zero_rate.m_Value = new RightValueExpression(zr);

            // Add to the project the created zero rate curve.
            rov1.Symbols.Add(zero_rate);

            RiskFreeInfo rfi = rov1.GetDiscountingModel() as RiskFreeInfo;
            rfi.ActualizationType = EActualizationType.ZeroCoupond;
            rfi.m_deterministicRF.m_Value = (RightValue)"exp( -ZeroRate(t)*t)";

            // Create the short rate process.
            HW1 hw1 = new HW1(a1, sigma1, "@ZeroRate");

            StochasticProcessExtendible hw = new StochasticProcessExtendible(rov1, hw1);

            rov1.Processes.AddProcess(hw);

            OptionTree ot = new OptionTree(rov1);
            ot.European = true;
            ot.PayoffInfo.PayoffExpression = payoff;
            ot.PayoffInfo.Timing.EndingTime.m_Value = (RightValue)maturity;

            rov1.Map.Root = ot;

            return rov1;
        }
    }

    /// <summary>
    /// Tests HW1 dynamics calculating a caplet price through simulation
    /// and compares it with the price calculated through closed formula.
    /// </summary>
    [TestFixture]
    public class TestHW1Caplet
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
            doc.DefaultProject.NMethods.m_UseAntiteticPaths = true;

            int n_sim = 20000;
            int n_steps = 1024;
            double a = 0.2;
            double DR = 0.02;
            double r0 = 0.015;
            double a1 = 0.02;
            double sigma1 = 0.01;
            double maturityOpt = 5.0;
            double strike = 0.005;
            double tau = 0.5;
            double strike2 = 1.0 / (1.0 + strike * tau);

            ModelParameter PT = new ModelParameter(maturityOpt, "TT");
            PT.VarName = "TT";
            rov.Symbols.Add(PT);

            ModelParameter Ptau = new ModelParameter(tau, "tau");
            Ptau.VarName = "tau";
            rov.Symbols.Add(Ptau);

            ModelParameter Pa = new ModelParameter(a, "a");
            Pa.VarName = "a";
            rov.Symbols.Add(Pa);

            ModelParameter PDR = new ModelParameter(DR, "PDR");
            PDR.VarName = "DR";
            rov.Symbols.Add(PDR);

            ModelParameter Pr0 = new ModelParameter(r0, "r0");
            Pr0.VarName = "r0";
            rov.Symbols.Add(Pr0);

            ModelParameter Pstrike = new ModelParameter(strike, "strike");
            Pstrike.VarName = "strike";
            rov.Symbols.Add(Pstrike);

            AFunction zerorate = new AFunction(rov);
            zerorate.VarName = "zr";
            zerorate.m_IndependentVariables = 1;
            zerorate.m_Value = (RightValue)("(1-exp(-a*x1))*DR + r0");
            rov.Symbols.Add(zerorate);

            HW1 process = new HW1(a1, sigma1, "@zr");

            StochasticProcessExtendible s = new StochasticProcessExtendible(rov, process);
            rov.Processes.AddProcess(s);

            // Set the discounting.
            RiskFreeInfo rfi = rov.GetDiscountingModel() as RiskFreeInfo;
            rfi.ActualizationType = EActualizationType.Stochastic;
            rfi.m_deterministicRF = (ModelParameter)"@V1";

            OptionTree op = new OptionTree(rov);

            // 1) RATE FUNCTION, with this the price is higher than the theoretical one
            // op.PayoffInfo.PayoffExpression = "tau*max(rate(TT;tau;@v1) - strike; 0)";
            // 2) OBTAIN RATE FROM bond = exp(-rate*t),
            // with this the price is higher than the theoretical one but it's more near than 1)
            // op.PayoffInfo.PayoffExpression = "tau*max(-ln(bond(TT;TT+tau;@v1))/tau - strike; 0)";
            // 3) CONVERT RATE from discrete to continuous through (1+r_d) = exp(r_c)
            // In this way the price is the same as the theoretical one.
            op.PayoffInfo.PayoffExpression = "tau*max(ln(1+rate(TT;tau;@v1)) - strike; 0)";

            op.PayoffInfo.Timing.EndingTime.m_Value = (RightValue)maturityOpt;
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
            double samplePrice = price.m_Value;

            double sampleDevSt = price.m_StdErr / Math.Sqrt(2.0 * (double)n_sim);

            // Calculation of the theoretical value of the caplet.
            CapHW1 cap = new CapHW1(zerorate);
            double theoreticalPrice = cap.HWCaplet(a1, sigma1, maturityOpt,
                                                   maturityOpt + tau, strike2);
            Console.WriteLine("Theoretical Price = " + theoreticalPrice.ToString());
            Console.WriteLine("Monte Carlo Price = " + samplePrice);
            Console.WriteLine("Standard Deviation = " + sampleDevSt.ToString());
            double tol = 4.0 * sampleDevSt;

            Assert.Less(Math.Abs(theoreticalPrice - samplePrice), tol);
        }
    }

    [TestFixture]
    public class TestHW1BondCall
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void Test()
        {
            // Tests HW1 dynamics comparing the price of a call option on a bond
            // calculated through simulation and the theoretical one.
            Engine.MultiThread = true;
            Document doc = new Document();
            ProjectROV rov = new ProjectROV(doc);
            doc.Part.Add(rov);
            doc.DefaultProject.NMethods.m_UseAntiteticPaths = true;

            int n_sim = 20000;
            int n_steps = 1024;
            double a = 0.2;
            double DR = 0.02;
            double r0 = 0.015;
            double a1 = 0.02;
            double sigma1 = 0.01;
            double maturityOpt = 5.0;
            double strike = 0.98192;
            double tau = 1.0;

            ModelParameter PT = new ModelParameter(maturityOpt, "TT");
            PT.VarName = "TT";
            rov.Symbols.Add(PT);

            ModelParameter Ptau = new ModelParameter(tau, "tau");
            Ptau.VarName = "tau";
            rov.Symbols.Add(Ptau);

            ModelParameter Pa = new ModelParameter(a, "a");
            Pa.VarName = "a";
            rov.Symbols.Add(Pa);

            ModelParameter PDR = new ModelParameter(DR, "PDR");
            PDR.VarName = "DR";
            rov.Symbols.Add(PDR);

            ModelParameter Pr0 = new ModelParameter(r0, "r0");
            Pr0.VarName = "r0";
            rov.Symbols.Add(Pr0);

            ModelParameter Pstrike = new ModelParameter(strike, "strike");
            Pstrike.VarName = "strike";
            rov.Symbols.Add(Pstrike);

            AFunction zerorate = new AFunction(rov);
            zerorate.VarName = "zr";
            zerorate.m_IndependentVariables = 1;
            zerorate.m_Value = (RightValue)("(1-exp(-a*x1))*DR + r0");
            rov.Symbols.Add(zerorate);

            HW1 process = new HW1(a1, sigma1, "@zr");

            StochasticProcessExtendible s = new StochasticProcessExtendible(rov, process);
            rov.Processes.AddProcess(s);

            // Set the discounting.
            RiskFreeInfo rfi = rov.GetDiscountingModel() as RiskFreeInfo;
            rfi.ActualizationType = EActualizationType.Stochastic;
            rfi.m_deterministicRF = (ModelParameter)"@V1";

            OptionTree op = new OptionTree(rov);
            op.PayoffInfo.PayoffExpression = "Max(bond(TT;TT+tau;@v1)-strike;0)";

            op.PayoffInfo.Timing.EndingTime.m_Value = (RightValue)maturityOpt;
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
            double samplePrice = price.m_Value;

            double sampleDevSt = price.m_StdErr / Math.Sqrt((double)n_sim);

            // Calculation of the theoretical value of the call.
            CapHW1 cap = new CapHW1(zerorate);

            double d1 = cap.D1(a1, sigma1, maturityOpt, maturityOpt + tau, strike);
            double d2 = cap.D2(a1, sigma1, maturityOpt, maturityOpt + tau, strike);
            double theoreticalPrice = ZCB(zerorate, maturityOpt + tau) * SpecialFunctions.NormCdf(d1) - strike * ZCB(zerorate, maturityOpt) * SpecialFunctions.NormCdf(d2);

            Console.WriteLine("Theoretical Price = " + theoreticalPrice.ToString());
            Console.WriteLine("Monte Carlo Price = " + samplePrice);
            Console.WriteLine("Standard Deviation = " + sampleDevSt.ToString());
            double tol = 4.0 * sampleDevSt;

            Assert.Less(Math.Abs(theoreticalPrice - samplePrice), tol);
        }

        private double ZCB(AFunction zr, double t)
        {
            return Math.Exp(-zr.Evaluate(t) * t);
        }
    }
}
