/* Copyright (C) 2010-2012 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
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
using DVPLDOM;
using DVPLI;
using DVPLSolver;
using Fairmat.Statistics;
using NUnit.Framework;

namespace TestVari
{
    /// <summary>
    /// This test compares the price of a call option on a zero coupon
    /// obtained through Monte Carlo simulation and the one
    /// obtained through analytical results.
    /// </summary>
    [TestFixture]
    public class TestHW2
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

            int n_sim = 10000;
            int n_steps = 512;
            double a = 0.2;
            double DR = 0.02;
            double r0 = 0.015;
            double a1 = 1.0;
            double sigma1 = 0.01;
            double a2 = 0.1;
            double sigma2 = 0.0165;
            double correlation = 0.6;
            double MaturityOpt = 5.0;
            double strike = 0.927;
            double tau = 2.0;

            ModelParameter PT = new ModelParameter(MaturityOpt, "TT");
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
            HW2ProcessType pt = new HW2ProcessType();

            // Set the short rate process.
            pt = HW2ProcessType.ShortRate;
            StocasticProcessHW2 process1 = new StocasticProcessHW2(rov, pt);

            process1.zero_rate_curve = "@zr";
            process1._a = (ModelParameter)a1;
            process1._b = (ModelParameter)sigma1;
            rov.Processes.AddProcess(process1);

            // Set the mean reversion process.
            pt = HW2ProcessType.Unobservable;
            StocasticProcessHW2 process2 = new StocasticProcessHW2(rov, pt);

            process2._a = (ModelParameter)a2;
            process2._b = (ModelParameter)sigma2;
            rov.Processes.AddProcess(process2);

            // Set the correlation.
            rov.Processes.r[0, 1] = (RightValue)correlation;

            // Set the discounting.
            RiskFreeInfo rfi = rov.GetDiscountingModel() as RiskFreeInfo;
            rfi.ActualizationType = EActualizationType.Stochastic;
            rfi.m_deterministicRF = (ModelParameter)"@v1";

            OptionTree op = new OptionTree(rov);
            op.PayoffInfo.PayoffExpression = "max(bond(TT;TT+tau;@v1)-strike;0)";
            op.PayoffInfo.Timing.EndingTime.m_Value = (RightValue)MaturityOpt;
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
            double SampleMean = price.m_Value;
            double SampleDevSt = price.m_StdErr / Math.Sqrt((double)n_sim);
            double ThPrice = HW2BondCall(zerorate, MaturityOpt, MaturityOpt + tau, strike, a1, sigma1, a2, sigma2, correlation);

            Console.WriteLine("\nTheoretical Price = " + ThPrice.ToString());
            Console.WriteLine("Monte Carlo Price = " + SampleMean.ToString());
            Console.WriteLine("Standard Deviation = " + SampleDevSt.ToString());

            bool result;
            double fact = 4.0;
            result = (Math.Abs(SampleMean - ThPrice) < fact * SampleDevSt);

            Assert.IsTrue(result);
        }

        private double HW2BondCall(AFunction zr, double t, double s, double K, double a, double sigma1, double b, double sigma2, double rho)
        {
            double a1 = SpecialFunctions.NormCdf(d1(zr, a, sigma1, b, sigma2, rho, t, s, K));
            double a2 = SpecialFunctions.NormCdf(d2(zr, a, sigma1, b, sigma2, rho, t, s, K));
            return ZCB(zr, s) * a1 - K * ZCB(zr, t) * a2;
        }

        private double d1(AFunction zr, double a, double sigma1, double b, double sigma2, double rho, double t, double s, double K)
        {
            double sp = SigmaP2(a, sigma1, b, sigma2, rho, t, s);
            return (Math.Log(ZCB(zr, s) / (ZCB(zr, t) * K)) / sp + sp / 2.0);
        }

        private double d2(AFunction zr, double a, double sigma1, double b, double sigma2, double rho, double t, double s, double K)
        {
            double sp = SigmaP2(a, sigma1, b, sigma2, rho, t, s);
            return (d1(zr, a, sigma1, b, sigma2, rho, t, s, K) - sp);
        }

        private double SigmaP(double a, double sigma1, double b, double sigma2, double rho, double t, double s)
        {
            // This is calculated by following the Technical note number 14 by
            // J.Hull "Options, futures and other derivatives"
            // note: there is some error, it gives negative values for sigma squared.
            double SigmaSquared, U, V, term1, term2, term3;
            U = (Math.Exp(-a * s) - Math.Exp(-a * t)) / (a * (a - b));
            V = (Math.Exp(-b * s) - Math.Exp(-b * t)) / (b * (a - b));
            term1 = sigma1 * sigma1 * Math.Pow(B(s - t, a), 2.0) * (1.0 - Math.Exp(-2.0 * a * t)) / (2.0 * a);
            term2 = sigma2 * sigma2 * (U * U * (Math.Exp(2.0 * a * t) - 1.0) / (2.0 * a) + V * V * (Math.Exp(2.0 * b * t) - 1.0) / (2.0 * b) - 2.0 * U * V * (Math.Exp(2.0 * (a + b) * t) - 1.0) / (a + b));
            term3 = 2.0 * rho * sigma1 * sigma2 * (Math.Exp(-a * t) - Math.Exp(-a * s)) * (U * (Math.Exp(2.0 * a * t) - 1.0) / (2.0 * a) - V * (Math.Exp((a + b) * t) - 1.0) / (a + b)) / a;
            SigmaSquared = term1 + term2 + term3;
            return Math.Sqrt(SigmaSquared);
        }

        private double SigmaP2(double a, double sigma1, double b, double sigma2, double rho, double t, double s)
        {
            // This is calculated by following formula 4.2 of Brigo-Mercurio "Interest rate models"
            double SigmaSquared, sigma3, sigma4, term1, term2, term3;
            sigma3 = Math.Sqrt(sigma1 * sigma1 + sigma2 * sigma2 / Math.Pow(a - b, 2.0) + 2.0 * rho * sigma1 * sigma2 / (b - a));
            sigma4 = sigma2 / (a - b);
            term1 = sigma3 * sigma3 * Math.Pow(1.0 - Math.Exp(-a * (s - t)), 2.0) * (1.0 - Math.Exp(-2.0 * a * t)) / (2 * a * a * a);
            term2 = sigma2 * sigma2 * Math.Pow(1.0 - Math.Exp(-b * (s - t)), 2.0) * (1.0 - Math.Exp(-2.0 * b * t)) / (2.0 * b * b * b * (a - b) * (a - b));
            term3 = 2.0 * (sigma1 * rho - sigma4) * sigma4 * (1.0 - Math.Exp(-a * (s - t))) * (1.0 - Math.Exp(-b * (s - t))) * (1.0 - Math.Exp(-(a + b) * t)) / (a * b * (a + b));
            SigmaSquared = term1 + term2 + term3;
            return Math.Sqrt(SigmaSquared);
        }

        private double B(double t, double a)
        {
            return (1.0 - Math.Exp(-a * t)) / a;
        }

        private double ZCB(AFunction zr, double time)
        {
            return Math.Exp(-zr.Evaluate(time) * time);
        }
    }
}
