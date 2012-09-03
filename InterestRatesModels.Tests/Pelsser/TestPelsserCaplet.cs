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

namespace Pelsser.Calibration
{
    [TestFixture]
    public class TestPelsserCaplet
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        /// <summary>
        /// Calculate the price of a Caplet through a Monte Carlo simulation and
        /// compare it with the theoretical value.
        /// </summary>
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

            int n_sim = 4000;
            double MaturityOpt = 6.5;

            // Simulation steps for a year. With StepPerYear = 150 the test will be passed.
            // But notice that the price calculated through Monte Carlo is unstable when
            // changing this value, even till 1000 steps per year.
            int StepsPerYear = 150;
            int n_steps = StepsPerYear * ((int)MaturityOpt);
            double strike = 0.01;
            double tau = 0.5;

            SquaredGaussianModel process = new SquaredGaussianModel();
            process.a1 = (ModelParameter)0.1;
            process.sigma1 = (ModelParameter)0.01;
            process.zr = (ModelParameter)"@zr";

            StochasticProcessExtendible s = new StochasticProcessExtendible(rov, process);
            rov.Processes.AddProcess(s);

            ModelParameter PT = new ModelParameter(MaturityOpt, "TT");
            PT.VarName = "TT";
            rov.Symbols.Add(PT);

            ModelParameter Ptau = new ModelParameter(tau, "tau");
            Ptau.VarName = "tau";
            rov.Symbols.Add(Ptau);

            ModelParameter Pstrike = new ModelParameter(strike, "strike");
            Pstrike.VarName = "strike";
            rov.Symbols.Add(Pstrike);

            // Set the discounting.
            RiskFreeInfo rfi = rov.GetDiscountingModel() as RiskFreeInfo;
            rfi.ActualizationType = EActualizationType.Stochastic;
            rfi.m_deterministicRF = (ModelParameter)"@V1";

            // Set the payoff.
            OptionTree op = new OptionTree(rov);
            op.PayoffInfo.PayoffExpression = "tau*max(rate(TT;tau;@v1) - strike; 0)";
            op.PayoffInfo.Timing.EndingTime.m_Value = (RightValue)(MaturityOpt + tau);
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
                rov.DisplayErrors();
            }

            Assert.IsFalse(rov.HasErrors);

            ResultItem price = rov.m_ResultList[0] as ResultItem;
            double MCPrice = price.m_Value;
            double MCDevST = price.m_StdErr / Math.Sqrt((double)n_sim);

            Caplet cplt = new Caplet();
            Vector Mat, fwd, Rk, CapMatV;
            double delta_k, CapMat;
            delta_k = 0.5;
            CapMat = MaturityOpt + tau;
            int nmat = 2 * ((int)CapMat) + 1;
            Mat = new Vector(nmat);
            fwd = new Vector(nmat);
            Mat[0] = 0;
            fwd[0] = zerorate.Evaluate(0);
            for (int k = 1; k < nmat; k++)
            {
                Mat[k] = tau * ((double)k);
                fwd[k] = zerorate.Evaluate(Mat[k]) * Mat[k] - zerorate.Evaluate(Mat[k - 1]) * Mat[k - 1];
            }

            fwd = fwd / tau;
            Rk = new Vector(1);
            Rk[0] = strike;
            CapMatV = new Vector(2);
            CapMatV[0] = MaturityOpt;
            CapMatV[1] = MaturityOpt + tau;
            Matrix Cap = cplt.PGSMCaplets(process, Mat, fwd, Rk, delta_k, CapMatV);

            double ThPrice = Cap[1, 0] - Cap[0, 0];

            Console.WriteLine("\nTheoretical Price = " + ThPrice.ToString());
            Console.WriteLine("Monte Carlo Price = " + MCPrice);
            Console.WriteLine("Standard Deviation = " + MCDevST);

            double tol = 4.0 * MCDevST;
            Assert.Less(Math.Abs(ThPrice - MCPrice), tol);
        }
    }
}
