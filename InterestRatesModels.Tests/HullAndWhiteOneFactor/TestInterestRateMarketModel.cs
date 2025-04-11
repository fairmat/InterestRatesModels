/* Copyright (C) 2009-2012 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
 * Author(s): 
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

using DVPLI.Enums;
using Fairmat.Finance;

namespace HullAndWhiteOneFactor;
using System;
using DVPLI;
using NUnit.Framework;
using DVPLDOM;
using DVPLI;
using DVPLSolver;
using Fairmat.Math;
using Mono.Addins;

/// <summary>
/// This test compares the price of a call calculated in c# with a benchmark value.
/// </summary>
[TestFixture]
public class TestIRMarketModelPricer
{
    static double atol = 1e-4;


    [SetUp]
    public void Init()
    {
        TestCommon.TestInitialization.CommonInitialization();
    }

    [Test]
    public void TestIRCap()
    {
        var maturity = 7;

        Document doc = new Document();
        ProjectROV rov1 = new ProjectROV(doc);
        int nSim = 20000;
        int nSteps = 1024;
        double a = 0.2;
        double DR = 0.02;
        double r0 = 0.015;
        double a1 = 0.002;
        double sigma1 = 0.01;
        double strike = 0.005;
        double tau = 0.5;

        ModelParameter PT = new ModelParameter(maturity, "TT");
        PT.VarName = "TT";
        rov1.Symbols.Add(PT);

        ModelParameter Ptau = new ModelParameter(tau, "tau");
        Ptau.VarName = "tau";
        rov1.Symbols.Add(Ptau);

        ModelParameter Pa = new ModelParameter(a, "a");
        Pa.VarName = "a";
        rov1.Symbols.Add(Pa);

        ModelParameter PDR = new ModelParameter(DR, "PDR");
        PDR.VarName = "DR";
        rov1.Symbols.Add(PDR);

        ModelParameter Pr0 = new ModelParameter(r0, "r0");
        Pr0.VarName = "r0";
        rov1.Symbols.Add(Pr0);



        var tenors = (Vector)(new double[] { 0.5, 1, 2, 3, 4, 5, 6, 10 });
        var rates = (Vector)(new double[] { 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045 });
        Function zerorate = new PFunction(tenors, rates);
        zerorate.VarName = "ZeroRate";
        rov1.Symbols.Add(zerorate);

        // Add to the project the created zero rate curve.

        var rfi = rov1.GetDiscountingModel();

        if (rfi is RiskFreeInfo rfiRiskInfo)
        {
            rfiRiskInfo.ActualizationType = EActualizationType.ZeroCoupond;
            rfiRiskInfo.m_deterministicRF.m_Value =  (RightValue)("exp( -ZeroRate(t)*t)");
        }



        ModelParameter Pstrike = new ModelParameter(strike, "strike");
        Pstrike.VarName = "strike";
        rov1.Symbols.Add(Pstrike);
        var payoff = "tau * max(@V1 - strike; 0)";

        doc.Part.Add(rov1);
        
        // Create the short rate process.
        var hw1 = new HW1(a1, sigma1, "@ZeroRate");

        StochasticProcessExtendible hw = new StochasticProcessExtendible(rov1, hw1);
        
        rov1.Processes.AddProcess(hw);

        OptionTree ot = new OptionTree(rov1);
        
        ot.PayoffInfo.PayoffExpression = payoff;
        
        ot.PayoffInfo.Timing.EndingTime.m_Value = (RightValue)maturity;
        ot.PayoffInfo.European = true;        
        rov1.Map.Root = ot;
        

        rov1.NMethods.Technology = ETechType.T_SIMULATION;
        rov1.NMethods.PathsNumber = nSim;
        rov1.NMethods.SimulationSteps = nSteps;


        rov1.Initialize();
        rov1.CreateSymbols();
        rov1.Parse();
        Project.ActiveProject = rov1;
        Option.CurrentSolving = ot;

        var solver = new ROVSolver();
        solver.BindToProject(rov1);
        //solver.Solve(ot);
        //Project.r_Project = 2


        //ResultItem price = rov1.m_ResultList[0] as ResultItem;
        //double monteCarloPrice = price.value;
        rfi = rov1.GetDiscountingModel();
        var df = rfi.GetDeterministicDiscountFactor(rov1, 0, maturity);
        var df_Expected = Math.Exp((-maturity * zerorate.Evaluate(maturity)));
        Assert.AreEqual(df_Expected, df, 1e-5);

        // Calculate the caplet price using the Bachelier model
        var bachelierPrice = hw1.Caplet(0, strike, tau, maturity, maturity, ResetType.Arrears, InterestRateMarketModel.BachelierNormalModel, null);

        var theoreticalPrice = new CapHW1(zerorate).HWCap(a1, sigma1, maturity, maturity + tau, strike);



        // Assert the monte carlo price
        Assert.AreEqual(theoreticalPrice, bachelierPrice.MarkToMarket, atol);
    }
}