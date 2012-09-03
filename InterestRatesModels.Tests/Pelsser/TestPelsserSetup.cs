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

using DVPLDOM;
using DVPLI;
using DVPLSolver;
using NUnit.Framework;

namespace Pelsser
{
    [TestFixture]
    public class TestPelsserSetup
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void Test()
        {
            Document doc = new Document();
            ProjectROV rov = new ProjectROV(doc);
            doc.Part.Add(rov);

            AFunction zerorate = new AFunction(rov);
            zerorate.VarName = "zr";
            zerorate.m_IndependentVariables = 1;
            zerorate.m_Value = (RightValue)"0.05";

            rov.Symbols.Add(zerorate);

            int n_sim = 5000;
            int n_steps = 900;
            SquaredGaussianModel process = new SquaredGaussianModel();
            process.a1 = (ModelParameter)0.0704877770928728;
            process.sigma1 = (ModelParameter)0.116251516678772;
            process.zr = (ModelParameter)"@zr";

            StochasticProcessExtendible s = new StochasticProcessExtendible(rov, process);
            rov.Processes.AddProcess(s);

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
        }
    }
}
