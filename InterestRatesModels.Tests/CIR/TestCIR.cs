/* Copyright (C) 2009-2012 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
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
using DVPLI;
using NUnit.Framework;

namespace CIRProcess
{
    /// <summary>
    /// Tests a bond call through CIR.
    /// </summary>
    [TestFixture]
    public class TestBondCall
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void Test()
        {
            double k = 1.0;
            double theta = 0.02;
            double sigma = 0.08;
            double r0 = 0.01;
            Vector par = new Vector(3);
            par[0] = k;
            par[1] = theta;
            par[2] = sigma;

            double strike = 0.98;
            double noz = 100.0;

            double T = 1.0;
            double S = 2.0;

            double callFairmat = noz * CIRCap.BondCall(r0, 0.0, T, S, strike, par);
            double callBenchmark = 0.317304505282290;

            Console.WriteLine("CallFairmat   = " + callFairmat);
            Console.WriteLine("CallBenchmark = " + callBenchmark);

            double maxError = 1e-10;
            Assert.Less(Math.Abs(callFairmat - callBenchmark), maxError);
        }
    }
}
