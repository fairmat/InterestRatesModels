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
using DVPLUtils;
using Mono.Addins;
using NUnit.Framework;

namespace HullAndWhiteOneFactor
{
    /// <summary>
    /// Implementation of Caps pricing with the HW1 model.
    /// </summary>
    [TestFixture]
    public class TestHW1Cap
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void Test()
        {
            double[,] values = { { .5, .015 },
                                 { 1, .0165 },
                                 { 2, .0168 },
                                 { 3, .0172 },
                                 { 5, .0182 },
                                 { 8, .0210 },
                                 { 10, .025 },
                                 { 15, .031 },
                                 { 20, .035 },
                                 { 30, .037 },
                                 { 40, .038 },
                               };

            Function zeroratecurve = new PFunction(null);
            zeroratecurve.Expr = values;
            (zeroratecurve as PFunction).m_Function.iType = EInterpolationType.LINEAR;

            // Execute the test.
            CapHW1 hwc = new CapHW1(zeroratecurve);

            // CAP and FLOOR Tests.
            double[] results = new double[2];
            results[0] = hwc.HWCap(0.14, 0.02, 0.02, 0.5, 1);
            results[1] = hwc.HWCap(0.14, 0.02, 0.02, 0.5, 2);

            Console.WriteLine("CAP 2 col 1 row:" + results[0]);
            Console.WriteLine("CAP 2 col 2 row:" + results[1]);

            // Check the results with previously manually calculated values.
            double[] targets = { 0.00214717607719883, 0.0084939015243779 };

            double eps = 10e-6;

            for (int r = 0; r < results.Length; r++)
            {
                Assert.LessOrEqual(Math.Abs(targets[r] - results[r]), eps);
            }
        }
    }
}
