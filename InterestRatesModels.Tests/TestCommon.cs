/* Copyright (C) 2012 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
 * Author(s): Francesco Biondi (francesco.biondi@fairmat.com)
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

namespace TestCommon
{
    /// <summary>
    /// Contains some methods used to initialize unit tests for Fairmat.
    /// </summary>
    public static class TestInitialization
    {
        /// <summary>
        /// Generic system initialization.
        /// </summary>
        public static void CommonInitialization()
        {
            DVPLI.PluginsManager.Init();
            Mono.Addins.AddinManager.Registry.ResetConfiguration();
            Mono.Addins.AddinManager.Registry.Update(new Mono.Addins.ConsoleProgressStatus(0));

            // Since some of the variables of the tests are thread specific reinitialize the parser
            // on each test (NUnit might launch different tests on different threads)
            // This should also ensure that the parser is clean on each test
            DVPLI.Engine.Parser.NewContext();
        }
    }
}
