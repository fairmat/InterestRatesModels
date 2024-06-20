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

using System.Reflection;
using Mono.Addins;

// The following lines tell that the assembly is an addin.
[assembly: Addin("Pelsser Squared Gaussian Model", "1.0.22", Category = "Stochastic Process")]
[assembly: AddinDependency("Fairmat", "1.0")]
[assembly: AddinAuthor("Fairmat SRL")]
[assembly: AddinDescription("The one-factor squared Gaussian model. " +
                            "This model assumes that the spot interest rate is a " +
                            "quadratic function of the underlying process, and it " +
                            "provides the advantage that the interest rates never become " +
                            "negative. The squared Gaussian is a no-arbitrage model, so it " +
                            "can be fitted to the initial term-structure of interest rates.")]

[assembly: AssemblyTrademark("Fairmat")]
[assembly: AssemblyCulture("")]
