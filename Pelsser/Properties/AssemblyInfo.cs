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
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using Mono.Addins;

// The following lines tell that the assembly is an addin.
[assembly: Addin("Pelsser squared gaussian model", "1.0.21", Category = "Stochastic Process")]
[assembly: AddinDependency("Fairmat", "1.0")]
[assembly: AddinAuthor("Fairmat SRL")]
[assembly: AddinDescription("The one-factor squared Gaussian model. " +
                            "This model assumes that the spot interest rate is a " +
                            "quadratic function of the underlying process, and it " +
                            "provides the advantage that the interest rates never become " +
                            "negative. The squared Gaussian is a no-arbitrage model, so it " +
                            "can be fitted to the initial term-structure of interest rates.")]

// General Information about an assembly is controlled through the following
// set of attributes. Change these attribute values to modify the information
// associated with an assembly.
[assembly: AssemblyTitle("Pelsser")]
[assembly: AssemblyDescription("The one-factor squared Gaussian model. " +
                               "This model assumes that the spot interest rate is a " +
                               "quadratic function of the underlying process, and it " +
                               "provides the advantage that the interest rates never become " +
                               "negative. The squared Gaussian is a no-arbitrage model, so it " +
                               "can be fitted to the initial term-structure of interest rates.")]
[assembly: AssemblyConfiguration("")]
[assembly: AssemblyProduct("Pelsser")]
[assembly: AssemblyCompany("Fairmat SRL")]
[assembly: AssemblyCopyright("Copyright © Fairmat SRL 2009-2015")]
[assembly: AssemblyTrademark("Fairmat")]
[assembly: AssemblyCulture("")]

// Setting ComVisible to false makes the types in this assembly not visible
// to COM components.  If you need to access a type in this assembly from
// COM, set the ComVisible attribute to true on that type.
[assembly: ComVisible(false)]

// The following GUID is for the ID of the typelib if this project is exposed to COM
[assembly: Guid("54b15e16-8ee4-4f00-8d4c-356caaee45da")]

// Version information for an assembly consists of the following four values:
//
//      Major Version
//      Minor Version
//      Build Number
//      Revision
//
// You can specify all the values or you can default the Build and Revision Numbers
// by using the '*' as shown below:
// [assembly: AssemblyVersion("1.0.*")]
[assembly: AssemblyVersion("1.0.21")]
[assembly: AssemblyFileVersion("1.0.21")]
