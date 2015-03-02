/* Copyright (C) 2009-2012 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
 * Author(s): Matteo Tesser (matteo.tesser@fairmat.com)
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

// The following lines tell that the assembly is an addin
[assembly: Addin("Hull-White one factor model", "1.1.4", Category = "Stochastic Process")]
[assembly: AddinDependency("Fairmat", "1.0")]
[assembly: AddinAuthor("Fairmat SRL")]
[assembly: AddinDescription("A no-arbitrage model which is the industry standard " +
                               "for modeling the future interest rate dynamic.")]

// General Information about an assembly is controlled through the following
// set of attributes. Change these attribute values to modify the information
// associated with an assembly.
[assembly: AssemblyTitle("Hull-WhiteOneFactorModel")]
[assembly: AssemblyDescription("A no-arbitrage model which is the industry standard " +
                               "for modeling the future interest rate dynamic.")]
[assembly: AssemblyConfiguration("")]
[assembly: AssemblyProduct("Hull-WhiteOneFactorModel")]
[assembly: AssemblyCompany("Fairmat SRL")]
[assembly: AssemblyCopyright("Copyright © Fairmat SRL 2009-2015")]
[assembly: AssemblyTrademark("Fairmat")]
[assembly: AssemblyCulture("")]

// Setting ComVisible to false makes the types in this assembly not visible
// to COM components.  If you need to access a type in this assembly from
// COM, set the ComVisible attribute to true on that type.
[assembly: ComVisible(false)]

// The following GUID is for the ID of the typelib if this project is exposed to COM
[assembly: Guid("82a890f6-c19a-42bb-a810-2d84c312e77e")]

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
[assembly: AssemblyVersion("1.1.4")]
[assembly: AssemblyFileVersion("1.1.4")]
