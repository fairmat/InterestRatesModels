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

using DVPLDOM;
using DVPLI;
using Mono.Addins;

namespace CIRProcess
{
    /// <summary>
    /// Implements the IEditableChoice interface for <see cref="CIR"/>,
    /// in order to allow displaying an user interface.
    /// </summary>
    [Extension("/Fairmat/ProcessTypeChoice")]
    public class CIRChoice : IEditableChoice
    {
        #region IEditableOption Members

        /// <summary>
        /// Gets the name of the model which will be shown to the user.
        /// </summary>
        public string Description
        {
            get
            {
                return "Interest Rate Models/" + CIR.CIRLabel;
            }
        }

        /// <summary>
        /// Creates an IEditable instance from a StochasticProcessExtendible,
        /// which will handle the CIR plugin (<see cref="CIR"/> part).
        /// </summary>
        /// <returns>A reference to a new IEditable instance.</returns>
        public IEditable CreateInstance()
        {
            return new StochasticProcessExtendible(null, new CIR());
        }

        #endregion
    }
}
