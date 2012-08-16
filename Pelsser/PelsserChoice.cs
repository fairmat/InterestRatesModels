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
using System.Collections.Generic;
using System.Runtime.Serialization;
using System.Text;
using DVPLDOM;
using DVPLI;
using Mono.Addins;
using ParallelVectors;

namespace Pelsser
{
    /// <summary>
    /// Implements the interfaces which allows the model to
    /// be instanced from the user interface.
    /// </summary>
    [Extension("/Fairmat/ProcessTypeChoice")]
    public class PelsserChoice : IEditableChoice
    {
        #region IEditableOption Members

        /// <summary>
        /// Gets the name of the model which will be shown to the user.
        /// </summary>
        public string Description
        {
            get
            {
                return "Pelsser Squared Gaussian Model";
            }
        }

        /// <summary>
        /// Creates an IEditable instance from a StochasticProcessExtendible,
        /// which will handle the Pelsser plug-in.
        /// </summary>
        /// <returns>A reference to a new IEditable instance.</returns>
        public IEditable CreateInstance()
        {
            return new StochasticProcessExtendible(null, new SquaredGaussianModel());
        }

        #endregion
    }
}
