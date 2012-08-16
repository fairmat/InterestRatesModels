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
using System.Text;

namespace Pelsser
{
    /// <summary>
    /// Dictionary Key for accessing <see cref="PelsserCache"/>,
    /// which is optimized to memorize integrals parametrized by two dates.
    /// </summary>
    public struct PelsserKey
    {
        /// <summary>
        /// The first date parameter which will act as index togheter with s.
        /// </summary>
        public double t;

        /// <summary>
        /// The second date parameter which will act as index togheter with t.
        /// </summary>
        public double s;

        /// <summary>
        /// Default constructor to create a new key for accessing the <see cref="PelsserCache"/>.
        /// </summary>
        /// <param name="t">The first date parameter.</param>
        /// <param name="s">The second date parameter.</param>
        public PelsserKey(double t, double s)
        {
            this.t = t;
            this.s = s;
        }

        /// <summary>
        /// Gets the HashCode of this <see cref="PelsserKey"/> object.
        /// </summary>
        /// <returns>The calculated HashCode for this object.</returns>
        public override int GetHashCode()
        {
            return this.s.GetHashCode() + 1234 * this.t.GetHashCode();
        }

        /// <summary>
        /// Checks if the provided object is equivalent this object.
        /// </summary>
        /// <param name="obj">The object to check for equivalence.</param>
        /// <returns>True if the two objects are equivalent.</returns>
        public override bool Equals(object obj)
        {
            return Equals((PelsserKey)obj);
        }

        /// <summary>
        /// Checks if the provided <see cref="PelsserKey"/> object is equivalent
        /// to this object.
        /// </summary>
        /// <param name="obj">True if the two objects are equivalent.</param>
        /// <returns>True if the two <see cref="PelsserKey"/> objects are equivalent.</returns>
        public bool Equals(PelsserKey obj)
        {
            return obj.t == this.t && obj.s == this.s;
        }
    }

    /// <summary>
    /// Caches some data for use by Pelsser, in order to avoid recalculating it several times.
    /// </summary>
    public class PelsserCache
    {
        /// <summary>
        /// A reference to the <see cref="SquaredGaussianModel"/> which is used with this cache.
        /// </summary>
        private SquaredGaussianModel instance;

        /// <summary>
        /// Gets or sets value to be stored for A(t,s).
        /// </summary>
        public double A { get; set; }

        /// <summary>
        /// Gets or sets value to be stored for B(t,s).
        /// </summary>
        internal double B { get; set; }

        /// <summary>
        /// Gets or sets value to be stored for C(s,t).
        /// </summary>
        internal double CtT0 { get; set; }

        /// <summary>
        /// Constructor to create a new instance of the PelsserCache.
        /// </summary>
        /// <param name="t">The first date parameter this element is referenced to.</param>
        /// <param name="s">The second date parameter this element is referenced to.</param>
        /// <param name="p_instance">
        /// The instance to the <see cref="SquaredGaussianModel"/> this cache references to.
        /// </param>
        public PelsserCache(double t, double s, SquaredGaussianModel p_instance)
        {
            this.instance = p_instance;

            double[] p;
            double[] btT;

            // Integral calculation is always done with daily intervals.
            int ipy = 252;

            double cached_dt = this.instance.CacheDates[1] - this.instance.CacheDates[0];

            // Indices representing t and s on the mDates discretization.
            if (this.instance.CacheDates[this.instance.CacheDates.Length - 1] < s || 1.0 / cached_dt < ipy)
            {
                double dt = 1.0 / ipy;
                double[] newDates = new double[(int)(1 + s * ipy)];
                for (int j = 0; j < newDates.Length; j++)
                    newDates[j] = j * dt;

                this.instance.CalculateValueForCache(newDates);
            }

            int ti = DVPLDOM.AdaptiveTimeDiscretization.DiscreteTime(t, this.instance.CacheDates);
            int si = DVPLDOM.AdaptiveTimeDiscretization.DiscreteTime(s, this.instance.CacheDates);
            double delta = this.instance.CacheDates[1] - this.instance.CacheDates[0];

            p = new double[si - ti + 1 + 1];
            p[0] = t;
            for (int j = 1; j < p.Length; j++)
            {
                p[j] = p[j - 1] + delta;
            }

            btT = this.instance.B(ti, si, delta);
            this.A = this.instance.A(ti, si, delta, btT);
            if (btT.Length > 0)
            {
                this.B = btT[0];
            }

            this.CtT0 = this.instance.C(s - t);
        }
    }
}
