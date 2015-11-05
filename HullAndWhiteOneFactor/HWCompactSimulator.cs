/* Copyright (C) 2015-2015 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
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

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using DVPLI;

namespace HullAndWhiteOneFactor
{
    /// <summary>
    /// Implementation of a simple Monte Carlo Simulator for the Hull and White model.
    /// </summary>
    internal class HWCompactSimulator
    {
        internal IFunction zr;
        internal double a;
        internal double sigma;

        int I = 12 * 40;
        int P = 600;
        Matrix epsilon;
        internal unsafe HWCompactSimulator()
        {
            epsilon = new Matrix(P, I);
            Engine.Generator.Normal(epsilon.Buffer, P * I);
        }

        unsafe internal void Simulate(double m, out List<double> dates, out List<double> finalRate, out List<double> avgRate)
        {
            List<double> avg = new List<double>();
            dates = new List<double>();
            finalRate = new List<double>();
            avgRate = new List<double>();
            Vector r = new Vector(P);
            double* _r = r.Buffer;

            double dt = m / I;
            double rdt = Math.Sqrt(dt);
            for (int i = 0; i < I; i++)
            {
                double theta = Theta(i * dt, dt);
                for (int p = 0; p < P; p++)
                    _r[i] += (theta - a * _r[p]) * dt + sigma * epsilon[p, i] * rdt;

                finalRate.Add(r.Mean());
                avgRate.Add(((Vector)(finalRate.ToArray())).Mean());
                dates.Add(i * dt);
            }

            return;
        }


        protected double Theta(double t, double dt)
        {
            return Ft(t, dt) + a * F(t, dt) + (1.0 - Math.Exp(-2.0 * a * t)) * sigma * sigma / (2.0 * a);
        }
        double F(double t, double dt)
        {
            double zrT = zr.Evaluate(t);
            return t * (zr.Evaluate(t + dt) - zrT) / dt + zrT;
        }
        double Ft(double t, double dt)
        {
            return (F(t + dt, dt) - F(t - dt, dt)) / (2 * dt);
        }
        internal double ExpectedShortRate(double t)
        {
            double term1 = Math.Exp(-a * t) * zr.Evaluate(0);
            double ds = 0.001;
            double term2 = 0;
            for (double s = 0; s <= t; s += ds)
                term2 += Math.Exp(a * (s - t)) * Theta(s, ds) * ds;
            return term1 + term2;
        }

        internal double ExpectedAverageRate(double t)
        {
            List<double> avg = new List<double>();
            double ds = 0.001;
            double term2 = 0;
            for (double s = 0; s <= t; s += ds)
            {
                double term1 = Math.Exp(-a * s) * zr.Evaluate(0);
                term2 += Math.Exp(a * (s - t)) * Theta(s, ds) * ds;
                avg.Add(term1 + term2);
            }
            var v = (Vector)(avg.ToArray());
            return v.Mean();
        }

    }
}
