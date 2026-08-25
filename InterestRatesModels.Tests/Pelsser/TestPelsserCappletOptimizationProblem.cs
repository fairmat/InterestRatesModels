using System;
using DVPLDOM;
using DVPLI;
using Fairmat.Finance;
using Fairmat.Optimization;
using NUnit.Framework;

namespace Pelsser.Calibration
{
    [TestFixture]
    public class TestPelsserCappletOptimizationProblem
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        private static ProjectROV CreateProjectWithModel(out Function zr)
        {
            Document doc = new Document();
            ProjectROV prj = new ProjectROV(doc);
            doc.Part.Add(prj);

            zr = new PFunction(null);
            zr.VarName = "zr";
            zr.Expr = new double[,] { { 0, 0.02 }, { 10, 0.02 } };
            prj.Symbols.Add(zr);

            SquaredGaussianModel model = new SquaredGaussianModel();
            model.a1 = (ModelParameter)0.1;
            model.sigma1 = (ModelParameter)0.01;
            model.zr = (ModelParameter)"@zr";
            StochasticProcessExtendible iex = new StochasticProcessExtendible(prj, model);
            prj.Processes.AddProcess(iex);

            prj.Parse();

            return prj;
        }

        private static PelsserCappletOptimizationProblem CreateProblem(ProjectROV prj, Function zr,
            out Vector capMaturity, out Vector fwd, out Vector capK, out double deltaK, out Vector capMat, out Matrix blackCaps)
        {
            deltaK = 0.5;
            capMat = new Vector(new double[] { 1.0 });
            capK = new Vector(new double[] { 0.02 });

            capMaturity = new Vector((int)(1.0 + capMat[capMat.Length - 1] / deltaK));
            for (int l = 0; l < capMaturity.Length; l++)
                capMaturity[l] = deltaK * l;

            BlackModel bm = new BlackModel(zr);
            fwd = new Vector(capMaturity.Length - 1);
            for (int i = 0; i < fwd.Length; i++)
                fwd[i] = bm.Fk(capMaturity[i + 1], deltaK);

            blackCaps = new Matrix(new double[,] { { 0.001 } });

            Caplet cp = new Caplet();
            return new PelsserCappletOptimizationProblem(prj, cp, capMaturity, fwd, capK, deltaK, capMat, blackCaps);
        }

        [Test]
        public void BoundsAreFixed()
        {
            Function zr;
            ProjectROV prj = CreateProjectWithModel(out zr);
            Vector capMaturity, fwd, capK, capMat;
            Matrix blackCaps;
            double deltaK;
            PelsserCappletOptimizationProblem problem = CreateProblem(prj, zr, out capMaturity, out fwd, out capK, out deltaK, out capMat, out blackCaps);

            Bounds b = problem.Bounds;

            CollectionAssert.AreEqual(new double[] { 1e-8, 1e-6 }, (double[])b.Lb.ToArray());
            CollectionAssert.AreEqual(new double[] { 1, .1 }, (double[])b.Ub.ToArray());
        }

        [Test]
        public void HasNoNonLinearOrLinearConstraints()
        {
            Function zr;
            ProjectROV prj = CreateProjectWithModel(out zr);
            Vector capMaturity, fwd, capK, capMat;
            Matrix blackCaps;
            double deltaK;
            PelsserCappletOptimizationProblem problem = CreateProblem(prj, zr, out capMaturity, out fwd, out capK, out deltaK, out capMat, out blackCaps);

            Assert.IsFalse(problem.HasNonLinearConstraints);
            Assert.IsNull(problem.LinearIneqConstraints);
        }

        [Test]
        public void GradIsNotImplemented()
        {
            Function zr;
            ProjectROV prj = CreateProjectWithModel(out zr);
            Vector capMaturity, fwd, capK, capMat;
            Matrix blackCaps;
            double deltaK;
            PelsserCappletOptimizationProblem problem = CreateProblem(prj, zr, out capMaturity, out fwd, out capK, out deltaK, out capMat, out blackCaps);

            Assert.Throws<NotImplementedException>(() => problem.Grad(new Vector(new double[] { 0.1, 0.01 })));
        }

        [Test]
        public void GIsNotImplemented()
        {
            Function zr;
            ProjectROV prj = CreateProjectWithModel(out zr);
            Vector capMaturity, fwd, capK, capMat;
            Matrix blackCaps;
            double deltaK;
            PelsserCappletOptimizationProblem problem = CreateProblem(prj, zr, out capMaturity, out fwd, out capK, out deltaK, out capMat, out blackCaps);

            Assert.Throws<NotImplementedException>(() => problem.G(new Vector(new double[] { 0.1, 0.01 })));
        }

        [Test]
        public void AssignSetsA1AndSigma1OnTheProjectModel()
        {
            Function zr;
            ProjectROV prj = CreateProjectWithModel(out zr);

            SquaredGaussianModel model = PelsserCappletOptimizationProblem.Assign(prj, new Vector(new double[] { 0.2, 0.05 }));

            Assert.AreEqual(0.2, model.a1.fV(), 1e-12);
            Assert.AreEqual(0.05, model.sigma1.fV(), 1e-12);
        }

        [Test]
        public void PelsserConstraintReturnsFiniteSingleValueVector()
        {
            Function zr;
            ProjectROV prj = CreateProjectWithModel(out zr);

            Vector result = PelsserCappletOptimizationProblem.PelsserConstraint(prj, new Vector(new double[] { 0.1, 0.01 }), 1.0);

            Assert.AreEqual(1, result.Length);
            Assert.IsTrue(double.IsFinite(result[0]));
        }

        [Test]
        public void ObjReturnsFiniteValue()
        {
            Function zr;
            ProjectROV prj = CreateProjectWithModel(out zr);
            Vector capMaturity, fwd, capK, capMat;
            Matrix blackCaps;
            double deltaK;
            PelsserCappletOptimizationProblem problem = CreateProblem(prj, zr, out capMaturity, out fwd, out capK, out deltaK, out capMat, out blackCaps);

            double result = problem.Obj(new Vector(new double[] { 0.1, 0.01 }));

            Assert.IsTrue(double.IsFinite(result));
        }
    }
}
