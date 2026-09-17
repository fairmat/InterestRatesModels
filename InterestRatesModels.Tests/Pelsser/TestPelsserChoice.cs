using DVPLDOM;
using DVPLI;
using NUnit.Framework;

namespace Pelsser
{
    /// <summary>
    /// Tests for <see cref="PelsserChoice"/>: a Mono.Addins plugin-registration descriptor
    /// with no real model logic - only Description and CreateInstance() are worth asserting on.
    /// </summary>
    [TestFixture]
    public class TestPelsserChoice
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void Description_ReturnsExpectedText()
        {
            PelsserChoice choice = new PelsserChoice();

            Assert.That(choice.Description, Is.EqualTo("Pelsser Squared Gaussian Model"));
        }

        [Test]
        public void CreateInstance_ReturnsStochasticProcessExtendibleWrappingSquaredGaussianModel()
        {
            PelsserChoice choice = new PelsserChoice();

            IEditable result = choice.CreateInstance();

            Assert.That(result, Is.InstanceOf<StochasticProcessExtendible>());
            Assert.That(((StochasticProcessExtendible)result).process, Is.InstanceOf<SquaredGaussianModel>());
        }
    }
}
