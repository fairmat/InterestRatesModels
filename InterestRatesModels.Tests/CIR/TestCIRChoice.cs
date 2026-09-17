using DVPLDOM;
using DVPLI;
using NUnit.Framework;

namespace CIRProcess
{
    /// <summary>
    /// Tests for <see cref="CIRChoice"/>: a Mono.Addins plugin-registration descriptor with
    /// no real model logic - only Description and CreateInstance() are worth asserting on.
    /// </summary>
    [TestFixture]
    public class TestCIRChoice
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void Description_ReturnsExpectedText()
        {
            CIRChoice choice = new CIRChoice();

            Assert.That(choice.Description, Is.EqualTo("Interest Rate Models/CIR"));
        }

        [Test]
        public void CreateInstance_ReturnsStochasticProcessExtendibleWrappingCIR()
        {
            CIRChoice choice = new CIRChoice();

            IEditable result = choice.CreateInstance();

            Assert.That(result, Is.InstanceOf<StochasticProcessExtendible>());
            Assert.That(((StochasticProcessExtendible)result).process, Is.InstanceOf<CIR>());
        }
    }
}
