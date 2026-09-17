using DVPLDOM;
using DVPLI;
using NUnit.Framework;

namespace CIRProcess
{
    /// <summary>
    /// Tests for <see cref="CIR2Choice"/>: a Mono.Addins plugin-registration descriptor with
    /// no real model logic - only Description and CreateInstance() are worth asserting on.
    /// </summary>
    [TestFixture]
    public class TestCIR2Choice
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void Description_ReturnsExpectedText()
        {
            CIR2Choice choice = new CIR2Choice();

            Assert.That(choice.Description, Is.EqualTo("Interest Rate Models/CIR (Two Factors)"));
        }

        [Test]
        public void CreateInstance_ReturnsStochasticProcessExtendibleWrappingCIR2()
        {
            CIR2Choice choice = new CIR2Choice();

            IEditable result = choice.CreateInstance();

            Assert.That(result, Is.InstanceOf<StochasticProcessExtendible>());
            Assert.That(((StochasticProcessExtendible)result).process, Is.InstanceOf<CIR2>());
        }
    }
}
