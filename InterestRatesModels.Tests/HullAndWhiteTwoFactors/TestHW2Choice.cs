using DVPLDOM;
using DVPLI;
using NUnit.Framework;

namespace HullAndWhiteTwoFactors
{
    /// <summary>
    /// Tests for <see cref="HW2Choice"/>: a Mono.Addins plugin-registration descriptor with
    /// no real model logic - only Description and CreateInstance() are worth asserting on.
    /// </summary>
    [TestFixture]
    public class TestHW2Choice
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void Description_ReturnsExpectedText()
        {
            HW2Choice choice = new HW2Choice();

            Assert.That(choice.Description, Is.EqualTo("Hull and White (Two Factors) model"));
        }

        [Test]
        public void CreateInstance_ReturnsStochasticProcessExtendibleWrappingHW2()
        {
            HW2Choice choice = new HW2Choice();

            IEditable result = choice.CreateInstance();

            Assert.That(result, Is.InstanceOf<StochasticProcessExtendible>());
            Assert.That(((StochasticProcessExtendible)result).process, Is.InstanceOf<HW2>());
        }
    }
}
