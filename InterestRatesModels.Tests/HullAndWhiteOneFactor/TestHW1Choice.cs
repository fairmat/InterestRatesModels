using DVPLDOM;
using DVPLI;
using NUnit.Framework;

namespace HullAndWhiteOneFactor
{
    /// <summary>
    /// Tests for <see cref="HW1Choice"/>: a Mono.Addins plugin-registration descriptor with
    /// no real model logic - only Description and CreateInstance() are worth asserting on.
    /// </summary>
    [TestFixture]
    public class TestHW1Choice
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void Description_ReturnsExpectedText()
        {
            HW1Choice choice = new HW1Choice();

            Assert.That(choice.Description, Is.EqualTo("Hull and White (One Factor) model"));
        }

        [Test]
        public void CreateInstance_ReturnsStochasticProcessExtendibleWrappingHW1()
        {
            HW1Choice choice = new HW1Choice();

            IEditable result = choice.CreateInstance();

            Assert.That(result, Is.InstanceOf<StochasticProcessExtendible>());
            Assert.That(((StochasticProcessExtendible)result).process, Is.InstanceOf<HW1>());
        }
    }
}
