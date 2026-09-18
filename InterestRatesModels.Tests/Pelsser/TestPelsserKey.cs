using System;
using NUnit.Framework;

namespace Pelsser
{
    /// <summary>
    /// Tests for <see cref="PelsserKey"/>. Equals(object) unconditionally casts its argument
    /// to PelsserKey (a value type), so it throws instead of returning false for null or a
    /// wrong type - these tests document that real (and arguably buggy) behavior.
    /// </summary>
    [TestFixture]
    public class TestPelsserKey
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void Equals_SameTAndS_ReturnsTrue()
        {
            PelsserKey key1 = new PelsserKey(1.5, 2.5);
            PelsserKey key2 = new PelsserKey(1.5, 2.5);

            Assert.That(key1.Equals(key2), Is.True);
            Assert.That(key1.Equals((object)key2), Is.True);
        }

        [Test]
        public void Equals_DifferentT_ReturnsFalse()
        {
            PelsserKey key1 = new PelsserKey(1.5, 2.5);
            PelsserKey key2 = new PelsserKey(1.6, 2.5);

            Assert.That(key1.Equals(key2), Is.False);
        }

        [Test]
        public void Equals_DifferentS_ReturnsFalse()
        {
            PelsserKey key1 = new PelsserKey(1.5, 2.5);
            PelsserKey key2 = new PelsserKey(1.5, 2.6);

            Assert.That(key1.Equals(key2), Is.False);
        }

        [Test]
        public void GetHashCode_EqualKeys_AreConsistent()
        {
            PelsserKey key1 = new PelsserKey(1.5, 2.5);
            PelsserKey key2 = new PelsserKey(1.5, 2.5);

            Assert.That(key1.GetHashCode(), Is.EqualTo(key2.GetHashCode()));
        }

        [Test]
        public void Equals_Null_ThrowsNullReferenceException()
        {
            PelsserKey key = new PelsserKey(1.5, 2.5);

            Assert.Throws<NullReferenceException>(() => key.Equals(null));
        }

        [Test]
        public void Equals_WrongType_ThrowsInvalidCastException()
        {
            PelsserKey key = new PelsserKey(1.5, 2.5);

            Assert.Throws<InvalidCastException>(() => key.Equals("not a key"));
        }
    }
}
