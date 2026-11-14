"""
qsqrt5.py
Exact arithmetic in the field Q(sqrt 5) = { a + b sqrt5 : a, b in Q }.

The order is the one inherited from R.  Sign of a + b sqrt5:
  a, b of the same sign (or one of them 0): that sign;
  otherwise compare a^2 with 5 b^2 (they are never equal unless a = b = 0,
  because sqrt5 is irrational).
Conjugation sqrt5 -> -sqrt5 preserves + - * / but NOT the order, so it is
never used to decide a comparison.
"""

from fractions import Fraction


def _sgn(q):
    return (q > 0) - (q < 0)


class R5:
    __slots__ = ("a", "b")

    def __init__(self, a=0, b=0):
        self.a = Fraction(a)
        self.b = Fraction(b)

    @staticmethod
    def _co(x):
        if isinstance(x, R5):
            return x
        if isinstance(x, (int, Fraction)):
            return R5(x)
        return NotImplemented

    # ------------------------------------------------------------ field
    def __add__(self, o):
        o = R5._co(o)
        return NotImplemented if o is NotImplemented else R5(self.a + o.a, self.b + o.b)

    __radd__ = __add__

    def __neg__(self):
        return R5(-self.a, -self.b)

    def __sub__(self, o):
        o = R5._co(o)
        return NotImplemented if o is NotImplemented else R5(self.a - o.a, self.b - o.b)

    def __rsub__(self, o):
        o = R5._co(o)
        return NotImplemented if o is NotImplemented else o - self

    def __mul__(self, o):
        o = R5._co(o)
        if o is NotImplemented:
            return o
        return R5(self.a * o.a + 5 * self.b * o.b, self.a * o.b + self.b * o.a)

    __rmul__ = __mul__

    def inverse(self):
        n = self.a * self.a - 5 * self.b * self.b      # nonzero unless self == 0
        if n == 0:
            raise ZeroDivisionError("R5 division by zero")
        return R5(self.a / n, -self.b / n)

    def __truediv__(self, o):
        o = R5._co(o)
        return NotImplemented if o is NotImplemented else self * o.inverse()

    def __rtruediv__(self, o):
        o = R5._co(o)
        return NotImplemented if o is NotImplemented else o * self.inverse()

    def __pow__(self, k):
        if k < 0:
            return self.inverse() ** (-k)
        r = R5(1)
        for _ in range(k):
            r = r * self
        return r

    # ------------------------------------------------------------ order
    def sign(self):
        sa, sb = _sgn(self.a), _sgn(self.b)
        if sb == 0:
            return sa
        if sa == 0 or sa == sb:
            return sb
        # opposite signs: |a| vs sqrt5 |b|
        return sa if self.a * self.a > 5 * self.b * self.b else sb

    def _cmp(self, o):
        o = R5._co(o)
        if o is NotImplemented:
            return o
        return (self - o).sign()

    def __eq__(self, o):
        o = R5._co(o)
        return NotImplemented if o is NotImplemented else (self.a == o.a and self.b == o.b)

    def __lt__(self, o): return self._cmp(o) < 0
    def __le__(self, o): return self._cmp(o) <= 0
    def __gt__(self, o): return self._cmp(o) > 0
    def __ge__(self, o): return self._cmp(o) >= 0

    def __bool__(self):
        return bool(self.a) or bool(self.b)

    def __hash__(self):
        return hash(self.a) if self.b == 0 else hash((self.a, self.b))

    def __float__(self):
        return float(self.a) + float(self.b) * 5 ** 0.5

    def __repr__(self):
        return str(self)

    def __str__(self):
        if self.b == 0:
            return str(self.a)
        return "%s%+s*sqrt5" % (self.a, self.b) if self.a else "%s*sqrt5" % self.b


SQRT5 = R5(0, 1)
PHI = R5(Fraction(1, 2), Fraction(1, 2))       # golden ratio (1 + sqrt5)/2


def selftest():
    import random
    rng = random.Random(1)
    assert PHI * PHI == PHI + 1
    assert PHI ** -1 == PHI - 1
    assert SQRT5 * SQRT5 == 5
    for _ in range(20000):
        x = R5(Fraction(rng.randint(-50, 50), rng.randint(1, 9)),
               Fraction(rng.randint(-50, 50), rng.randint(1, 9)))
        y = R5(Fraction(rng.randint(-50, 50), rng.randint(1, 9)),
               Fraction(rng.randint(-50, 50), rng.randint(1, 9)))
        # sign agrees with floating point whenever the float is not tiny
        fx = float(x)
        if abs(fx) > 1e-9:
            assert x.sign() == (fx > 0) - (fx < 0), x
        if x:
            assert x * x.inverse() == 1
        assert (x < y) == (float(x) < float(y)) or abs(float(x) - float(y)) < 1e-9
        assert (x - y) + y == x
    # Pell convergents of sqrt5 miss it by 1/(2 q^2 sqrt5) and alternate in side:
    # 161^2 = 5*72^2 + 1,  682^2 = 5*305^2 - 1,  2889^2 = 5*1292^2 + 1
    assert R5(Fraction(161, 72)) > SQRT5
    assert R5(Fraction(682, 305)) < SQRT5
    assert R5(Fraction(2889, 1292)) > SQRT5
    assert (R5(Fraction(682, 305)) - SQRT5).sign() == -1
    print("qsqrt5 selftest ok")


if __name__ == "__main__":
    selftest()
