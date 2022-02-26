class Bits:
    def __init__(self, bits: str = '1001') -> None:
        self._bits = bits.replace(' ', '').replace('_', '')

    def to_int(self) -> int:
        return self.int(self._bits)

    def to_f32(self, verbose: bool = True) -> float:
        assert len(self._bits) == 32
        sign, exponent, fraction = self._bits[0], self._bits[1:9], self._bits[9:]
        if all(e=='0' for e in exponent):
            if all(f=='0' for f in fraction):
                ans = +0.0 if sign=='0' else -0.0
            else:
                ans = (1.0 if sign=='0' else -1.0) * 2**-126 * self.float('0.'+fraction)
            if verbose:
                print(f'  (-1)^{sign} x 2^-126 x (0.{fraction})₂')
                sign, fraction = self.int(sign), self.float('0.'+fraction)
                print(f'= (-1)^{sign} x 2^-126 x {fraction}')
                print(f'= {(-1)**sign} x {2**-126} x {fraction}')
                print(f'= {ans}')
        elif all(e=='1' for e in exponent):
            if all(f=='0' for f in fraction):
                ans = float('+inf' if sign=='0' else '-inf')
            else:
                ans = float('nan')
            if verbose:
                print(f'= {ans}')
        else:
            ans = (1.0 if sign=='0' else -1.0) * 2**(self.int(exponent)-127) * self.float('1.'+fraction)
            if verbose:
                print(f'  (-1)^{sign} x 2^(({exponent})₂-127) x (1.{fraction})₂')
                sign, exponent, fraction = self.int(sign), self.int(exponent), self.float('1.'+fraction)
                print(f'= (-1)^{sign} x 2^({exponent}-127) x {fraction}')
                print(f'= {(-1)**sign} x {2**(exponent-127)} x {fraction}')
                print(f'= {ans}')
        return ans

    @classmethod
    def int(cls, bits: str) -> int:
        return int(bits, base=2)

    @classmethod
    def float(cls, bits: str) -> float:
        integer, fraction = bits.split('.')
        ans, cof = float(cls.int(integer)), 1.0
        for bit in fraction:
            cof /= 2.0
            if bit == '1':
                ans += cof
        return ans


if __name__ == '__main__':
    import math

    api = lambda bits: Bits(bits).to_f32(verbose=False)
    err = lambda x, y: abs(x-y) / y
    eps = 1e-10

    # Examples from https://en.wikipedia.org/wiki/Single-precision_floating-point_format
    assert math.isnan(api('0 11111111 00000000000000000000001'))
    assert api('0 11111111 00000000000000000000000') == float('+inf')
    assert api('1 11111111 00000000000000000000000') == float('-inf')
    assert api('0 01111100 01000000000000000000000') == 0.15625
    assert api('0 10000010 10001100000000000000000') == 12.375
    assert api('0 01111111 00000000000000000000000') == 1.0
    assert api('1 10000000 00000000000000000000000') == -2.0
    assert err(api('0 00000000 11111111111111111111111'), 1.1754942107e-38) <= eps
    assert err(api('0 00000001 00000000000000000000000'), 1.1754943508e-38) <= eps
    assert err(api('0 11111110 11111111111111111111111'), 3.4028234664e+38) <= eps
    assert err(api('0 01111110 11111111111111111111111'), 0.999999940395355225) <= eps
    assert err(api('0 01111111 00000000000000000000001'), 1.00000011920928955) <= eps
