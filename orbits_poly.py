from __future__ import print_function
from functools import total_ordering
from sage.all import *
Illegal = Exception("Illegal Operation")
import pickle

# No need to define gcd since it's already a built in function in sage??
# def gcd(x, y):
#     if x == 0:
#         if y == 0:
#             raise ValueError("gcd(0,0) is undefined.")
#         else:
#             return abs(y)
#     x = abs(x)
#     y = abs(y)
#     while y != 0:
#         r = x%y
#         x = y
#         y = r
#     return x

@total_ordering
class Interval:
    """
    A finite subinterval of the integers.
    """
    def __init__(self, a, b):
        self.start = min(a,b)
        self.end = max(a,b)
        self.width = self.end - self.start + 1

    def __repr__(self):
        return '[{}, {}]'.format(self.start, self.end)

    def __eq__(self, other):
        return (self.start, self.end) == (other.start, other.end)

    def __lt__(self, other):
        return (self.start, self.end) < (other.start, other.end)

    def __contains__(self, x):
        """
        True if the Interval contains the integer or Interval argument.
        """
        if isinstance(x, sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular) or isinstance(x, int):
            return self.start <= x <= self.end
        else:
            return self.start <= x.start and x.end <= self.end

    def __xor__(self, other):
        """
        Intersection of two intervals
        """
        start = max(self.start, other.start)
        end = min(self.end, other.end)
        if end < start:
            return None
        else:
            return Interval(start, end)

    def set_end(self, end):
        self.end = end
        self.width = self.end - self.start + 1

    def set_start(self, start):
        self.start = start
        self.width = self.end - self.start + 1

def ToInterval(x):
    """
    Converts an integer or a 2-tuple to an Interval.
    """
    if x.__class__ == Interval:
        return x
    if isinstance(x, sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular) or isinstance(x, int):
        return Interval(x,x)
    else:
        return Interval(x[0],x[1])

class Isometry:
    """
    An element of the infinite dihedral group acting on the integers.
    """
    def __init__(self, shift, flip=0):
        # shift = integer amount(?) rotated, flip = whether or not reflected (either True or 0?)
        self.shift, self.flip = shift, flip

    def __eq__(self, other):
        return self.shift == other.shift and self.flip == other.flip

    def __repr__(self):
        if self.flip:
            return 'x -> -x + {}'.format(self.shift)
        else:
            return 'x -> x + {}'.format(self.shift)

    def __mul__(self, other):
        """
        Composition operator for Isometries.
        """
        flip = self.flip ^ other.flip
        if self.flip:
            shift = self.shift - other.shift
        else:
            shift = other.shift + self.shift
        return Isometry(shift, flip)

    def __pow__(self, n):
        """
        Power operator for Isometries.
        """
        if self.flip:
            if n%2 != 0:
                return Isometry(self.shift, self.flip)
            else:
                return Isometry(0,0)
        else:
            return Isometry(n*self.shift, self.flip)

    def __invert__(self):
        """
        Inversion operator for Isometries.
        """
        if self.flip:
            return Isometry(self.shift, self.flip)
        else:
            return Isometry(-self.shift, self.flip)

    def __call__(self, x):
        """
        An Isometry as a mapping (of an integer or an interval).
        """
        if isinstance(x, sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular) or isinstance(x, int):
            if self.flip:
                return -x + self.shift
            else:
                return x + self.shift
        if x.__class__ == Interval:
            return Interval(self(x.start), self(x.end))

@total_ordering
class Pairing:
    """
    The restriction of an isometry to a finite interval.
    """
    def __init__(self, domain, isometry):
        # domain is an 'Interval' class. isometry is an 'Isometry' class
        self.domain, self.isometry = domain, isometry
        self.range = Interval(self.isometry(domain.start), self.isometry(domain.end))
        self.complexity = (self.range.end, self.domain.width,
                           self.domain.start, self.isometry.flip)

    def __repr__(self):
        # <Used temporarily to print usable text>
        # if self.isometry.flip:
        #     return 'Flip(' + str(self.domain) + ', ' + str(self.range) + ')'
        # else:
        #     return 'Shift(' + str(self.domain) + ', ' + str(self.range) + ')'

        if self.isometry.flip:
            op = ' ~> '
        else:
            op = ' -> '
        return str(self.domain) + op + str(self.range)

    def __call__(self, x):
        """
        A Pairing as a mapping.
        """
        if not x in self.domain:
            raise Illegal("Operand is not contained in domain.")
        else:
            return self.isometry(x)

    def __eq__(self, other):
        return self.domain == other.domain and self.isometry == other.isometry

    def __lt__(self, other):
        # pairings are ordered in reverse(from large to small) by 1. the end of their range
        # 2. width
        # 3. the start of their domain
        if self.range.end != other.range.end:
            return self.range.end > other.range.end
        if self.domain.width != other.domain.width:
            return self.domain.width > other.domain.width
        if self.domain.start != other.domain.start:
            return self.domain.start < other.domain.start
        return self.isometry.flip < other.isometry.flip

    def __contains__(self,x):
        """
        True if the argument is contained in either the domain or range.
        """
        return x in self.domain or x in self.range

    def is_preserving(self):
        """
        True if the Pairing is orientation preserving.
        """
        return self.isometry.flip == 0 or self.domain.width == 1

    def is_periodic(self):
        """
        True if the Pairing is orientation preserving, and
        its domain and range meet.
        """
        return self.is_preserving() and self.domain ^ self.range

    def is_trivial(self):
        """
        True if the Pairing is a restriction of the identity map.
        """
        if self.is_preserving:
            return self.isometry.shift == 0
        else:
            return self.domain.width == 1 and self.domain == self.range

    def contract(self,I):
        """
        Adjust the Pairing to account for removal of a static interval.
        """
        I = ToInterval(I)
        if I ^ self.domain or I ^ self.range:
            raise Illegal("Contraction interval is not static.")
        shift = Isometry(-I.width)
        if I.end < self.domain.start:
            return Pairing(shift(self.domain), shift * self.isometry * ~shift)
        elif I.end < self.range.start:
            return Pairing(self.domain, shift * self.isometry)
        else:
            return self

    def trim(self):
        """
        Trim an orientation reversing pairing so that its domain and
        range become disjoint.
        """
        if self.is_preserving():
            return self
        else:
            intersection = self.domain ^ self.range
            if intersection:
                middle = (self.domain.start + self.range.end - 1) // 2
                domain = Interval(self.domain.start, middle)
                return Pairing(domain, self.isometry)
            else:
                return self

    def merge(self, other):
        """
        Merge a periodic Pairing with an overlapping orientation
        preserving Pairing.
        """
        if self.is_periodic() and other.is_preserving():
            R = Interval(self.domain.start, self.range.end)
            I = R ^ other.domain
            shift = gcd(self.isometry.shift, other.isometry.shift)
            if (other(I) ^ R).width >= self.isometry.shift :
                domain = Interval(R.start, R.end - shift)
                isometry = Isometry(shift)
                return Pairing(domain, isometry)
            else:
                return None
        else:
            raise Illegal("Pairing cannot be merged.")

    def transmit(self, other):
        """
        Left shift the domain and range of another Pairing as far as possible.
        """
        trim = self.trim()
        if other.range not in trim.range:
            return other
        domain = other.domain
        if not trim.is_preserving():
            isometry = trim.isometry * other.isometry
            if domain in trim.range:
                isometry = isometry * trim.isometry
                domain = trim.isometry(domain)
        else:
            shift = trim.isometry.shift
            post = -max(1, (1 + (other.range.start - trim.range.start) // shift))
            isometry = (trim.isometry**post) * other.isometry
            if domain in trim.range:
                pre = max(1, 1 + (other.domain.start - trim.range.start) // shift)
                isometry = isometry * (trim.isometry**pre)
                domain = (trim.isometry**(-pre))(domain)
        range = isometry(domain)
        if range.start < domain.start:
            isometry = isometry**(-1)
            domain = range
        return Pairing(domain, isometry)

def Shift(domain, range):
    """
    Constructor for an orientation preserving pairing, given the domain
    and range.
    """
    if domain.__class__ != Interval:
        domain = Interval(domain[0], domain[1])
    if range.__class__ != Interval:
        range = Interval(range[0], range[1])
    if domain.width != range.width:
        raise Illegal("The domain and range must have the same width.")
    if range.start < domain.start:
        domain, range = range, domain
    isometry = Isometry(range.start - domain.start)
    return Pairing(domain, isometry)

def Flip(domain, range):
    """
    Constructor for an orientation reversing pairing, given the domain
    and range.
    """
    if domain.__class__ != Interval:
        domain = Interval(domain[0], domain[1])
    if range.__class__ != Interval:
        range = Interval(range[0], range[1])
    if domain.width != range.width:
        raise Illegal("The domain and range must have the same width.")
    if range.start < domain.start:
        domain, range = range, domain
    isometry = Isometry(range.end + domain.start, 1)
    return Pairing(domain, isometry)

class Pseudogroup:
    """
    Pseudogroup(P,U) is the pseudogroup of maps of the interval U which
    is generated by the Pairings in the list P.
    """
    def __init__(self, pairings, universe=None):
        self.pairings = pairings
        start = min([p.domain.start for p in self.pairings])
        end = max([p.range.end for p in self.pairings])
        if universe:
            universe = ToInterval(universe)
            if start < universe.start or end > universe.end:
                raise ValueError("Universe must contain all domains and ranges.")
            self.universe = universe
        else:
            self.universe = Interval(start, end)

    def __repr__(self):
        result = 'Pseudogroup on %s:\n'%str(self.universe)
        if self.pairings:
            self.pairings.sort()
            for pairing in self.pairings:
                result += str(pairing) + '\n'
        return result

    def clean(self):
        """
        Get rid of trivial Pairings.
        """
        self.pairings = [p for p in self.pairings if not p.is_trivial()]

    def trim(self):
        """
        Trim all orientation reversing pairings.
        """
        self.pairings = [p.trim() for p in self.pairings]

    def static(self):
        """
        Find a static interval.
        """
        if len(self.pairings) == 0:
            return self.universe
        intervals = [p.domain for p in self.pairings]
        intervals += [p.range for p in self.pairings]
        intervals.sort()
        I = intervals.pop(0)
        start, end = I.start, I.end
        if start > 1:
            return Interval(1, start - 1)
        for interval in intervals:
            if end < interval.start - 1:
                return Interval(end+1, interval.start-1)
            end = max(end,interval.end)
        if end < self.universe.end:
            return Interval(end + 1, self.universe.end)
        return None

    def contract(self):
        """
        Remove all static intervals.  Return the total size.
        """
        result = 0
        I = self.static()
        while I:
            result += I.width
            if I.end != self.universe.end:
                self.pairings =  [p.contract(I) for p in self.pairings]
            self.universe.set_end(self.universe.end - I.width)
            I = self.static()
        return result

    def merge(self):
        """
        Merge periodic pairing whenever possible.
        """
        if len(self.pairings) < 2:
             return
        done=0
        while not done:
            periodics = [p for p in self.pairings if p.is_periodic()]
            done=1
            for p in periodics[:-1]:
                for q in periodics[1+periodics.index(p):]:
                    g = None
                    try:
                        g = p.merge(q)
                    except: pass
                    if g:
                        self.pairings.remove(p)
                        self.pairings.remove(q)
                        self.pairings.append(g)
                        done=0
                        break

    def transmit(self):
        """
        Use the largest Pairing to transmit others.
        """
        self.pairings.sort()
        g = self.pairings[0]
        self.pairings = [g] + [ g.transmit(p) for p in self.pairings[1:] ]

    def truncate(self):
        """
        Truncate the largest pairing.
        """
        self.pairings.sort()
        g = self.pairings.pop(0)
        if len(self.pairings) > 0:
            support_end = self.pairings[0].range.end
        else:
            support_end = g.range.start - 1
        if support_end < g.range.start:
            self.universe.set_end(support_end)
            return
        if not g.is_preserving():
            g.trim()
        self.universe.set_end(support_end)
        range = Interval(g.range.start, support_end)
        domain = (~g.isometry)(range)
        self.pairings = [Pairing(domain, g.isometry)] + self.pairings

    def simplify(self):
        """
        Do one cycle of the orbit counting reduction algorithm due
        to Agol, Hass and Thurston.
        """
        self.clean()
        if len(self.pairings) == 0:
            self.pairings = None
            return self.universe.width
        print("cleaned\n", self)
        count = self.contract()
        print("contracted\n", self)
        self.trim()
        print("trimmed\n", self)
        self.merge()
        print("merged\n", self)
        self.transmit()
        print("transmitted\n", self)
        self.truncate()
        print("truncated\n", self)
        print('count = ', count)
        return count

    def reduce(self):
        """
        Reduce the pseudogroup to nothing. Return the number of orbits.
        """
        count = 0
        while self.pairings and len(self.pairings) != 0:
            count += self.simplify()
        return count



def main():
    # f = open('example.pickle', 'rb')
    # test = pickle.load(f)
    # test_interval = Interval(test[0].start, test[0].end)
    # test_pairings = test[1]

    x0, x1 = PolynomialRing(ZZ, 2, 'x').gens()
    # test_interval = [1, 24 * x0 + 24 * x1]
    # test_pairings = [Shift([1, 2 * x1], [2 * x1 + 1, 4 * x1]),
    #                  Flip([1, 2 * x1], [6 * x0 + 8 * x1 + 1, 6 * x0 + 10 * x1]),
    #                  Shift([x0 + 3 * x1 + 1, 2 * x0 + 6 * x1], [5 * x0 + 7 * x1 + 1, 6 * x0 + 10 * x1]),
    #                  Flip([2 * x1 + 1, x0 + 3 * x1], [11 * x0 + 13 * x1 + 1, 12 * x0 + 14 * x1]),
    #                  Shift([2 * x0 + 6 * x1 + 1, 5 * x0 + 7 * x1], [8 * x0 + 12 * x1 + 1, 11 * x0 + 13 * x1]),
    #                  Shift([6 * x0 + 10 * x1 + 1, 7 * x0 + 11 * x1], [8 * x0 + 12 * x1 + 1, 9 * x0 + 13 * x1]),
    #                  Shift([7 * x0 + 11 * x1 + 1, 8 * x0 + 12 * x1], [11 * x0 + 13 * x1 + 1, 12 * x0 + 14 * x1]),
    #                  Flip([8 * x0 + 12 * x1 + 1, 11 * x0 + 13 * x1], [9 * x0 + 13 * x1 + 1, 12 * x0 + 14 * x1]),
    #                  Shift([12 * x0 + 14 * x1 + 1, 13 * x0 + 16 * x1], [14 * x0 + 16 * x1 + 1, 15 * x0 + 18 * x1]),
    #                  Shift([11 * x0 + 14 * x1 + 1, 12 * x0 + 14 * x1], [13 * x0 + 16 * x1 + 1, 14 * x0 + 16 * x1]),
    #                  Flip([8 * x0 + 12 * x1 + 1, 11 * x0 + 14 * x1], [15 * x0 + 18 * x1 + 1, 18 * x0 + 20 * x1]),
    #                  Shift([2 * x1 + 1, x0 + 4 * x1], [12 * x0 + 14 * x1 + 1, 13 * x0 + 16 * x1]),
    #                  Flip([6 * x0 + 10 * x1 + 1, 7 * x0 + 10 * x1], [13 * x0 + 16 * x1 + 1, 14 * x0 + 16 * x1]),
    #                  Shift([x0 + 4 * x1 + 1, 2 * x0 + 6 * x1], [7 * x0 + 10 * x1 + 1, 8 * x0 + 12 * x1]),
    #                  Flip([2 * x1 + 1, x0 + 5 * x1], [17 * x0 + 17 * x1 + 1, 18 * x0 + 20 * x1]),
    #                  Flip([9 * x0 + 13 * x1 + 1, 12 * x0 + 14 * x1], [14 * x0 + 16 * x1 + 1, 17 * x0 + 17 * x1]),
    #                  Flip([x0 + 5 * x1 + 1, 2 * x0 + 6 * x1], [8 * x0 + 12 * x1 + 1, 9 * x0 + 13 * x1]),
    #                  Shift([2 * x1 + 1, 2 * x0 + 6 * x1], [4 * x0 + 6 * x1 + 1, 6 * x0 + 10 * x1]),
    #                  Flip([2 * x1 + 1, x0 + 5 * x1], [x0 + 3 * x1 + 1, 2 * x0 + 6 * x1]),
    #                  Shift([2 * x1 + 1, x0 + 3 * x1], [20 * x0 + 20 * x1 + 1, 21 * x0 + 21 * x1]),
    #                  Shift([x0 + 5 * x1 + 1, 2 * x0 + 6 * x1], [21 * x0 + 21 * x1 + 1, 22 * x0 + 22 * x1]),
    #                  Flip([1, 2 * x1], [14 * x0 + 14 * x1 + 1, 14 * x0 + 16 * x1]),
    #                  Flip([2 * x0 + 6 * x1 + 1, 4 * x0 + 8 * x1], [12 * x0 + 14 * x1 + 1, 14 * x0 + 16 * x1]),
    #                  Shift([4 * x0 + 8 * x1 + 1, 6 * x0 + 10 * x1], [22 * x0 + 22 * x1 + 1, 24 * x0 + 24 * x1]),
    #                  Flip([10 * x0 + 14 * x1 + 1, 12 * x0 + 14 * x1], [18 * x0 + 20 * x1 + 1, 20 * x0 + 20 * x1]),
    #                  Flip([8 * x0 + 12 * x1 + 1, 10 * x0 + 14 * x1], [22 * x0 + 22 * x1 + 1, 24 * x0 + 24 * x1]),
    #                  Shift([1, x1], [22 * x0 + 22 * x1 + 1, 22 * x0 + 23 * x1]),
    #                  Flip([x1 + 1, 2 * x1], [20 * x0 + 20 * x1 + 1, 20 * x0 + 21 * x1]),
    #                  Shift([20 * x0 + 21 * x1 + 1, 22 * x0 + 22 * x1], [22 * x0 + 23 * x1 + 1, 24 * x0 + 24 * x1]),
    #                  Flip([2 * x1 + 1, 3 * x1], [24 * x0 + 23 * x1 + 1, 24 * x0 + 24 * x1]),
    #                  Flip([16 * x0 + 19 * x1 + 1, 18 * x0 + 20 * x1], [22 * x0 + 22 * x1 + 1, 24 * x0 + 23 * x1]),
    #                  Flip([3 * x1 + 1, 2 * x0 + 6 * x1], [14 * x0 + 16 * x1 + 1, 16 * x0 + 19 * x1]),
    #                  Shift([20 * x0 + 20 * x1 + 1, 20 * x0 + 21 * x1], [22 * x0 + 22 * x1 + 1, 22 * x0 + 23 * x1]),
    #                  Shift([10 * x0 + 13 * x1 + 1, 12 * x0 + 14 * x1], [20 * x0 + 21 * x1 + 1, 22 * x0 + 22 * x1]),
    #                  Flip([8 * x0 + 12 * x1 + 1, 10 * x0 + 13 * x1], [22 * x0 + 23 * x1 + 1, 24 * x0 + 24 * x1]),
    #                  Flip([7 * x0 + 10 * x1 + 1, 8 * x0 + 12 * x1], [12 * x0 + 14 * x1 + 1, 13 * x0 + 16 * x1]),
    #                  Shift([13 * x0 + 16 * x1 + 1, 14 * x0 + 16 * x1], [19 * x0 + 20 * x1 + 1, 20 * x0 + 20 * x1]),
    #                  Shift([6 * x0 + 10 * x1 + 1, 7 * x0 + 10 * x1], [18 * x0 + 20 * x1 + 1, 19 * x0 + 20 * x1]),
    #                  Shift([1, 2 * x1], [2 * x0 + 4 * x1 + 1, 2 * x0 + 6 * x1]),
    #                  Shift([2 * x1 + 1, 2 * x0 + 4 * x1], [12 * x0 + 14 * x1 + 1, 14 * x0 + 16 * x1])]

    # test_interval = [1, 2*x0 + 2*x1]
    # test_pairings = [Shift([x1 + 1, 2*x0], [3*x1 + 1, 2*x0 + 2*x1]),
    #                  Flip([1, x0 + x1], [x0 + x1 + 1, 2*x0 + 2*x1])]

    test_interval = [1, 2 * x0]
    test_pairings = [Flip([1, x0], [x0 - 2 * x1 + 1, 2 * x0 - 2 * x1]),
                     Shift([x1 + 1, 2 * x0 - 2 * x1], [3 * x1 + 1, 2 * x0])]


    # test_PG = Pseudogroup(test_pairings, test_interval)
    # print(test_PG.reduce())

if __name__ == '__main__':
    main()

# probably the last pairing Flip([9 * x0 + 13 * x1 + 1, 12 * x0 + 14 * x1], [14 * x0 + 16 * x1 + 1, 17 * x0 + 17 * x1]) is causing some issue
# the code runs forever for some reason when this pairing is used
# the problem is reduced to
# interval: [1, 11*x0 + 14*x1]
# remaining pairings:
# [9*x0 + 13*x1 + 1, 11*x0 + 12*x1] -> [9*x0 + 15*x1 + 1, 11*x0 + 14*x1]
# [9*x0 + 12*x1 + 1, 10*x0 + 13*x1] ~> [10*x0 + 13*x1 + 1, 11*x0 + 14*x1]
# [2*x0 + 6*x1 + 1, 5*x0 + 7*x1] -> [8*x0 + 12*x1 + 1, 11*x0 + 13*x1]
# [2*x1 + 1, x0 + 5*x1] -> [8*x0 + 12*x1 + 1, 9*x0 + 15*x1]
# [2*x1 + 1, x0 + 4*x1] -> [8*x0 + 12*x1 + 1, 9*x0 + 14*x1]
# [2*x1 + 1, x0 + 3*x1] -> [8*x0 + 12*x1 + 1, 9*x0 + 13*x1]
# [6*x0 + 10*x1 + 1, 7*x0 + 11*x1] -> [8*x0 + 12*x1 + 1, 9*x0 + 13*x1]
# [7*x0 + 11*x1 + 1, 8*x0 + 12*x1] ~> [8*x0 + 12*x1 + 1, 9*x0 + 13*x1]
# [6*x0 + 10*x1 + 1, 7*x0 + 10*x1] -> [8*x0 + 12*x1 + 1, 9*x0 + 12*x1]
# [x0 + 4*x1 + 1, 2*x0 + 6*x1] -> [7*x0 + 10*x1 + 1, 8*x0 + 12*x1]
# [x0 + 3*x1 + 1, 2*x0 + 6*x1] -> [5*x0 + 7*x1 + 1, 6*x0 + 10*x1]
# [1, 2*x1] ~> [6*x0 + 8*x1 + 1, 6*x0 + 10*x1]
# [1, 2*x1] -> [2*x1 + 1, 4*x1]