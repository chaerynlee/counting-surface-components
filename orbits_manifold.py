from __future__ import print_function
from functools import total_ordering
from sage.all import *
import itertools
from Orbits import orbits

Illegal = Exception("Illegal Operation")

class Polynomial:
    def __init__(self, constant=0, **poly):
        self.poly = poly  # dictionary with keys: variables, values: coefficient
        self.num_var = len(poly)
        self.constant = constant

    def __repr__(self):
        if self.poly == dict():
            return str(self.constant)

        name = ''
        for var in self.poly.keys():
            if self.poly[var] == 0:
                continue
            elif self.poly[var] == 1:
                name += str(var) + ' + '
            elif self.poly[var] == -1:
                name += '-' + str(var) + ' + '
            else:
                name += str(self.poly[var]) + str(var) + ' + '
        name += str(self.constant)
        return name

    def __eq__(self, other):
        if isinstance(other, int):
            other = Polynomial(other)
        new_self_poly = dict()
        for var in self.poly.keys():
            new_self_poly[var] = self.poly[var]
        for var in set(other.poly.keys()) - set(self.poly.keys()):
            new_self_poly[var] = 0

        new_other_poly = dict()
        for var in other.poly.keys():
            new_other_poly[var] = other.poly[var]
        for var in set(self.poly.keys()) - set(other.poly.keys()):
            new_other_poly[var] = 0

        new_self = Polynomial(self.constant, **new_self_poly)
        new_other = Polynomial(other.constant, **new_other_poly)
        return (new_self.poly, new_self.constant) == (new_other.poly, new_other.constant)

    def __neg__(self):
        neg_constant = - self.constant
        neg_poly = dict()
        for var in self.poly.keys():
            neg_poly[var] = -self.poly[var]
        return Polynomial(neg_constant, **neg_poly)

    def __add__(self, other):
        if isinstance(other, int):
            other = Polynomial(other)
        sum_constant = self.constant + other.constant
        sum_poly = dict()
        for var in self.poly:
            sum_poly[var] = self.poly[var]
        for var in other.poly:
            if var in self.poly:
                sum_poly[var] = self.poly[var] + other.poly[var]
            else:
                sum_poly[var] = other.poly[var]
        return Polynomial(sum_constant, **sum_poly)

    def __sub__(self, other):
        if isinstance(other, int):
            other = Polynomial(other)
        sub_constant = self.constant - other.constant
        sub_poly = dict()
        for var in self.poly:
            sub_poly[var] = self.poly[var]
        for var in other.poly:
            if var in self.poly:
                sub_poly[var] = self.poly[var] - other.poly[var]
            else:
                sub_poly[var] = -other.poly[var]
        return Polynomial(sub_constant, **sub_poly)

    # def __lt__(self, other):
    #     # IMPORTANT: this polynomial is a weak inequality by design
    #     if set(self.poly.keys()) == set(other.poly.keys()):
    #         # case where polynomials have same variables
    #         comparison = []
    #         if self.constant == other.constant:
    #             comparison.append('S')
    #         else:
    #             comparison.append(self.constant < other.constant)
    #         for var in self.poly.keys():
    #             if self.poly[var] == other.poly[var]:
    #                 comparison.append('S')
    #             else:
    #                 comparison.append(self.poly[var] < other.poly[var])
    #         comparison = [e for e in comparison if e != 'S']
    #         if comparison == None:
    #             return 'equal'
    #         elif set(comparison) == {True}:
    #             return True
    #         elif set(comparison) == {False}:
    #             return False
    #         else:
    #             return 'unknown'
    #
    #     elif set(self.poly.keys()) < set(other.poly.keys()):
    #         # case where the other polynomial has more variables than the given polynomial, can only return True or unknown
    #         comparison = []
    #         comparison_ne = []
    #         nonexistent_var = set(other.poly.keys()).difference(self.poly.keys())
    #         comparison.append(self.constant <= other.constant)
    #         for var in self.poly.keys():
    #             comparison.append(self.poly[var] <= other.poly[var])
    #
    #         for var in nonexistent_var:
    #             comparison_ne.append(other.poly[var] >= 0)
    #
    #         if False not in comparison and False not in comparison_ne:
    #             return True
    #         else:
    #             return 'unknown'
    #
    #     elif set(other.poly.keys()) < set(self.poly.keys()):
    #         # case where the our polynomial has more variables than the other polynomial, can only return False or unknown
    #         comparison = []
    #         comparison_ne = []
    #         nonexistent_var = set(self.poly.keys()).difference(other.poly.keys())
    #         comparison.append(other.constant <= self.constant)
    #         for var in other.poly.keys():
    #             comparison.append(other.poly[var] <= self.poly[var])
    #
    #         for var in nonexistent_var:
    #             comparison_ne.append(self.poly[var] >= 0)
    #         if False not in comparison and False not in comparison_ne:
    #             return False
    #         else:
    #             return 'unknown'
    #
    #     else:
    #         # when variables differ we cannot make a comparison
    #         return 'unknown'

    def __lt__(self, other):
        if isinstance(other, int):
            other = Polynomial(other)
        new_self_poly = dict()
        for var in self.poly.keys():
            new_self_poly[var] = self.poly[var]
        for var in set(other.poly.keys()) - set(self.poly.keys()):
            new_self_poly[var] = 0

        new_other_poly = dict()
        for var in other.poly.keys():
            new_other_poly[var] = other.poly[var]
        for var in set(self.poly.keys()) - set(other.poly.keys()):
            new_other_poly[var] = 0

        new_self = Polynomial(self.constant, **new_self_poly)
        new_other = Polynomial(other.constant, **new_other_poly)

        comparison = []
        if new_self.constant == new_other.constant:
            comparison.append('S')
        else:
            comparison.append(new_self.constant < new_other.constant)
        for var in new_self.poly.keys():
            if new_self.poly[var] == new_other.poly[var]:
                comparison.append('S')
            else:
                comparison.append(new_self.poly[var] < new_other.poly[var])
        comparison = [e for e in comparison if e != 'S']
        if comparison == []:
            return False
        elif set(comparison) == {True}:
            return True
        elif set(comparison) == {False}:
            return False
        else:
            return 'unknown'

    def __gt__(self, other):
        if isinstance(other, int):
            other = Polynomial(other)
        return other < self

    def __mul__(self, other):
        if not isinstance(other, Polynomial):
            other = Polynomial(constant=other, **{})
        # for now assume that we only need multiplication of integers not of polynomials themselves
        if other.poly != dict():
            raise Exception('Can only multiply polynomials with numbers')
        mul_constant = self.constant * other.constant
        mul_poly = dict()
        for var in self.poly.keys():
            mul_poly[var] = self.poly[var] * other.constant
        return Polynomial(mul_constant, **mul_poly)

    def __truediv__(self, other):
        # for now assume that we only need division of integers not of polynomials themselves
        if other.poly != dict():
            raise Exception('Can only multiply polynomials with numbers')
        div_constant = self.constant * other.constant
        div_poly = dict()
        for var in self.poly.keys():
            div_poly[var] = self.poly[var] / other.constant
        return Polynomial(div_constant, **div_poly)

    def evaluate(self, **kwargs):
        assert set(self.poly.keys()).issubset(set(kwargs.keys()))
        num = self.constant
        for var in self.poly.keys():
            num += kwargs[var] * self.poly[var]
        return num


def linear_transform_poly(poly, transform):
    """
    For now, for simplicity's sake poly must not have a constant()
    transform: list of dictionaries with the same format as polynomials, note that the keys in transform must
    match that of poly
    eg. [{'y0':1, 'y1':2}, {'y0':3, 'y1':-1}] means x0 = 1y0 + 2y1, x1 = 3y0 - 1y1
    """
    if poly == Polynomial(0):
        return Polynomial(0)
    assert poly.constant == 0
    trans_poly = Polynomial(0)
    var = list(poly.poly.keys())
    var.sort()
    for i, eq in enumerate(transform):
        coeff = poly.poly[var[i]]
        trans_poly += Polynomial(**eq) * coeff
    return trans_poly

def linear_transform_pairing(pairing, transform):
    domain = Interval(linear_transform_poly(pairing.domain.start, transform), linear_transform_poly(pairing.domain.end, transform))
    iso = Isometry(linear_transform_poly(pairing.isometry.shift, transform), pairing.isometry.flip)
    return Pairing(domain, iso, pairing.domain_index, pairing.range_index)

def gcd(x, y):
    if x == 0:
        if y == 0:
            raise ValueError("gcd(0,0) is undefined.")
        else:
            return abs(y)
    x = abs(x)
    y = abs(y)
    while y != 0:
        r = x%y
        x = y
        y = r
    return x


class Interval:
    """
    A finite subinterval of the integers.
    """
    def __init__(self, a, b):
        if (a < b) == 'unknown':
            self.start = a
            self.end = b
        else:
            self.start = min(a, b)
            self.end = max(a, b)
        self.width = self.end - self.start
        # self.width = Polynomial(1) + self.end - self.start

    def __repr__(self):
        return '[{}, {}]'.format(self.start, self.end)

    def __eq__(self, other):
        if not isinstance(other, Interval):
            return False
        if self.start == other.start and self.end == other.end:
            return True
        else:
            return False
    def __lt__(self, other):
        if (self.start).__lt__(other.start) == 'unknown':
            return 'unknown'
        if self.start == other.start:
            return (self.end).__lt__(other.end)
        if (self.start).__lt__(other.start) == True:
            return True
        else:
            return False

    def __contains__(self, x):
        """
        True if the Interval contains the integer or Interval argument.
        """
        if isinstance(x, int):
            x = Polynomial(x)

        if isinstance(x, Polynomial):
            if x == self.start or x == self.end:
                return True
            comp_start = (self.start < x)
            comp_end = (x < self.end)
            if comp_start == 'unknown' or comp_end == 'unknown':
                return 'unknown'
            elif comp_start == True and comp_end == True:
                return True
            else:
                return False
        elif isinstance(x, Interval):
            comp_start = (self.start < x.start)
            comp_end = (x.end < self.end)
            if comp_start == 'unknown' or comp_end == 'unknown':
                return 'unknown'
            elif x.start == self.start and x.end == self.end:
                return True
            elif x.start == self.start and comp_end == True:
                return True
            elif comp_start == True and x.end == self.end:
                return True
            elif comp_start == True and comp_end == True:
                return True
            else:
                return False
        else:
            raise TypeError(f'{x} is not a Polynomial or Interval')

    def __xor__(self, other):
        """
        Intersection of two intervals
        """
        if self.start.__lt__(other.start) == 'unknown':
            return 'unknown'
        else:
            start = max(self.start, other.start)

        if self.end.__lt__(other.start) == 'unknown':
            return 'unknown'
        else:
            end = min(self.end, other.end)

        if start.__lt__(end) == 'unknown':
            return 'unknown'
        elif start.__lt__(end) == True:
            return Interval(start, end)
        else:
            return None

    def set_end(self, end):
        self.end = end
        self.width = self.end - self.start

    def set_start(self, start):
        self.start = start
        self.width = self.end - self.start

    def union(self, other):
        """
        Attempted to use this for find_significant_edges but turns out to be unhelpful.

        Takes the union of self with other.
        Can only be performed in the case where self.start is smaller than other.start and where the two intervals intersect.
        Returns False if the union cannot be performed.
        """
        if (self.start).__lt__(other.start) == True or self.start == other.start:
            start = self.start
        else:
            return False

        if (self.end).__lt__(other.start) == True or (self.end).__lt__(other.start) == 'unknown':
            return False

        if (self.end).__lt__(other.end) == True or self.end == other.end:
            end = other.end
            return Interval(start, end)
        else:
            end = self.end
            return Interval(start, end)

def ToInterval(x):
    """
    Converts an integer or a 2-tuple to an Interval.
    """
    if x.__class__ == Interval:
        return x
    if isinstance(x, Polynomial) or isinstance(x, int):
        return Interval(x,x)
    if x == []:
        return x
    else:
        return Interval(x[0],x[1])

class Isometry:
    """
    An element of the infinite dihedral group acting on the integers.
    """
    def __init__(self, shift, flip=0):
        # shift = integer amount rotated, flip = whether or not reflected written in binary (0: False, 1: True)
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
        if isinstance(x, Polynomial) or isinstance(x, int):
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
    def __init__(self, domain, isometry, domain_index, range_index):
        # domain is an 'Interval' class. isometry is an 'Isometry' class
        self.domain, self.isometry = domain, isometry
        self.range = Interval(self.isometry(domain.start), self.isometry(domain.end))
        self.complexity = (self.range.end, self.domain.width,
                           self.domain.start, self.isometry.flip)
        self.domain_index = domain_index
        self.range_index = range_index

    def __repr__(self):
        if self.isometry.flip:
            op = ' ~> '
        else:
            op = ' -> '
        return f'{self.domain}{op}{self.range}, {self.domain_index}{op}{self.range_index}, {self.isometry}'

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
        # pairings are ordered in reverse (from large to small) by
        # 1. the end of their range
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
        return self.isometry.flip == 0

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
        if self.domain.width == Polynomial(0):
            return True
        if self.is_preserving():
            return self.isometry.shift == 0
        else:
            return self.domain.width == Polynomial(0) and self.domain == self.range

    def contract(self,I):
        """
        Adjust the Pairing to account for removal of a static interval.
        Any returned pairing has domain/range_index set to 0. Only used in Pseudogroup_comparable and not used in Pseudogroup
        """
        I = ToInterval(I)
        if I.__xor__(self.domain) == 'unknown' or I.__xor__(self.range) == 'unknown':
            raise ValueError(f'Unknown whether the interval can be contracted. {self}, {I}')
        if I ^ self.domain or I ^ self.range:
            raise ValueError(f'Contraction interval is not static. {self}, {I}')

        shift = Isometry(-I.width)
        if I.end < self.domain.start or I.end == self.domain.start:
            if I.end < self.range.end or I.end == self.range.start:
                return Pairing(shift(self.domain), shift * self.isometry * ~shift, 0, 0)
            elif I.start > self.range.end or I.start == self.range.end:
                return Pairing(shift(self.domain), self.isometry * ~shift, 0, 0)
        elif I.end < self.range.start or I.end == self.range.start:
            return Pairing(self.domain, shift * self.isometry, 0, 0)
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
            if self.domain.start == self.range.start:
                intersection = True
            domain_left = self.domain.start.__lt__(self.range.start)
            if domain_left == True:
                intersection = self.range.start.__lt__(self.domain.end)
            elif domain_left == False:
                intersection = self.domain.start.__lt__(self.range.end)
            elif domain_left == 'unknown':
                return self

            if intersection == True:
                if domain_left == True:
                    middle = (self.domain.start + self.range.end) / Polynomial(2)
                    domain = Interval(self.domain.start, middle)
                    # print('trimmed', self.__repr__(), 'to', Pairing(domain, self.isometry, self.domain_index, self.range_index))
                    return Pairing(domain, self.isometry, self.domain_index, self.range_index)
                elif domain_left == False:
                    middle = (self.domain.end + self.range.start) / Polynomial(2)
                    domain = Interval(middle, self.domain.end)
                    # print('trimmed', self.__repr__(), 'to', Pairing(domain, self.isometry, self.domain_index, self.range_index))
                    return Pairing(domain, self.isometry, self.domain_index, self.range_index)
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
        Just shift other by self once.
        Any returned pairing has domain/range_index set to 0.
        Only used in Pseudogroup_comparable and not used in Pseudogroup.
        """
        trim = self.trim()
        range_inclusion = (trim.range).__contains__(other.range)
        domain_inclusion = (trim.range).__contains__(other.domain)
        if range_inclusion in [False, 'unknown']:
            # print(self)
            # print(other)
            # print('not transmitted')
            # print()
            return other
        if domain_inclusion == True:
            inv_trim = ~trim.isometry
            domain = inv_trim(other.domain)
            isometry = ~trim.isometry * other.isometry * trim.isometry
        else:
            domain = other.domain
            isometry = ~trim.isometry * other.isometry
        range = isometry(domain)
        p = Pairing(domain, isometry, 0, 0)
        # print(self)
        # print(other)
        # print(p)
        # print()
        return p.switch_domain_range()

    def transmit_inclusive(self, other):
        """
        Transmit other pairing by self if either the domain or range lies in the domain or range of self
        """
        if other.domain_index == other.range_index:
            if self.domain.__contains__(other.domain)  == True:
                new_domain = self.isometry(other.domain)
                new_iso = self.isometry * other.isometry * ~self.isometry
                return Pairing(new_domain, new_iso, self.range_index, self.range_index)
            elif self.range.__contains__(other.domain) == True:
                new_domain = (~self.isometry)(other.domain)
                new_iso = ~self.isometry * other.isometry * self.isometry
                return Pairing(new_domain, new_iso, self.domain_index, self.domain_index)
            else:
                raise Illegal('This pairing cannot be transmitted')

        if self.domain.__contains__(other.domain) == True:
            assert self.domain_index == other.domain_index
            new_iso = self.isometry * ~other.isometry
            return Pairing(other.range, new_iso, other.range_index, self.range_index)
        elif self.domain.__contains__(other.range) == True:
            assert self.domain_index == other.range_index
            new_iso = self.isometry * other.isometry
            return Pairing(other.domain, new_iso, other.domain_index, self.range_index)
        elif self.range.__contains__(other.domain) == True:
            assert self.range_index == other.domain_index
            new_iso = ~self.isometry * ~other.isometry
            return Pairing(other.range, new_iso, other.range_index, self.domain_index)
        elif self.range.__contains__(other.range) == True:
            assert self.range_index == other.range_index
            new_iso = ~self.isometry * other.isometry
            return Pairing(other.domain, new_iso, other.domain_index, self.domain_index)
        else:
            raise Illegal('This pairing cannot be transmitted')

    def switch_domain_range(self):
        if self.domain_index > self.range_index:
            return Pairing(self.range, ~self.isometry, self.range_index, self.domain_index)
        elif self.domain_index == self.range_index:
            if (self.range.start).__lt__(self.domain.start) == True:
                return Pairing(self.range, ~self.isometry, self.range_index, self.domain_index)
        return self


def Shift(domain, range, domain_index, range_index):
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
        domain_index, range_index = range_index, domain_index
    isometry = Isometry(range.start - domain.start)
    return Pairing(domain, isometry, domain_index, range_index)

def Flip(domain, range, domain_index, range_index):
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
        domain_index, range_index = range_index, domain_index
    isometry = Isometry(range.end + domain.start, 1)
    return Pairing(domain, isometry, domain_index, range_index)

class Pseudogroup:
    """
    Pseudogroup(P,U) is the pseudogroup of maps of the interval U which
    is generated by the Pairings in the list P.
    """
    def __init__(self, pairings, universe, divided_universe):
        self.pairings = pairings
        if isinstance(universe, Interval):
            self.universe = universe
        else:
            self.universe = ToInterval(universe)
        self.divided_universe = [ToInterval(I) for I in divided_universe]

    def __repr__(self):
        result = 'Pseudogroup on %s:\n'%str(self.universe)
        if self.pairings:
            self.pairings.sort()
            for pairing in self.pairings:
                result += str(pairing) + '\n'
        return result

    def clean(self):
        """
        Make sure domain_index < range_index for all pairings,
        for pairings whose domain_index = range_index make sure the domain is to the left of the range,
        and remove duplicate pairings and trivial pairings.
        """
        self.pairings = [p.switch_domain_range() for p in self.pairings if not p.is_trivial()]
        remove_dup = []
        [remove_dup.append(p) for p in self.pairings if p not in remove_dup]
        self.pairings = remove_dup

    def trim(self):
        """
        Trim all orientation reversing pairings.
        """
        self.pairings = [p.trim() for p in self.pairings]

    def identify_edge_containing(self, interval):
        """
        For the given interval, find a pairing that has that as its domain or range.
        If no such pairing exists return False.
        """
        if not isinstance(interval, Interval):
            return False
        for pairing in self.pairings:
            if pairing.domain == interval and pairing.domain_index != pairing.range_index:
                return pairing
            elif pairing.range == interval and pairing.domain_index != pairing.range_index:
                return pairing
        return False

    def peel_edge(self, edge_index):
        """
        Removes subinterval of given edge index from divided universe (makes it into an empty interval,
        the corresponding edge index is maintained)
        Translates parings to accommodate for this removal.
        """
        interval = self.divided_universe[edge_index]

        shifted_pairings = []
        for p in self.pairings:
            # print('original pairing', p)
            assert p.domain_index <= p.range_index
            if edge_index < p.domain_index:
                new_domain = Interval(p.domain.start - interval.width, p.domain.end - interval.width)
                # possible fix?
                if p.isometry.flip:
                    new_iso = p.isometry * Isometry(interval.width * Polynomial(2), 0)
                    # print('original iso', p.isometry)
                    # print('changed to', Pairing(new_domain, new_iso, p.domain_index, p.range_index))
                    shifted_pairings.append(Pairing(new_domain, new_iso, p.domain_index, p.range_index))
                else:
                    shifted_pairings.append(Pairing(new_domain, p.isometry, p.domain_index, p.range_index))
                # shifted_pairings.append(Pairing(new_domain, p.isometry, p.domain_index, p.range_index))
            elif p.domain_index < edge_index < p.range_index:
                new_iso = Isometry(-interval.width, 0) * p.isometry
                # print('new iso', new_iso)
                # print('changed to', Pairing(p.domain, new_iso, p.domain_index, p.range_index))
                shifted_pairings.append(Pairing(p.domain, new_iso, p.domain_index, p.range_index))
            else:
                shifted_pairings.append(p)
        self.pairings = shifted_pairings

        new_end = self.universe.end - interval.width
        self.universe = Interval(Polynomial(0), new_end)

        for i in range(len(self.divided_universe)):
            if i < edge_index:
                continue
            elif i == edge_index:
                self.divided_universe[edge_index] = []
            else:
                iso = Isometry(-interval.width, 0)
                new_interval = iso(self.divided_universe[i])
                self.divided_universe[i] = new_interval

    def find_redundant_edges(self):
        """
        For each subinterval, checks if it is redundant i.e. if there exist pairings from this subinterval to another whose domains/ranges
        cover the entire subinterval (note that an empty subinterval is naturally redundant).
        Returns a list of edge indices corresponding to redundant edges.
        This check is not 100% accurate by design, it only checks to see if the subinterval can be split in two pieces at a certain endpoint.
        """
        redundant_edges = []
        for i in range(len(self.divided_universe)):
            subinterval = self.divided_universe[i]
            if not isinstance(subinterval, Interval):
                redundant_edges.append(i)
            else:
                intervals_on_edge = []
                for p in self.pairings:
                    if p.domain_index == i and p.range_index != i:
                        intervals_on_edge.append(p.domain)
                    elif p.domain_index != i and p.range_index == i:
                        intervals_on_edge.append(p.range)
                intervals_on_edge.sort(key=lambda I: I.start)

                endpoints = []
                for I in intervals_on_edge:
                    if I.start not in endpoints: endpoints.append(I.start)
                    if I.end not in endpoints: endpoints.append(I.end)

                for pt in endpoints:
                    if Interval(subinterval.start, pt) in intervals_on_edge and Interval(pt, subinterval.end) in intervals_on_edge:
                        redundant_edges.append(i)
                        break

                # Tried to take the union of intervals in order but didn't work
                # U = intervals_on_edge.pop(0)
                # while U:
                #     I = intervals_on_edge.pop(0)
                #     U = U.union(I)
                #
                # if U:
                #    if subinterval != U:
                #        significant_edges.append(i)
                # else:
                #     significant_edges.append(i)
        return redundant_edges

    def gcd_candidates_rev_type1(self, edge_index):
        """
        Looks for a pairing described below, this makes up the orientation reversing gcd picture combined with gcd_candidates_rev_type2.
        : finds an orientation reversing pairing that identifies the two halves of the subinterval at given index.
        Returns False if no such pairing exists.
        """
        subinterval = self.divided_universe[edge_index]
        pairings_on_edge = [p for p in self.pairings if p.domain_index == edge_index and p.range_index == edge_index and not p.is_preserving()]
        for p in pairings_on_edge:
            endpoints = [p.domain.start,  p.domain.end, p.range.start, p.range.end]
            if subinterval.start in endpoints and subinterval.end in endpoints:
                return p
        return False

    def gcd_candidates_rev_type2(self, edge_index):
        """
        Looks for a pairing described below, this makes up the orientation reversing  gcd picture combined with gcd_candidates_rev_type1.
        Finds pairs of orientation reversing pairings whose domains and ranges are all disjoint and union make up the entire subinterval.
        Returns False if no such pairs of pairings exist.
        """
        subinterval = self.divided_universe[edge_index]
        pairings_on_edge = []
        for p in self.pairings:
            if p.domain_index == edge_index and p.range_index == edge_index and not p.is_preserving():
                if p.domain.end == p.range.start or p.domain.start == p.range.end:
                    pairings_on_edge.append(p)

        candidates = []
        for i, p in enumerate(pairings_on_edge):
            endpoints_p = [p.domain.start, p.domain.end, p.range.start, p.range.end]
            p_start = min(endpoints_p)
            p_end = max(endpoints_p)
            for q in pairings_on_edge[i+1:]:
                endpoints_q = [q.domain.start, q.domain.end, q.range.start, q.range.end]
                q_start = min(endpoints_q)
                q_end = max(endpoints_q)
                if p_end == q_start and p_start == subinterval.start and q_end == subinterval.end:
                    candidates.append((p, q))
                elif q_end == p_start and q_start == subinterval.start and p_end == subinterval.end:
                    candidates.append((p, q))

        if len(candidates) == 0:
            return False
        else:
            return candidates

    def gcd_candidates_pre(self, edge_index):
        """
        Looks for 3 pairings described below that make up the orientation preserving gcd picture
        : finds pairs of orientation preserving pairings whose domains are disjoint and union make up the entire subinterval and whose range
        have the same property.
        Returns False if no such pairs of pairings exist.
        """
        subinterval = self.divided_universe[edge_index]
        pairings_on_edge = []
        for p in self.pairings:
            if p.domain_index == edge_index and p.range_index == edge_index and p.is_preserving():
                endpoints = [p.domain.start, p.domain.end, p.range.start, p.range.end]
                if subinterval.start in endpoints and subinterval.end in endpoints:
                    pairings_on_edge.append(p)

        candidates = []
        for i, p in enumerate(pairings_on_edge):
            for q in pairings_on_edge[i + 1:]:
                if p.domain.end == q.range.start and q.range.end == subinterval.end:
                    if p.range.start == q.domain.end and q.domain.start == subinterval.start:
                        candidates.append((p, q))

        if len(candidates) == 0:
            return False
        else:
            return candidates

    def simplify_transmit(self):
        """
        Simplifies the pseudogroup by transmitting pairings on edges where there exists a pairing sending the entire
        edge to another one.
        """
        # print statements to check pseudogroup
        # print('universe:', self.universe)
        # print('divided universe:', self.divided_universe)
        # print('pairings:', len(self.pairings))
        # for p in self.pairings:
        #     print(p)
        # print()
        # for i, I in enumerate(self.divided_universe):
        #     print(i, I.width, I, self.identify_edge_containing(I))

        transmitted_edges = []
        for i, I in enumerate(self.divided_universe):
            if self.identify_edge_containing(I):
                transmitted_edges.append(i)
                base_pairing = self.identify_edge_containing(I)

                pairings_on_edge = []
                for p in self.pairings:
                    if p.domain_index == i or p.range_index == i:
                        pairings_on_edge.append(p)

                self.pairings = [p for p in self.pairings if p not in pairings_on_edge]

                for p in pairings_on_edge:
                    new_p = base_pairing.transmit_inclusive(p)
                    if new_p.domain_index == i or new_p.range_index == i:
                        new_p = base_pairing.transmit_inclusive(new_p)
                        # weird case where domain/range of pairing are both included in domain/range of base_pairing and does the transmission
                        # the other way round, i.e. base_pairing: 2 -> 5, pairing: 5 -> 2, unluckily got 5 -> 5 so transmit again to get 2 -> 2
                    self.pairings.append(new_p)
                self.clean()

        for i in range(len(self.divided_universe)):
            if i in transmitted_edges:
                self.peel_edge(i)
        self.clean()
        self.pairings.sort()

        # print statements to check pseudogroup
        # print()
        # print('universe:', self.universe)
        # print('divided universe:', self.divided_universe)
        # print('peeled pairings:', len(self.pairings))
        # for p in self.pairings:
        #     print(p.domain.width, '/', p.isometry, '/', p)

    def find_candidates(self):
        """
        Looks for the gcd picture amongst the remaining pairings.
        If there are such pairings returns those pairings, otherwise returns False.
        """
        self.trim()
        self.clean()
        candidates_exist = False
        redundant_edges = self.find_redundant_edges()
        # print('redundant_edges:', redundant_edges)
        for i in range(len(self.divided_universe)):
            if i not in redundant_edges:
                if self.gcd_candidates_rev_type1(i) and self.gcd_candidates_rev_type2(i):
                    candidates_exist = True
                    # print('candidate ori rev 1', self.gcd_candidates_rev_type1(i))
                    # print('candidate ori rev 2')
                    # for pairs in self.gcd_candidates_rev_type2(i):
                        # print(pairs)
                    return [self.gcd_candidates_rev_type1(i), self.gcd_candidates_rev_type2(i)]
                if self.gcd_candidates_pre(i):
                    candidates_exist = True
                    # print('candidate ori pre')
                    # for pairs in self.gcd_candidates_pre(i):
                        # print(pairs)
                    return self.gcd_candidates_pre(i)
        if not candidates_exist:
            # print('No candidates')
            return False

        # endpoints_dup = []
        # endpoints = []
        # for p in self.pairings:
        #     endpoints_dup.append(p.domain.start)
        #     endpoints_dup.append(p.domain.end)
        #     endpoints_dup.append(p.range.start)
        #     endpoints_dup.append(p.range.end)
        # for pt in endpoints_dup:
        #     if pt not in endpoints:
        #         endpoints.append(pt)
        # endpoints = sorted(endpoints)
        #
        # print()
        # print(len(endpoints))
        # for pt in endpoints:
        #     print(pt)

        # # Check how many subintervals these points lie in to find an appropriate subinterval to transmit pairings into (not helpful)
        # count_dict = dict()
        # for pt in endpoints:
        #     count, unknown = self.num_pairings_including(pt, self.pairings)
        #     # print(pt, count, unknown)
        #     if count in count_dict:
        #         if pt not in count_dict[count]:
        #             count_dict[count].append(pt)
        #     else:
        #         count_dict[count] = [pt]
        # count_dict = dict(sorted(count_dict.items()))
        # # for key in count_dict:
        #     # print(key, count_dict[key])

    def reduce(self):
        """
        Does simplify_transmit() as much as possible until no subintervals can be transmitted.
        Then attempts to find the gcd picture with find_candidates().
        """
        transmission_possible = True
        while transmission_possible:
            self.simplify_transmit()
            transmission_possible = False
            for I in self.divided_universe:
                if not self.identify_edge_containing(I):
                    possible_transmissions = True
                    break
        return self.find_candidates()

    def reduce_amap(self):
        """
        Does simplify_transmit() as much as possible until no subintervals can be transmitted.
        Returns the interval(universe) and pairings as a result.
        """
        transmission_possible = True
        while transmission_possible:
            self.simplify_transmit()
            transmission_possible = False
            for I in self.divided_universe:
                if not self.identify_edge_containing(I):
                    possible_transmissions = True
                    break
        return self.universe, self.pairings


class Pseudogroup_comparable():
    """
    Pseudogroup where all polynomials are comparable.
    """
    def __init__(self, pairings, universe):
        self.pairings = pairings
        if isinstance(universe, Interval):
            self.universe = universe
        else:
            self.universe = ToInterval(universe)
        # self.check_comparable()

        # self.evaluations = []
        # num_var = universe.end.num_var
        # for comb in itertools.product(range(1, 40), repeat=num_var):
        #     assign = dict()
        #     for n in range(num_var):
        #         var_name = 'x' + str(n)
        #         assign[var_name] = comb[n]
        #     eval = evaluate_pseudogroup(self.universe, self.pairings, assign)
        #     self.evaluations.append(eval)

    def __repr__(self):
        result = 'Pseudogroup on %s:\n' % str(self.universe)
        if self.pairings:
            self.pairings.sort()
            for pairing in self.pairings:
                result += str(pairing) + '\n'
        return result

    def check_comparable(self):
        all_poly = []
        for p in self.pairings:
            all_poly.append(p.domain.start)
            all_poly.append(p.domain.end)
            all_poly.append(p.range.start)
            all_poly.append(p.range.start)

        for i, x in enumerate(all_poly):
            for y in all_poly[i+1:]:
                if (x < y) == 'unknown':
                    raise ValueError(f'Some polynomials are not comparable:{x}, {y}')

    # def evaluate(self):
    #     current_evaluations = []
    #     num_var = self.universe.end.num_var
    #     for comb in itertools.product(range(1, 40), repeat=num_var):
    #         assign = dict()
    #         for n in range(num_var):
    #             var_name = 'x' + str(n)
    #             assign[var_name] = comb[n]
    #         eval = evaluate_pseudogroup(self.universe, self.pairings, assign)
    #         current_evaluations.append(eval)
    #     if self.evaluations != current_evaluations:
    #         raise ValueError(f'The pseudogroup evaluations have changed:\n{self.evaluations}\n{current_evaluations}')

    def transform(self, transform):
        self.universe = Interval(linear_transform_poly(self.universe.start, transform), linear_transform_poly(self.universe.end, transform))
        self.pairings = [linear_transform_pairing(p, transform) for p in self.pairings]

    def clean(self):
        """
        Get rid of trivial Pairings and remove duplicates
        """
        self.pairings = [p.switch_domain_range() for p in self.pairings if not p.is_trivial()]
        remove_dup = []
        [remove_dup.append(p) for p in self.pairings if p not in remove_dup]
        self.pairings = remove_dup

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
        # for i in range(len(intervals)-1):
        #     print(intervals[i], '<', intervals[i+1], intervals[i] < intervals[i+1])
        I = intervals.pop(0)
        start, end = I.start, I.end
        if Polynomial(0).__lt__(start) == True:
            # print('1', Interval(Polynomial(0), start))
            return Interval(Polynomial(0), start)
        for interval in intervals:
            if end.__lt__(interval.start) == True:
                # print('2', Interval(end, interval.start))
                return Interval(end, interval.start)
            end = max(end, interval.end)
        if end.__lt__(self.universe.end) == True:
            # print('3', Interval(end, self.universe.end))
            return Interval(end, self.universe.end)
        # print('4', None)
        return None

    def contract(self):
        """
        Remove all static intervals. Return the total size.
        """
        result = Polynomial(0)
        I = self.static()
        while I:
            # print(I)
            # print(self.pairings)
            # print(self.universe)
            result += I.width
            if I.end != self.universe.end:
                self.pairings = [p.contract(I) for p in self.pairings]
            self.universe.set_end(self.universe.end - I.width)
            I = self.static()
        return result

    def merge(self):
        """
        Merge periodic pairing whenever possible.
        """
        if len(self.pairings) < 2:
            return
        done = 0
        while not done:
            periodics = [p for p in self.pairings if p.is_periodic()]
            done = 1
            for p in periodics[:-1]:
                for q in periodics[1 + periodics.index(p):]:
                    g = None
                    try:
                        g = p.merge(q)
                    except:
                        pass
                    if g:
                        self.pairings.remove(p)
                        self.pairings.remove(q)
                        self.pairings.append(g)
                        done = 0
                        break
                    # add this to break out of both loops when a periodic merger is performed
                    # (probably was intentioned since self.pairings continues to change but wasn't coded that way)
                    else:
                        continue
                break

    def transmit(self):
        """
        Use the largest Pairing to transmit others.
        """
        self.pairings.sort()
        g = self.pairings[0]
        self.pairings = [g] + [g.transmit(p) for p in self.pairings[1:]]

    def truncate(self):
        """
        Truncate the largest pairing.
        """
        self.pairings.sort()
        g = self.pairings.pop(0)
        if len(self.pairings) > 0:
            support_end = self.pairings[0].range.end
        else:
            support_end = g.range.start
        if (support_end < g.range.start) == True:
            self.universe.set_end(g.range.start)
            return
        if not g.is_preserving():
            g.trim()
        self.universe.set_end(support_end)
        range = Interval(g.range.start, support_end)
        domain = (~g.isometry)(range)
        self.pairings = [Pairing(domain, g.isometry, 0, 0)] + self.pairings

    def simplify(self):
        """
        Do one cycle of the orbit counting reduction algorithm due
        to Agol, Hass and Thurston.
        """
        self.check_comparable()
        self.clean()
        # self.evaluate()
        if len(self.pairings) == 0:
            self.pairings = None
            return self.universe.width
        # print("cleaned\n", self)
        count = self.contract()
        # print("contracted\n", self)
        self.clean()
        # self.evaluate()
        self.trim()
        # print("trimmed\n", self)
        self.clean()
        # self.evaluate()
        self.merge()
        # print("merged\n", self)
        self.clean()
        # self.evaluate()
        self.transmit()
        # print("transmitted\n", self)
        self.clean()
        # self.evaluate()
        self.truncate()
        # print("truncated\n", self)
        self.clean()
        # self.evaluate()
        # print('count = ', count)
        return count

    def reduce(self):
        """
        Reduce the pseudogroup to nothing. Return the number of orbits.
        """
        count = Polynomial(0)
        while self.pairings != None:
            try:
                count += self.simplify()
            except Exception as error:
                print(error)
                return self
        return count

def evaluate_interval(interval, assign):
    """
    Evaluates given interval with the variables in assign (must be a dictionary with the right variable names)
    and returns an orbits.Interval.
    """
    return orbits.Interval(interval.start.evaluate(**assign) + 1, interval.end.evaluate(**assign))

def evaluate_pairing(pairing, assign):
    """
    Evaluates given pairing with the variables in assign (must be a dictionary with the right variable names)
    and returns an orbits.Pairing.
    """
    domain_start = pairing.domain.start.evaluate(**assign)
    domain_end = pairing.domain.end.evaluate(**assign)
    iso_shift = pairing.isometry.shift.evaluate(**assign)
    iso_flip = pairing.isometry.flip
    # if pairing.domain.width == Polynomial(0) and pairing.domain.start == Polynomial(0):
    #     if pairing.is_trivial():
    #         return orbits.Pairing(orbits.Interval(1, 1), orbits.Isometry(0, iso_flip))
    #     if iso_flip:
    #         return orbits.Pairing(orbits.Interval(1, 1), orbits.Isometry(iso_shift + 1, iso_flip))
    #     else:
    #         return orbits.Pairing(orbits.Interval(1, 1), orbits.Isometry(iso_shift - 1, iso_flip))
    # if pairing.domain.width == Polynomial(0):
    #     return orbits.Pairing(orbits.Interval(domain_start, domain_end), orbits.Isometry(iso_shift, iso_flip))
    # else:
    if not iso_flip:
        return orbits.Pairing(orbits.Interval(domain_start + 1, domain_end), orbits.Isometry(iso_shift, iso_flip))
    if iso_flip:
        range_start = pairing.range.start.evaluate(**assign)
        range_end = pairing.range.end.evaluate(**assign)
        if domain_start < range_start:
            return orbits.Pairing(orbits.Interval(domain_start + 1, domain_end), orbits.Isometry(iso_shift + 1, iso_flip))
        else:
            return orbits.Pairing(orbits.Interval(range_start + 1, range_end), orbits.Isometry(iso_shift + 1, iso_flip))

def evaluate_pseudogroup(interval, pairings, assign):
    """
    Evaluates given interval and pairings of polynomials with the variables in assign (must be a dictionary with the right variable names)
    and returns the number of orbits.
    """
    assigned_interval = orbits.Interval(interval.start.evaluate(**assign) + 1, interval.end.evaluate(**assign))
    assigned_pairings = [evaluate_pairing(p, assign) for p in pairings]
    G = orbits.Pseudogroup(assigned_pairings, assigned_interval)
    return G.reduce()



# Below are some functions that are no longer used for Pseudogroup
# Saved in case of future use
class Pseudogroup_notused:
    def __init__(self, pairings, universe, divided_universe):
        self.pairings = pairings
        if isinstance(universe, Interval):
            self.universe = universe
        else:
            self.universe = ToInterval(universe)
        self.divided_universe = [ToInterval(I) for I in divided_universe]

    def static_left(self):
        """
        Find a leftmost static interval if there is one.
        """
        if len(self.pairings) == 0:
            return self.universe

        endpoints = []
        endpoints_dup = []
        for p in self.pairings:
            endpoints_dup.append(p.domain.start)
            endpoints_dup.append(p.domain.end)
            endpoints_dup.append(p.range.start)
            endpoints_dup.append(p.range.end)
        for pt in endpoints_dup:
            if pt not in endpoints:
                endpoints.append(pt)
        endpoints = sorted(endpoints)

        if Polynomial(0) in endpoints:
            print('No leftmost static interval, at least one starting at 0')
            return False
        else:
            print('0 not in any of the pairings')
            # Find the smallest point in the given endpoints if any
            for i in range(len(endpoints)):
                comp = set()
                for j in [x for x in range(len(endpoints)) if x != i]:
                    comp.add(endpoints[i].__lt__(endpoints[j]))
                if comp == {True}:
                    print(f'static interval: [0, {endpoints[i]}]')
                    return Interval(Polynomial(0), endpoints[i])
            print('No pairing starting at 0 but unknown whether there is a static interval')
            return False

    def contract_left(self):
        """
        Remove leftmost static interval if there is one. Return the total size.
        """
        I = self.static_left()
        if I:
            result = I.width
            self.universe.end = self.universe.end - I.end
            shifted_pairings = []
            for p in self.pairings:
                if p.isometry.flip == 0:
                    shifted_pairings.append(Pairing(Interval(p.domain.start - I.end, p.domain.end - I.end), p.isometry))
                else:
                    shifted_pairings.append(Pairing(Interval(p.domain.start - I.end, p.domain.end - I.end),
                                                    p.isometry * Isometry(I.end * Polynomial(2), 0)))
            self.pairings = shifted_pairings
            print(f'Contracted by {I.end}')
            return result
        else:
            return 0

    def transmit_overlap(self):
        """
        Transmit parings whose domain and range are the same until there are no intervals that are the same
        """
        # TO-DO:
        # -get all pairs of pairings (put this into a list of pairs?)
        # -check ones that have overlapping domain/range and perform appropriate modification

        pairings = [p for p in self.pairings]
        for i, pairing in enumerate(pairings):
            domain = pairing.domain
            range = pairing.range
            for other_pairing in (pairings[:i] + pairings[i+1:]):
                if domain == other_pairing.domain:
                    pairings.remove(pairing)
                    pairings.remove(other_pairing)
                    count_start, unknown_start = self.num_pairings_including(domain.start, pairings)
                    count_end, unknown_end = self.num_pairings_including(domain.end, pairings)
                    if [count_start, count_end, unknown_start, unknown_end] == [0, 0, 0, 0]:
                        new_iso = other_pairing.isometry * ~pairing.isometry
                        assert Interval(new_iso(range.start), new_iso(range.end)) == other_pairing.range
                        print('one overlapping pairing transmitted')
                        self.pairings = pairings
                        return
                elif domain == other_pairing.range:
                    new_iso = ~other_pairing.isometry * ~pairing.isometry
                    assert Interval(new_iso(range.start), new_iso(range.end)) == other_pairing.domain
                    pairings.remove(pairing)
                    pairings.append(Pairing(range, new_iso))
                    print('one overlapping pairing transmitted')
                    self.pairings = pairings
                    return
                elif range == other_pairing.domain:
                        new_iso = other_pairing.isometry * pairing.isometry
                        assert Interval(new_iso(domain.start), new_iso(domain.end)) == other_pairing.range
                        pairings.remove(pairing)
                        pairings.append(Pairing(domain, new_iso))
                        print('one overlapping pairing transmitted')
                        self.pairings = pairings
                        return
                elif range == other_pairing.range:
                        new_iso = ~other_pairing.isometry * pairing.isometry
                        assert Interval(new_iso(domain.start), new_iso(domain.end)) == other_pairing.domain
                        pairings.remove(pairing)
                        pairings.append(Pairing(domain, new_iso))
                        print('one overlapping pairing transmitted')
                        self.pairings = pairings
                        return
        print('no overlaps')
        return

    def transmit_inclusive(self):
        """
        Transmit parings whose domain and range lie entirely in the domain or range of another pairing
        Is only performed to transmit to the right
        """
        pairings = [p for p in self.pairings]
        for i, pairing in enumerate(pairings):
            domain = pairing.domain
            range = pairing.range
            for other_pairing in (pairings[:i] + pairings[i+1:]):
                if other_pairing.domain.__contains__(domain) == True:
                    new_iso = other_pairing.isometry * ~pairing.isometry
                    new_pairing = Pairing(range, new_iso)
                    if pairing.domain.start.__lt__(new_pairing.range.start) == True:
                        pairings.remove(pairing)
                        pairings.append(new_pairing)
                        print('one pairing transmitted')
                        self.pairings = pairings
                        return
                elif other_pairing.range.__contains__(domain) == True:
                    new_iso = ~other_pairing.isometry * ~pairing.isometry
                    new_pairing = Pairing(range, new_iso)
                    if pairing.domain.start.__lt__(new_pairing.range.start) == True:
                        pairings.remove(pairing)
                        pairings.append(new_pairing)
                        print('one pairing transmitted')
                        self.pairings = pairings
                        return
                elif other_pairing.domain.__contains__(range) == True:
                    new_iso = other_pairing.isometry * pairing.isometry
                    new_pairing = Pairing(domain, new_iso)
                    if pairing.range.start.__lt__(new_pairing.range.start) == True:
                        pairings.remove(pairing)
                        pairings.append(new_pairing)
                        print('one pairing transmitted')
                        self.pairings = pairings
                        return
                elif other_pairing.range.__contains__(range) == True:
                    new_iso = ~other_pairing.isometry * pairing.isometry
                    new_pairing = Pairing(domain, new_iso)
                    if pairing.range.start.__lt__(new_pairing.range.start) == True:
                        pairings.remove(pairing)
                        pairings.append(new_pairing)
                        print('one pairing transmitted')
                        self.pairings = pairings
                        return
        print('no subintervals included')
        return

    def truncate_left(self):
        """
        Truncate the leftmost interval if it is only contained in one pairing
        """
        pairings = [p for p in self.pairings]
        first_int = []
        first_pairing = []
        for p in pairings:
            if p.domain.start == Polynomial(0):
                first_int.append(p.domain.end)
                first_pairing.append(p)
            elif p.range.start == Polynomial(0):
                first_int.append(p.range.end)
                first_pairing.append(p)

        for i, endpoint in enumerate(first_int):
            count, unknown = self.num_pairings_including(endpoint, pairings)
            if count == 1 and unknown == 0:
                self.universe.end = self.universe.end - endpoint
                pairings.remove(first_pairing[i])
                shifted_pairings = []
                for p in pairings:
                    if p.isometry.flip == 0:
                        shifted_pairings.append(Pairing(Interval(p.domain.start - endpoint, p.domain.end - endpoint), p.isometry))
                    else:
                        shifted_pairings.append(Pairing(Interval(p.domain.start - endpoint, p.domain.end - endpoint),
                                                        p.isometry * Isometry(endpoint * Polynomial(2), 0)))
                print('truncated leftmost subinterval')
                self.pairings = shifted_pairings
                return
        print('cannot truncate')
        return

    def num_pairings_including(self, x, list_pairings):
        count = 0
        unknown = 0
        for p in list_pairings:
            inc_domain = p.domain.__contains__(x)
            inc_range = p.range.__contains__(x)
            if inc_domain == True or inc_range == True:
                count += 1
            elif inc_domain == 'unknown' or inc_range == 'unknown':
                unknown += 1
        return count, unknown

    def simplify(self):
        print('universe')
        print(self.universe)
        print('pairings')
        print(len(self.pairings))
        for p in self.pairings:
            print(p, p.isometry)

        print()

        count = self.contract_left()

        self.clean()
        print('removed trivial pairings')
        for p in self.pairings:
            print(p, p.isometry)
        self.trim()
        print('trimmed necessary pairings')
        for p in self.pairings:
            print(p, p.isometry)

        print()
        # two variables used for polynomials in given example
        # u = list(self.pairings[1].range.start.poly.keys())[0]
        # v = list(self.pairings[1].range.start.poly.keys())[1]
        # x = Polynomial(**{u: 10, v: 14})

        # get rid of overlapping pairings with overlapping domain/range
        flag = True
        current_step = [p for p in self.pairings]
        while flag:
            self.transmit_overlap()
            next_step = [p for p in self.pairings]
            flag = (next_step != current_step)
            current_step = next_step
        print(len(self.pairings))
        for p in self.pairings:
            print(p, p.isometry)

        print()

        # transmit pairings whose domain/range are contained in others
        flag = True
        current_step = [p for p in self.pairings]
        while flag:
            self.transmit_inclusive()
            next_step = [p for p in self.pairings]
            flag = (next_step != current_step)
            current_step = next_step
        print(len(self.pairings))
        for p in self.pairings:
            print(p, p.isometry)

        print()

        # # Take all start, endpoints of domains and ranges
        # endpoints_dup = []
        # endpoints = []
        # for p in self.pairings:
        #     endpoints_dup.append(p.domain.start)
        #     endpoints_dup.append(p.domain.end)
        #     endpoints_dup.append(p.range.start)
        #     endpoints_dup.append(p.range.end)
        # for pt in endpoints_dup:
        #     if pt not in endpoints:
        #         endpoints.append(pt)
        # endpoints = sorted(endpoints)
        # # for pt in endpoints:
        # #     print(pt)
        #
        # # Check how many subintervals these points lie in to find an appropriate subinterval to transmit pairings into (not helpful)
        # count_dict = dict()
        # for pt in endpoints:
        #     count, unknown = self.num_pairings_including(pt, self.pairings)
        #     # print(pt, count, unknown)
        #     if count in count_dict:
        #         if pt not in count_dict[count]:
        #             count_dict[count].append(pt)
        #     else:
        #         count_dict[count] = [pt]
        # count_dict = dict(sorted(count_dict.items()))
        # # for key in count_dict:
        #     # print(key, count_dict[key])

        # truncate leftmost pairing if possible
        self.truncate_left()
        print(len(self.pairings))
        for p in self.pairings:
            print(p, p.isometry)


        # # Take all domains and ranges of pairings
        # subint = []
        # for p in self.pairings:
        #     subint.append(p.domain)
        #     subint.append(p.range)
        #
        # for i, interval in enumerate(subint):
        #     for nextinterval in subint[i+1:]:
        #         if nextinterval.__contains__(interval) == True:
        #             print(interval, nextinterval)
        #             print(type(nextinterval))

        return count