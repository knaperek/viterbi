"""
Simple implementation of class handling arbitrary long binary data sequences.
"""

class BinData:
    """ Wrapper of standard int/long type. Keeps track of the number of bits occupied. """

    def __init__(self, number=0, length=0):
        if isinstance(number, int):  # first try to deal with number as it is regular int
            self.num = number
            self.len = length
        else:  # ...or try dealing with it as it is a string (in binary representation)
            number = number.strip()
            self.num = int(number, 2)
            self.len = len(number)

    def __str__(self):
        return '{{:0{}b}}'.format(self.len).format(self.num)

    def __repr__(self):
        return '{}({}, {})'.format(self.__class__.__name__, self.num, self.len)

    def __eq__(self, other):
        return (self.num == other.num) and (self.len == other.len)

    def __ne__(self, other):
        return (self.num != other.num) or (self.len != other.len)

    def __lt__(self, other):
        return self.num < other.num

    def __le__(self, other):
        return self.num <= other.num

    def __gt__(self, other):
        return self.num > other.num

    def __ge__(self, other):
        return self.num >= other.num

    def __bool__(self):  # True if len != 0
        return bool(self.len)

    def __len__(self):
        return self.len

    def __getitem__(self, index):
        if index < 0:  # support negative indexes
            index += self.len
        if index >= self.len or index < 0:
            raise IndexError
        return 1 if (1 << index) & self.num else 0

    def __iter__(self):
        return ( self[i] for i in reversed(range(self.len)) )

    def __reversed__(self):
        return ( self[i] for i in range(self.len) )

    def __add__(self, other):
        if isinstance(other, self.__class__):
            return self.__class__((self.num << other.len) | other.num, self.len + other.len)
        else:
            bit = other & 1  # if other is not instance of this class, it's considered as single bit
            return self.__class__((self.num << 1) | bit, self.len + 1)

    def __xor__(self, other):
        if isinstance(other, self.__class__):
            return self.__class__(self.num ^ other.num, max(self.len, other.len))
        else:
            return self.__class__(self.num ^ other, self.len)
