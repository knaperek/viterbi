#!/usr/bin/env python3
"""
Viterbi algorithm for decoding convolutional codes.

"""
__author__ = "Jozef Knaperek"

import sys
import argparse
from collections import namedtuple
from bindata import BinData


def hamming_distance(a, b):
    """ Returns Hamming distance of two codes """
    ones = 0  # number of 1s in xorred value (equals to HD)
    xorred = a ^ b
    while(xorred):
        ones += xorred & 1
        xorred >>= 1
    return ones

def parity(number):
    """ xors all bits in number and returns the result """
    par = 0  # inicialize to even parity
    while number:
        par ^= number & 1
        number >>= 1
    return par

# Least significant bit in polynomial matches new bit (comming into window)
# Parity computed with first polynomial in the list will be stored in most significant bit of the resulting parity bitstream
def conv_parity(window, polynomials):
    """ Calculates and returns window parity bitstream according to given list of polynomials """
    par_bits = (parity(window & pol) for pol in polynomials)  # select only bits that are significant for given polynomial
    code = 0
    for bit in par_bits:
        code = (code << 1) | bit
    return code

def bits_count(number):
    """ Returns minimum number of bit digits required to store the number. """
    i = 0
    while (1 << i) <= number:
        i += 1
    return i or 1

def ones_count(number):
    """ Returns the number of binary 1s in the number. """
    number = int(number)  # make a copy
    ones = 0
    while number:
        ones += number & 1
        number >>= 1
    return ones


# stav[3][0/1] = (new_state, parity)  # named tuples used - Edge(new_state, parity)
Edge = namedtuple('Edge', 'new_state, parity')

class Transitions(object):
    """docstring for Transitions"""
    def __init__(self, polynomials):
        self.n_state_bits = bits_count(max(polynomials)) - 1
        self.n_states = 2**self.n_state_bits
        self.polynomials = polynomials
        self.parity_len = len(polynomials)  # length of the parity per one data bit
        self.states =   [
                            [
                                Edge(new_state=((i_state<<1)&(2**self.n_state_bits-1))|bit,
                                     parity=conv_parity(i_state<<1|bit, polynomials) )
                            for bit in range(2)
                            ]
                        for i_state in range(self.n_states)
                        ]

    def __str__(self):
        output = ''
        for i, state in enumerate(self.states):
            state_code = '{{:0{}b}}'.format(self.n_state_bits).format(i)
            lines = [' --{{}}/{{:0{}b}}--> {{:0{}b}}\n'.format(self.parity_len, self.n_state_bits).format(b, state[b].parity, state[b].new_state) for b in (0, 1)]
            output += state_code + lines[0] + ' '*len(state_code) + lines[1] + '\n'
        return output

    def generate_parities(self, bindata):
        """ Returns generator that gives (encoding) parity checks for input data """
        # It starts with most significant bits in data, going to least significant
        state = 0
        for i in reversed(range(bindata.len)):
            bit = bindata[i]
            yield self.states[state][bit].parity
            state = self.states[state][bit].new_state

    def encode(self, bindata):
        """ Encodes data using convolutional code. Public method (API). """
        parity_sequence = BinData(0, 0)
        for parity in self.generate_parities(bindata):
            parity_sequence += BinData(parity, self.parity_len)
        return parity_sequence

    def extract_parity_sequence(self, parity_sequence_bindata):
        """ Returns generator iterating through parities in parity sequence (encoded data). """
        parity_mask = (1 << self.parity_len) - 1
        parity_selector = parity_sequence_bindata.len  # number of least signifficant bits to be discarded (>>) for parity to be readable by parity_mask
        while parity_selector:
            parity_selector -= self.parity_len
            yield (parity_sequence_bindata.num & (parity_mask << parity_selector)) >> parity_selector

    def decode(self, parity_sequence_bindata):
        """ Decodes convolutional code using the Viterbi algorithm. Public method (API). """
        gen = self.extract_parity_sequence(parity_sequence_bindata)
        state = 0  # initial state

        INF = float('inf')  # constant definition
        class Node():
            def __init__(self, metric=INF, bindata=None):
                self.metric = metric
                self.bindata = bindata or BinData(0, 0)

        # init trellis
        olds = [ Node(INF, BinData(0, 0)) for i in range(self.n_states)]  # aktualna metrika, data bits
        news = [ Node(None, None) for i in range(self.n_states)]  # nova metrika, data bits (with added one new bit)
        olds[0].metric = 0  # set metrics of first state to 0

        # iterate through parities in encoded data (parity sequence)
        for parity in gen:
            # initialize news 
            for new in news:
                new.metric = INF  # set new PM to infinity

            # choose best paths for new step
            for i in range(self.n_states):
                for bit in (0, 1): 
                    t = self.states[i][bit].new_state
                    p = self.states[i][bit].parity
                    hd = hamming_distance(p, parity)

                    new_PM = olds[i].metric + hd  # compute candidate PM
                    if new_PM < news[t].metric:  # if this new candidate is better than existing new candidate
                        news[t].metric = new_PM
                        news[t].bindata = olds[i].bindata + bit

            # update "column" in trellis with best paths chosen in previous step and prepare for next iteration
            for i in range(self.n_states):
                olds[i].metric = news[i].metric
                olds[i].bindata = news[i].bindata

        # Finalization
        # Get state with best PM
        best_state, best_PM = None, INF
        for old in olds:
            if old.metric < best_PM:
                best_PM = old.metric
                best_state = old

        # Decoded databits:
        return best_state.bindata


def main():
    parser = argparse.ArgumentParser(description="-= Encoder/Decoder of convolutional codes.\nAuthor: Jozef Knaperek =-\n")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-e', '--encode', action='store_true', help='encode data with convolutional code')
    group.add_argument('-d', '--decode', action='store_true', help='decode data using Viterbi algorithm')

    def make_bindata(arg):
        return BinData(arg)
    make_bindata.__name__ = 'input data'
    parser.add_argument('-i', '--input', type=make_bindata, help='input data bit-stream (instead of using stdin)')

    def make_pols_list(arg):
        pols = [int(pol, 2) for pol in arg.split(',')]
        if min(map(ones_count, pols)) < 2:
            raise ValueError('Every valid polynomial must have at least two binary 1s')
        return pols
    make_pols_list.__name__ = 'polynomials list'
    parser.add_argument('polynomials', help='comma separated list of binnary polynomials (of at least two binary 1s in each)', type=make_pols_list)

    args = parser.parse_args()

    try:
        input_data = args.input or BinData(sys.stdin.read())
    except ValueError:
        sys.exit('Invalid input data: ' + stdin_input)

    if args.encode:  # encode
        print(Transitions(args.polynomials).encode(input_data))
    else:  # decode
        if len(input_data) % len(args.polynomials):
            sys.exit('Decoding error: The number of data bits ({}) is not multiple of the number of polynomials ({})!'.format(len(input_data), len(args.polynomials)))
        print(Transitions(args.polynomials).decode(input_data))


if __name__ == '__main__':
    main()
