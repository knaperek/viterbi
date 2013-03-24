#!/usr/bin/env python3
"""
Tests for convolutional encoder and Viterbi decoder routines.

"""
__author__ = "Jozef Knaperek"


import unittest
from random import randint, shuffle
from viterbi import *

class TestFunctions(unittest.TestCase):

    def setUp(self):
        self.n_state_bits = randint(2, 10)

        n_possible_polynomials = 2**(self.n_state_bits+1)-(self.n_state_bits+2)
        max_polynomials = min( # do not use more than 10 polynomials
            randint(3, 10),
            randint(1, n_possible_polynomials)
        )

        pols = set()  # set of polynomials
        for i in range(max_polynomials):
            ones = randint(2, self.n_state_bits + 1) # number of 1s in the polynomial
            p_bitlist = [1]*ones + [0]*(self.n_state_bits + 1 - ones)
            shuffle(p_bitlist)
            p = 0
            while p_bitlist:
                p = (p << 1) | p_bitlist.pop()
            pols.add(p)
        polynomials = list(pols)
        self.n_polynomials = len(polynomials)

        self.transitions = Transitions(polynomials)

        self.input = BinData()
        for bit in range(randint(10, 10000)):
            self.input += bit

        # print('Testing with {} state bits, {} polynomials and {} input data bits'.format(self.n_state_bits, self.n_polynomials, self.input.len))

    def test_a_formal(self):
        """ Tests formal side of algorithm. Any catched error is fatal and means bug in implementation! """
        data = self.input
        encoded = self.transitions.encode(data)
        encoded_2 = self.transitions.encode(data)
        decoded = self.transitions.decode(encoded)
        decoded_2 = self.transitions.decode(encoded)

        self.assertEqual(encoded.len, data.len * self.n_polynomials)  # Encoded data length check
        self.assertEqual(decoded.len, data.len)  # Decoded data length check

        self.assertEqual(encoded, encoded_2)  # Encoding should give always the same output
        self.assertEqual(decoded, decoded_2)  # Decoding should give always the same output


    def test_b_errorfree(self):
        """ Tests decoding after error-free transmission. """
        data = self.input
        encoded = self.transitions.encode(data)
        decoded = self.transitions.decode(encoded)

        hd = hamming_distance(data.num, decoded.num)
        self.assertEqual(data, decoded, 'hamming distance: {}'.format(hd))  # Following must be true: Decode(Encode(Data)) == Data

        self.assertEqual(hd, 0, 'hamming distance of "equal" bit-streams is not zero!')  # if data == decoded, then hamming distance must be 0!

    def test_c_single_error(self):
        """ Tests decoding after a single bit error introduced during "transmission". """
        data = self.input
        encoded = self.transitions.encode(data)

        bit_to_flip = randint(0, encoded.len)
        error_vector = 1 << bit_to_flip

        encoded_err = encoded ^ error_vector
        decoded = self.transitions.decode(encoded_err)

        hd = hamming_distance(data.num, decoded.num)
        self.assertTrue(data == decoded, 'hamming distance: {}'.format(hd))  # Following should be true: Decode(Encode(Data)) == Data

    def test_d_more_errors(self):
        """ Tests decoding after a few bit errors introduced during "transmission". """
        data = self.input
        encoded = self.transitions.encode(data)

        error_vector = 0
        for i in range(int(randint(1, 20) / 100. * encoded.len)):
            bit_to_flip = randint(0, encoded.len)
            error_vector |= 1 << bit_to_flip

        encoded_err = encoded ^ error_vector
        decoded = self.transitions.decode(encoded_err)

        n_errors = hamming_distance(encoded.num, encoded_err.num)
        print('{} of {} bit errors ({}%)'.format(n_errors, encoded.len, int(n_errors / encoded.len * 100)))

        hd = hamming_distance(data.num, decoded.num)
        self.assertTrue(data == decoded, 'hamming distance: {}'.format(hd))  # Following should/could be true: Decode(Encode(Data)) == Data


    def tearDown(self):
        print('Tested with {} state bits, {} polynomials and {} input data bits'.format(self.n_state_bits, self.n_polynomials, self.input.len))

if __name__ == '__main__':
    unittest.main()
