"""Verification that we can run tests using the pytest command."""

import math
import os
import sys

sys.path.append(os.path.abspath(os.path.join('..', 'mongo-rdkit')))


def test_sqrt():
    num = 25
    assert math.sqrt(num) == 5

def test_squareFail():
    num = 3
    assert num * num != 2

def test_equals():
    num = 5
    assert num == num
