# -*- coding: utf-8 -*-
"""
Created on Fri Apr 13 17:46:36 2018

@author: rturley
"""

import unittest
import numpy as np

def fun(x):
    return x+1

class MyTest(unittest.TestCase):
    def test(self):
        self.assertEqual(fun(3),4)
        
    def test2(self):
        self.assertLess(fun(4),6)
        self.assertGreater(fun(4),2)
        
    def test3(self):
        self.assertTrue(fun(3)>3)
        self.assertAlmostEqual(np.pi,3.14159,5)
        
if __name__ == '__main__':
    unittest.main()