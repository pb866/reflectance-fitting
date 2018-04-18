# -*- coding: utf-8 -*-
"""
Created on Fri Apr 13 17:46:36 2018

@author: rturley
"""

import unittest
import refl
import numpy as np

class refl_test(unittest.TestCase):
    def setUp(self):
        self.AlIndex = refl.Index('Al')
        self.SiIndex = refl.Index('Si')
        self.SiO2Index = refl.Index('SiO2')

    def test_Al(self):
        self.assertAlmostEqual(self.AlIndex.at(30.4),0.951+0.006j,3)
        
    def test_Si(self):
        self.assertAlmostEqual(self.SiIndex.at(30.4),0.9293+0.009j,3)
    
    def test_Al_refl(self):
        lam = 30.4
        n=np.array([1,self.AlIndex.at(lam), self.SiIndex.at(lam)])
        t=np.array([0,50,0])
        thetad=15
        sigma=0
        r=refl.matR(n, t, thetad, lam, sigma)
        self.assertAlmostEqual(r[0],0.724,3)
        
    def test_s_pol(self):
        self.assertAlmostEqual(refl.fracs(30.4),0.93,2)
        
    def test_matR_refl_comp(self):
        lam = 15
        n = np.array([1, self.AlIndex.at(lam), self.SiO2Index.at(lam)])
        t = np.array([0, 20, 0])
        thetad = 20
        rp = refl.Parratt(n, t, thetad, lam)
        rm = refl.matR(n, t, thetad, lam)[0]
        self.assertAlmostEqual(rm, 0.01260649028, 7)
        self.assertAlmostEqual(rp, rm, 8)

if __name__ == '__main__':
    unittest.main()