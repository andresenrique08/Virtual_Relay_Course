# -*- coding: utf-8 -*-
"""
Created on Sat Sep 26 13:37:02 2020

@author: uer
"""

import pyComtrade
comtradeObj = pyComtrade.ComtradeRecord()
#comtradeObj.read('./Comtrade TEST 1-1/Nasser/Sampling1M.cfg', './Comtrade TEST 1-1/Nasser/Sampling1M.dat')
#comtradeObj.read('../Comtrade TEST 1-1/L1L2L3_75%/BINAIRE/Test_1-1_ABC60_75%_BIN.cfg', '../Comtrade TEST 1-1/L1L2L3_75%/BINAIRE/Test_1-1_ABC60_75%_BIN.dat')
comtradeObj.read('C:/Users/uer/Documents/GitHub/Virtual_Relay_Course/Comtrade TEST 1-1/L1L2L3_75%/ASCII/Test_1-1_ABC60_75%_ASCII.cfg', 'C:/Users/uer/Documents/GitHub/Virtual_Relay_Course/Comtrade TEST 1-1/L1L2L3_75%/ASCII/Test_1-1_ABC60_75%_ASCII.dat')
ns = comtradeObj['endsamp'][-1]
v=len(comtradeObj['A'][2]['values'])
print(v)
time = comtradeObj.get_timestamps()
print(len(time))