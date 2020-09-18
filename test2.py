# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 23:58:18 2020

@author: anavr
"""

import matplotlib.pyplot as plt
import Comtrade

rec = Comtrade()
rec.load("sample_files/sample_ascii.cff")
print("Trigger time = {}s".format(rec.trigger_time))

plt.figure()
plt.plot(rec.time, rec.analog[0])
plt.plot(rec.time, rec.analog[1])
plt.legend([rec.analog_channel_ids[0], rec.analog_channel_ids[1]])
plt.show()