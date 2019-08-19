import pickle
from Ising import *

import numpy as np 
from matplotlib import pyplot as plt



pickleTest = Ising(10,10)

pickle.dump(pickleTest, open("pickleTest.p","wb"))

pickleTest_Result = pickle.load(open("pickleTest.p","rb"))


