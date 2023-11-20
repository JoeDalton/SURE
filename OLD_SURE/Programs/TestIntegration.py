#!/usr/bin/env python2.7
import Link
import math
import matplotlib.pyplot as plt
from Numerics.Integration import Integrate

def func0(x):
    return 1


def func1(x):
    return x


def func2(x):
    return x*x


def func3(x):
    return x*x*x


def func4(x):
    return math.exp(x)


def func5(x):
    return 1/math.sqrt(x+1)


def func6(x):
    return math.cos(x)


def func7(xList):
    return 1


def func8(xList):
    result = 1.0
    for i in range(0, len(xList)):
        result *= xList[i]
    result +=1
    return result


def func9(xList):
    result = 1.0
    for i in range(0, len(xList)):
        result *= func6(xList[i])
    result +=1
    return result

def func10(xList):
  return func4(xList[1])




print('#########################')
print('#  Test ND integration  #')
print('#########################')
func   = func10
func1D = func6
dim    = 3
OneDBounds = [-2,2]
bounds = [OneDBounds]*dim
nPoint = 10000
level  = 4


print("")
print('Reference')
print(Integrate(func1D, OneDBounds, method='SecondFejer', level=level+3))
print('Simpson')
print(Integrate(func1D, OneDBounds, method='Simpson', level=2))

#print(' ----------------------- ')
#print('|      New methods      |')
#print(' ----------------------- ')
#
print("")
print("Full Second Fejer")
print("Dimension = " + str(dim))
print('level = ' + str(level))
print(Integrate(func, bounds, samplingMethod='Cubature', quadratureMethod='SecondFejer', maxLevel=level, isSparse=False)) #, plot=True))

#print("")
#print("Full CC")
#print("Dimension = " + str(dim))
#print('level = ' + str(level))
#print(Integrate(func, bounds, samplingMethod='Cubature', quadratureMethod='Clenshaw-Curtis', maxLevel=level, isSparse=False)) #, plot=True))

print("")
print("Sparse Second Fejer")
print("Dimension = " + str(dim))
print('level = ' + str(level))
print(Integrate(func, bounds, samplingMethod='Cubature', quadratureMethod='SecondFejer', maxLevel=level, isSparse=True))

#print("")
#print("Sparse CC")
#print("Dimension = " + str(dim))
#print('level = ' + str(level))
#print(Integrate(func, bounds, samplingMethod='Cubature', quadratureMethod='Clenshaw-Curtis', maxLevel=level, isSparse=True))

