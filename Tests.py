from StewartPlatform import StewartPlatform

from math import *

def test1():

    stewie = StewartPlatform(2.0, sqrt(2), sqrt(2), pi / 2.0)
    stewie.setP(sqrt(5), sqrt(5), sqrt(5))

    print(stewie.f(-pi / 4.0))
    print(stewie.f(pi / 4.0))

    # assert(stewie.f(-pi / 4.0) == 0.0)
    # assert(stewie.f(pi / 4.0) == 0.0)

    return


print("Test 1 - Verify Activity 1")
test1()




