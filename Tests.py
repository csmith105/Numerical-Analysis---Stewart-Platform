from StewartPlatform import StewartPlatform

from math import *

def test1():

    print("Test 1 - Verify Activity 1")

    stewie = StewartPlatform(2.0, sqrt(2.0), sqrt(2.0), pi / 2.0, 4.0, 0.0, 0.0, 4.0)
    stewie.setP(sqrt(5.0), sqrt(5.0), sqrt(5.0))

    f1 = stewie.f(-pi / 4.0)
    f2 = stewie.f(pi / 4.0)

    if(f1 < 0.000001):
        print(str(f1) + " is acceptably close to zero")

    if(f2 < 0.000001):
        print(str(f2) + " is acceptably close to zero")

    if (f1 < 0.000001 and f2 < 0.000001):
        print("PASS")
    else:
        print("FAIL")
        assert()

def test2():

    print("Test 2 - Verify Activity 2")

    stewie = StewartPlatform(2.0, sqrt(2.0), sqrt(2.0), pi / 2.0, 4.0, 0.0, 0.0, 4.0)
    stewie.setP(sqrt(5.0), sqrt(5.0), sqrt(5.0))

    print("Plotting...")
    stewie.plotF(-pi / 2.0, pi / 2.0, "Activity 2 - Plotting f from -PI/2 to PI/2")

    print("PASS")

def test3():

    print("Test 3 - Verify Activity 3")

    stewie = StewartPlatform(2.0, sqrt(2.0), sqrt(2.0), pi / 2.0, 4.0, 0.0, 0.0, 4.0)
    stewie.setP(sqrt(5.0), sqrt(5.0), sqrt(5.0))

    print("Plotting...")
    stewie.setTheta(-pi / 4)
    stewie.plotTriangle("Activity 3 - Theta = -PI/4")

    print("Plotting...")
    stewie.setTheta(pi / 4)
    stewie.plotTriangle("Activity 3 - Theta = PI/4")

    print("PASS")

def test4():

    print("Test 4 - Verify Activity 4")

    stewie = StewartPlatform(3.0, 3 * sqrt(2.0), 3.0, pi / 4.0, 5, 0.0, 0.0, 6.0)
    stewie.setP(5.0, 5.0, 3.0)

    print("Plotting...")
    stewie.plotF(-pi, pi, "Activity 4 - Plotting f from -PI/2 to PI/2")

    print("Solving...")

    results = stewie.solve()

    print(results)

    for result in results:

        print("Plotting...")

        stewie.setTheta(result)
        stewie.plotTriangle("Activity 4 - Theta = " + str(result))

    print("PASS")

test1()
test2()
test3()
test4()




