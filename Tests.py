from StewartPlatform import StewartPlatform2D, StewartPlatform3D

from math import *
import matplotlib.pyplot as plt
import numpy as np

def test1():

    print("Test 1 - Verify Activity 1")

    stewie = StewartPlatform2D(2.0, sqrt(2.0), sqrt(2.0), pi / 2.0, 4.0, 0.0, 0.0, 4.0)
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

    stewie = StewartPlatform2D(2.0, sqrt(2.0), sqrt(2.0), pi / 2.0, 4.0, 0.0, 0.0, 4.0)
    stewie.setP(sqrt(5.0), sqrt(5.0), sqrt(5.0))

    print("Plotting...")
    stewie.plotF(-pi / 2.0, pi / 2.0, "Activity 2 - Plotting f from -PI/2 to PI/2")

    print("PASS")

def test3():

    print("Test 3 - Verify Activity 3")

    stewie = StewartPlatform2D(2.0, sqrt(2.0), sqrt(2.0), pi / 2.0, 4.0, 0.0, 0.0, 4.0)
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

    stewie = StewartPlatform2D(3.0, 3 * sqrt(2.0), 3.0, pi / 4.0, 5, 0.0, 0.0, 6.0)
    stewie.setP(5.0, 5.0, 3.0)

    print("Plotting...")
    stewie.plotF(-pi, pi, "Activity 4 - Plotting f from -PI to PI")

    print("Solving...")

    results = stewie.solve()

    print(results)

    for result in results:

        print("Plotting...")

        stewie.setTheta(result)
        stewie.plotTriangle("Activity 4 - Theta = " + str(result))

    print("PASS")

def test5():

    print("Test 5 - Verify Activity 5")

    stewie = StewartPlatform2D(3.0, 3 * sqrt(2.0), 3.0, pi / 4.0, 5, 0.0, 0.0, 6.0)
    stewie.setP(5.0, 7.0, 3.0)

    print("Plotting...")
    stewie.plotF(-pi, pi, "Activity 5 - Plotting f from -PI to PI")

    print("Solving...")

    results = stewie.solve()

    print(results)

    for result in results:

        print("Plotting...")

        stewie.setTheta(result)
        stewie.plotTriangle("Activity 5 - Theta = " + str(result))

    print("PASS")

def test6():

    print("Test 6 - Verify Activity 6")

    stewie = StewartPlatform2D(3.0, 3 * sqrt(2.0), 3.0, pi / 4.0, 5, 0.0, 0.0, 6.0)
    stewie.setP(5.0, 9.0, 3.0)

    print("Plotting...")
    stewie.plotF(-pi, pi, "Activity 6 - Plotting f from -PI to PI")

    print("Solving...")

    results = stewie.solve()

    print(results)

    for result in results:

        print("Plotting...")

        stewie.setTheta(result)
        stewie.plotTriangle("Activity 6 - Theta = " + str(result))

    print("PASS")

def test7():

    print("Test 7 - Verify Activity 7")

    stewie = StewartPlatform2D(3.0, 3 * sqrt(2.0), 3.0, pi / 4.0, 5, 0.0, 0.0, 6.0)

    xValues = []
    yValues = []

    for i in np.arange(0, 11, 0.1):

        stewie.setP(5.0, i, 3.0)
        xValues.append(i)
        yValues.append(len(stewie.solve()))

    print("Plotting...")

    plt.plot(xValues, yValues, 'k')

    # Setup plot

    plt.xlabel("x")
    plt.ylabel('Number of roots')

    title = "Activity 7 - Plotting number of roots for p2"

    plt.title(title)

    fig = plt.gcf()
    fig.canvas.set_window_title(title)

    plt.show()

    print("PASS")

def test8():

    stewie = StewartPlatform3D(5, 1, 7)

    stewie.setB(7)

    stewie.setPlatformPosition(0, 0, 5)
    print(stewie.inverseKinematics())
    stewie.plotPlatform()

    stewie.setPlatformRotation(30, 0, 0)
    print(stewie.inverseKinematics())
    stewie.plotPlatform()

    stewie.setPlatformRotation(0, 30, 0)
    print(stewie.inverseKinematics())
    stewie.plotPlatform()

    stewie.setPlatformRotation(0, 0, 30)
    print(stewie.inverseKinematics())
    stewie.plotPlatform()

# Run Tests
test8()
#test1()
#test2()
#test3()
#test4()
#test5()
#test6()
#test7()




