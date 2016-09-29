# Numerical Analysis 1
# Group project 1 - Stewart Platform
#
# Authors:
#   Cody Smith
#   William Rooney
#   David Andrews
#   Yunzhou Li

from math import *
import matplotlib.pyplot as plt
import numpy as np
from basic_units import radians

class StewartPlatform:

    def __init__(self, l1, l2, l3, gamma, x1, y1, x2, y2):

        # Set class members
        self.setL(l1, l2, l3)
        self.setX1(x1)
        self.setY1(y1)
        self.setX2(x2)
        self.setY2(y2)
        self.setGamma(gamma)

    ### Setters

    # Set x1
    def setX1(self, x1):
        self.x1 = x1

    # Set x2
    def setX2(self, x2):
        self.x2 = x2

    # Set y1
    def setY1(self, y1):
        self.y1 = y1

    # Set y2
    def setY2(self, y2):
        self.y2 = y2

    # Set side 1 length
    def setL1(self, l1):
        self.l1 = l1

    # Set side 2 length
    def setL2(self, l2):
        self.l2 = l2

    # Set side 3 length
    def setL3(self, l3):
        self.l3 = l3

    # Set all lengths (L1, L2, L3)
    def setL(self, l1, l2, l3):
        self.setL1(l1)
        self.setL2(l2)
        self.setL3(l3)

    # Set strut 1 length
    def setP1(self, p1):
        self.p1 = p1

    # Set strut 2 length
    def setP2(self, p2):
        self.p2 = p2

    # Set strut 3 length
    def setP3(self, p3):
        self.p3 = p3

    # Set all strut lengths
    def setP(self, p1, p2, p3):
        self.setP1(p1)
        self.setP2(p2)
        self.setP3(p3)

    def setGamma(self, gamma):
        self.gamma = gamma

    ### Solvers

    ### Public methods

    # Given L1, L2 and L3 as well as P1, P2 and P3, compute x, y and theta
    def forwardKinematics(self, p1, p2, p3):

        # Multiple solutions?

        return

    # Given x, y and theta, compute p1, p2 and p3
    def inverseKinematics(self, x, y, theta):

        # Run the squared methods
        P1S = self.P1S(x, y)
        P2S = self.P2S(x, y, theta)
        P3S = self.P3S(x, y, theta)

        # Take the square roots of those values
        P1 = sqrt(P1S)
        P2 = sqrt(P2S)
        P3 = sqrt(P3S)

        return P1, P2, P3

    ### Private methods

    def P1S(self, x, y):

        # p1 = x^2 + y^2

        return pow(x, 2) + pow(y, 2)

    def P2S(self, x, y, theta):

        # p2^2 = (x + A2)^2 + (y + B2)^2

        # Used for solving x and y
        # p2^2 = x^2 + y^2 + 2 A2 x + 2 B2 y + A2^2 + B2^2
        # p2^2 = p1^2 + 2 A2 x + 2 B2 y + A2^2 + B2^2

        return pow(x + self.A2(theta), 2) + pow(y + self.B2(theta), 2)

    def P3S(self, x, y, theta):

        # p3 = (x + A3)^2 + (y + B3)^2

        # Used for solving x and y
        # p3^2 = x^2 + y^2 + 2 A3 x + 2 B3 y + A3^2 + B3^2
        # p3^2 = = p1^2 + 2 A3 x + 2 B3 y + A3^2 + B3^2

        return pow(x + self.A3(theta), 2) + pow(y + self.B3(theta), 2)

    def A2(self, theta):

        # A2 = L3 cos theta - x1

        return self.l3 * cos(theta) - self.x1

    def A3(self, theta):

        # A3 = L2 cos(theta + gamma) - x2

        return self.l2 * cos(theta + self.gamma) - self.x2

    def B2(self, theta):

        # B2 = L3 sin theta

        return self.l3 * sin(theta)

    def B3(self, theta):

        # B3 = L2 sin(theta + gamma) - y2

        return self.l2 * sin(theta + self.gamma) - self.y2

    def D(self, theta):

        # D = 2(A2 * B3 - B2 * A3)

        return 2 * (self.A2(theta) * self.B3(theta) - self.B2(theta) * self.A3(theta))

    def N1(self, x, y, theta):

        # N1 = B3(P2^2 - P1^2 - A2^2 - B2^2) - B2(P3^2 - P1^2 - A3^2 - B3^2)

        # Store results for A
        A2 = self.A2(theta)
        A2S = pow(A2, 2)

        A3 = self.A3(theta)
        A3S = pow(A3, 2)

        # Store results for B
        B2 = self.B2(theta)
        B2S = pow(B2, 2)

        B3 = self.B3(theta)
        B3S = pow(B3, 2)

        # Store results for P
        P1S = self.P1S(x, y)
        P2S = self.P2S(x, y, theta)
        P3S = self.P3S(x, y, theta)

        return B3(P2S - P1S - A2S - B2S) - B2(P3S - P1S - A3S - B3S)

    def N2(self, x, y, theta):

        # N2 = -A3(P2^2 - P1^2 - A2^2 - B2^2) + A2(P3^2 - P1^2 - A3^2 - B3^2)

        # Store results for A
        A2 = self.A2(theta)
        A2S = pow(A2, 2)

        A3 = self.A3(theta)
        A3S = pow(A3, 2)

        # Store results for B
        B2 = self.B2(theta)
        B2S = pow(B2, 2)

        B3 = self.B3(theta)
        B3S = pow(B3, 2)

        # Store results for P
        P1S = self.P1S(x, y)
        P2S = self.P2S(x, y, theta)
        P3S = self.P3S(x, y, theta)

        return (-A3)(P2S - P1S - A2S - B2S) + A2(P3S - P1S - A3S - B3S)

    def x(self, x, y, theta):

        # x = N1 / D
        # As long as D =/= 0
        # TODO: Implement zero checking

        return self.N1(x, y, theta) / self.D(theta)

    def y(self, x, y, theta):

        # y = N2 / D
        # As long as D =/= 0
        # TODO: Implement zero checking

        return self.N2(x, y, theta) / self.D(theta)

    # Activity #1 - Evaluate f(theta)
    def f(self, theta):

        # Store results for A
        A2 = self.A2(theta)
        A2S = pow(A2, 2)

        A3 = self.A3(theta)
        A3S = pow(A3, 2)

        # Store results for B
        B2 = self.B2(theta)
        B2S = pow(B2, 2)

        B3 = self.B3(theta)
        B3S = pow(B3, 2)

        # Store results for P
        P1S = pow(self.p1, 2)
        P2S = pow(self.p2, 2)
        P3S = pow(self.p3, 2)

        # Store results for N
        N1 = B3 * (P2S - P1S - A2S - B2S) - B2 * (P3S - P1S - A3S - B3S)
        N1S = pow(N1, 2)

        N2 = (-A3) * (P2S - P1S - A2S - B2S) + A2 * (P3S - P1S - A3S - B3S)
        N2S = pow(N2, 2)

        # Store results for D
        DS = pow(self.D(theta), 2)

        # Compute f
        return N1S + N2S - P1S * DS

    # Activity #2 - Plot f(theta) from -PI to PI
    def plotF(self, min, max):

        x = []
        y = []

        for i in np.arange(min, max, 0.001):
            x.append(i*radians)
            y.append(self.f(i))

        plt.plot(x, y, xunits = radians)
        plt.ylabel('f(theta)')
        plt.xlabel("theta (Radians)")
        plt.show()