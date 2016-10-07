# Numerical Analysis 1
# Group project 1 - Stewart Platform
#
# Authors:
#   Cody Smith
#   William Rooney
#   David Andrews
#   Yunzhou Li

from basic_units import radians
from scipy.optimize import brentq, fsolve
from numpy import sin, cos, pi, sqrt, deg2rad

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import numpy as np

class StewartPlatform2D:

    def __init__(self, l1, l2, l3, gamma, x1, y1, x2, y2):

        # Set class members
        self.setL(l1, l2, l3)
        self.setX1(x1)
        self.setY1(y1)
        self.setX2(x2)
        self.setY2(y2)
        self.setGamma(gamma)
        self.debug = 0

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

    # Set theta
    def setTheta(self, theta):
        self.theta = theta

    # Set gamma
    def setGamma(self, gamma):
        self.gamma = gamma

    ### Private methods

    def P1S(self, x, y):

        # p1 = x^2 + y^2

        return pow(x, 2) + pow(y, 2)

    def P2S(self, x, y, theta = None):

        if theta is None:
            theta = self.theta

        # p2^2 = (x + A2)^2 + (y + B2)^2

        # Used for solving x and y
        # p2^2 = x^2 + y^2 + 2 A2 x + 2 B2 y + A2^2 + B2^2
        # p2^2 = p1^2 + 2 A2 x + 2 B2 y + A2^2 + B2^2

        return pow(x + self.A2(theta), 2) + pow(y + self.B2(theta), 2)

    def P3S(self, x, y, theta = None):

        if theta is None:
            theta = self.theta

        # p3 = (x + A3)^2 + (y + B3)^2

        # Used for solving x and y
        # p3^2 = x^2 + y^2 + 2 A3 x + 2 B3 y + A3^2 + B3^2
        # p3^2 = = p1^2 + 2 A3 x + 2 B3 y + A3^2 + B3^2

        return pow(x + self.A3(theta), 2) + pow(y + self.B3(theta), 2)

    def A2(self, theta = None):

        if theta is None:
            theta = self.theta

        # A2 = L3 cos theta - x1

        return self.l3 * cos(theta) - self.x1

    def A3(self, theta):

        if theta is None:
            theta = self.theta

        # A3 = L2 cos(theta + gamma) - x2

        return self.l2 * cos(theta + self.gamma) - self.x2

    def B2(self, theta = None):

        if theta is None:
            theta = self.theta

        # B2 = L3 sin theta

        return self.l3 * sin(theta)

    def B3(self, theta = None):

        if theta is None:
            theta = self.theta

        # B3 = L2 sin(theta + gamma) - y2

        return self.l2 * sin(theta + self.gamma) - self.y2

    def D(self, theta = None):

        if theta is None:
            theta = self.theta

        # D = 2(A2 * B3 - B2 * A3)

        return 2 * (self.A2(theta) * self.B3(theta) - self.B2(theta) * self.A3(theta))

    def N1(self, theta = None):

        if theta is None:
            theta = self.theta

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
        P1 = self.p1
        P2 = self.p2
        P3 = self.p3

        P1S = pow(P1, 2)
        P2S = pow(P1, 2)
        P3S = pow(P1, 2)

        return B3 * (P2S - P1S - A2S - B2S) - B2 * (P3S - P1S - A3S - B3S)

    def N2(self, theta = None):

        # N2 = -A3(P2^2 - P1^2 - A2^2 - B2^2) + A2(P3^2 - P1^2 - A3^2 - B3^2)

        if theta is None:
            theta = self.theta

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
        P1 = self.p1
        P2 = self.p2
        P3 = self.p3

        P1S = pow(P1, 2)
        P2S = pow(P1, 2)
        P3S = pow(P1, 2)

        return (-A3) * (P2S - P1S - A2S - B2S) + A2 * (P3S - P1S - A3S - B3S)

    def x(self, theta = None):

        # x = N1 / D
        # As long as D =/= 0
        # TODO: Implement zero checking

        if theta is None:
            theta = self.theta

        return self.N1(theta) / self.D(theta)

    def y(self, theta = None):

        # y = N2 / D
        # As long as D =/= 0

        if theta is None:
            theta = self.theta

        # TODO: Implement zero checking

        return self.N2(theta) / self.D(theta)

    ### Public methods

    # Activity #1 - Evaluate f(theta)
    def f(self, theta = None):

        if theta is None:
            theta = self.theta

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
        D = self.D(theta)
        DS = pow(D, 2)

        # Compute f
        return N1S + N2S - P1S * DS

    # Activity #2 - Plot f(theta) from -PI to PI
    def plotF(self, min, max, title = None):

        xValues = []
        yValues = []

        for i in np.arange(min, max, 0.001):
            xValues.append(i*radians)
            yValues.append(self.f(i))

        # Setup grid and units

        fig, ax = plt.subplots()

        ax.plot(xValues, yValues, xunits=radians)

        ax.grid(True, which='both')
        ax.axhline(y=0, color='k')
        ax.axvline(x=0, color='k')

        # Setup plot

        plt.xlabel("theta (Radians)")
        plt.ylabel('f(theta)')

        plt.title(title)

        fig = plt.gcf()
        fig.canvas.set_window_title(title)

        plt.show()

    # Activity #3 - Plot the triangle
    def plotTriangle(self, title = None):

        # Plot the triangle
        x = self.x()
        y = self.y()

        u1x = x
        u1y = y

        u2x = x + self.l3 * cos(self.theta)
        u2y = y + self.l3 * sin(self.theta)

        u3x = x + self.l2 * cos(self.theta + self.gamma)
        u3y = y + self.l2 * sin(self.theta + self.gamma)

        xValuesTriangle = [u1x, u2x, u3x, u1x]
        yValuesTriangle = [u1y, u2y, u3y, u1y]

        plt.plot(xValuesTriangle, yValuesTriangle, 'k')

        # Plot the anchor points

        xValuesAnchor = [0, self.x1, self.x2]
        yValuesAnchor = [0, self.y1, self.y2]

        plt.plot(xValuesAnchor, yValuesAnchor, 'ko')

        # Plot the struts

        plt.plot([0, u1x], [0, u1y], 'r--')
        plt.plot([self.x1, u2x], [self.y1, u2y], 'g--')
        plt.plot([self.x2, u3x], [self.y2, u3y], 'b--')

        # Setup plot

        plt.xlabel("x")
        plt.ylabel('y')

        plt.title(title)

        fig = plt.gcf()
        fig.canvas.set_window_title(title)

        plt.show()

        ### Solvers

    def solve(self):

        # From http://stackoverflow.com/questions/14878110/how-to-find-all-zeros-of-a-function-using-numpy-and-scipy

        U = np.linspace(-pi, pi, 500)
        c = self.f(U)
        s = np.sign(c)
        ret = []

        for i in range(500 - 1):

            if s[i] + s[i + 1] == 0:
                u = brentq(self.f, U[i], U[i + 1])
                z = self.f(u)

                if np.isnan(z) or abs(z) > 1e-3:
                    continue

                ret.append(u)

        return ret

class StewartPlatform3D:

    # Assume that the platform has legs of the same size
    # and that the platform is hexagonal

    def __init__(self, platformLegLength, strutMinLength, strutMaxLength):

        # Set class members

        self.setL(platformLegLength)

        self.setStrutMinLength(strutMinLength)
        self.setStrutMaxLength(strutMaxLength)

        # Set base points to default values

        self.setB(0)

        self.setPlatformPosition(0, 0, 0)
        self.setPlatformRotation(0, 0, 0)

        self.debug = 0

    ### Setters

    # Set min strut length
    def setStrutMinLength(self, min):
        self.min = min

    # Set max strut length
    def setStrutMaxLength(self, max):
        self.max = max

    # Set tbe position of the platform
    def setPlatformPosition(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    # Set Rotation of the platform
    def setPlatformRotation(self, pitch, roll, yaw):
        self.pitch = deg2rad(pitch)
        self.roll = deg2rad(roll)
        self.yaw = deg2rad(yaw)

    # Set base point 1
    def setB1(self, x, y, z):
        self.b1 = [x, y, z]

    # Set base point 2
    def setB2(self, x, y, z):
        self.b2 = [x, y, z]

    # Set base point 3
    def setB3(self, x, y, z):
        self.b3 = [x, y, z]

    # Set base point 4
    def setB4(self, x, y, z):
        self.b4 = [x, y, z]

    # Set base point 5
    def setB5(self, x, y, z):
        self.b5 = [x, y, z]

    # Set base point 6
    def setB6(self, x, y, z):
        self.b6 = [x, y, z]

    # Set platform length
    def setB(self, length):

        self.length = length

        # Generate the platform points based on length

        a = length
        b = length / 2
        c = a * sin(deg2rad(60.0))

        self.b1 = [a, 0, 0]
        self.b2 = [b, c, 0]
        self.b3 = [-b, c, 0]
        self.b4 = [-a, 0, 0]
        self.b5 = [-b, -c, 0]
        self.b6 = [b, -c, 0]

        # NOTE: These are sort of locally relative, will need to transform these per
        # the position and rotation values of the platform to get actual values

        #       p3     p2
        #         -----
        #       /       \
        #      /         \
        #  p4       +       p1
        #      \         /
        #       \       /
        #         -----
        #       p5     p6

    # Set platform length
    def setL(self, length):
        self.length = length

        # Generate the platform points based on length

        a = length
        b = length / 2
        c = a * sin(deg2rad(60.0))

        self.p1 = [a, 0, 0]
        self.p2 = [b, c, 0]
        self.p3 = [-b, c, 0]
        self.p4 = [-a, 0, 0]
        self.p5 = [-b, -c, 0]
        self.p6 = [b, -c, 0]

        # NOTE: These are sort of locally relative, will need to transform these per
        # the position and rotation values of the platform to get actual values

        #       p3     p2
        #         -----
        #       /       \
        #      /         \
        #  p4       +       p1
        #      \         /
        #       \       /
        #         -----
        #       p5     p6

    # Set strut 1 length
    def setP1(self, p1):
        self.p1 = p1

    # Set strut 2 length
    def setP2(self, p2):
        self.p2 = p2

    # Set strut 3 length
    def setP3(self, p3):
        self.p3 = p3

    # Set strut 4 length
    def setP4(self, p4):
        self.p4 = p4

    # Set strut 5 length
    def setP5(self, p5):
        self.p5 = p5

    # Set strut 6 length
    def setP6(self, p6):
        self.p6 = p6

    # Set all strut lengths
    def setP(self, p1, p2, p3, p4, p5, p6):
        self.setP1(p1)
        self.setP2(p2)
        self.setP3(p3)
        self.setP4(p4)
        self.setP5(p5)
        self.setP6(p6)

    ### Solvers

    ### Public methods

    # Given L, P1 ... P6, base position, etc. compute platform position and rotation
    def forwardKinematics(self):

        # Multiple solutions?

        return

    # Given L, base position, platform position, platform rotation, etc. compute P1 ... P6
    def inverseKinematics(self):

        PLAT_AUG_PTS = self.getPlatformAugmentedPoints()

        P1 = sqrt(
            pow(PLAT_AUG_PTS[0][0] - self.b1[0], 2) +
            pow(PLAT_AUG_PTS[0][1] - self.b1[1], 2) +
            pow(PLAT_AUG_PTS[0][2] - self.b1[2], 2))

        P2 = sqrt(
            pow(PLAT_AUG_PTS[1][0] - self.b2[0], 2) +
            pow(PLAT_AUG_PTS[1][1] - self.b2[1], 2) +
            pow(PLAT_AUG_PTS[1][2] - self.b2[2], 2))

        P3 = sqrt(
            pow(PLAT_AUG_PTS[2][0] - self.b3[0], 2) +
            pow(PLAT_AUG_PTS[2][1] - self.b3[1], 2) +
            pow(PLAT_AUG_PTS[2][2] - self.b3[2], 2))

        P4 = sqrt(
            pow(PLAT_AUG_PTS[3][0] - self.b4[0], 2) +
            pow(PLAT_AUG_PTS[3][1] - self.b4[1], 2) +
            pow(PLAT_AUG_PTS[3][2] - self.b4[2], 2))

        P5 = sqrt(
            pow(PLAT_AUG_PTS[4][0] - self.b5[0], 2) +
            pow(PLAT_AUG_PTS[4][1] - self.b5[1], 2) +
            pow(PLAT_AUG_PTS[4][2] - self.b5[2], 2))

        P6 = sqrt(
            pow(PLAT_AUG_PTS[5][0] - self.b6[0], 2) +
            pow(PLAT_AUG_PTS[5][1] - self.b6[1], 2) +
            pow(PLAT_AUG_PTS[5][2] - self.b6[2], 2))

        return [P1, P2, P3, P4, P5, P6]

    def plotPlatform(self):

        fig = plt.figure()

        ax = fig.add_subplot(111, projection = '3d')

        PLAT_AUG_PTS = self.getPlatformAugmentedPoints()

        # Plot struts

        ax.plot(
            [self.b1[0], PLAT_AUG_PTS[0][0]],
            [self.b1[1], PLAT_AUG_PTS[0][1]],
            [self.b1[2], PLAT_AUG_PTS[0][2]], 'b')
        ax.plot(
            [self.b2[0], PLAT_AUG_PTS[1][0]],
            [self.b2[1], PLAT_AUG_PTS[1][1]],
            [self.b2[2], PLAT_AUG_PTS[1][2]], 'b')
        ax.plot(
            [self.b3[0], PLAT_AUG_PTS[2][0]],
            [self.b3[1], PLAT_AUG_PTS[2][1]],
            [self.b3[2], PLAT_AUG_PTS[2][2]], 'b')
        ax.plot(
            [self.b4[0], PLAT_AUG_PTS[3][0]],
            [self.b4[1], PLAT_AUG_PTS[3][1]],
            [self.b4[2], PLAT_AUG_PTS[3][2]], 'b')
        ax.plot(
            [self.b5[0], PLAT_AUG_PTS[4][0]],
            [self.b5[1], PLAT_AUG_PTS[4][1]],
            [self.b5[2], PLAT_AUG_PTS[4][2]], 'b')
        ax.plot(
            [self.b6[0], PLAT_AUG_PTS[5][0]],
            [self.b6[1], PLAT_AUG_PTS[5][1]],
            [self.b6[2], PLAT_AUG_PTS[5][2]], 'b')

        # Plot strut endpoints

        ax.plot(
            [self.b1[0], PLAT_AUG_PTS[0][0]],
            [self.b1[1], PLAT_AUG_PTS[0][1]],
            [self.b1[2], PLAT_AUG_PTS[0][2]], 'ko')
        ax.plot(
            [self.b2[0], PLAT_AUG_PTS[1][0]],
            [self.b2[1], PLAT_AUG_PTS[1][1]],
            [self.b2[2], PLAT_AUG_PTS[1][2]], 'ko')
        ax.plot(
            [self.b3[0], PLAT_AUG_PTS[2][0]],
            [self.b3[1], PLAT_AUG_PTS[2][1]],
            [self.b3[2], PLAT_AUG_PTS[2][2]], 'ko')
        ax.plot(
            [self.b4[0], PLAT_AUG_PTS[3][0]],
            [self.b4[1], PLAT_AUG_PTS[3][1]],
            [self.b4[2], PLAT_AUG_PTS[3][2]], 'ko')
        ax.plot(
            [self.b5[0], PLAT_AUG_PTS[4][0]],
            [self.b5[1], PLAT_AUG_PTS[4][1]],
            [self.b5[2], PLAT_AUG_PTS[4][2]], 'ko')
        ax.plot(
            [self.b6[0], PLAT_AUG_PTS[5][0]],
            [self.b6[1], PLAT_AUG_PTS[5][1]],
            [self.b6[2], PLAT_AUG_PTS[5][2]], 'ko')

        # Plot Platform

        px = [PLAT_AUG_PTS[0][0], PLAT_AUG_PTS[1][0], PLAT_AUG_PTS[2][0], PLAT_AUG_PTS[3][0],
              PLAT_AUG_PTS[4][0], PLAT_AUG_PTS[5][0], PLAT_AUG_PTS[0][0]]
        py = [PLAT_AUG_PTS[0][1], PLAT_AUG_PTS[1][1], PLAT_AUG_PTS[2][1], PLAT_AUG_PTS[3][1],
              PLAT_AUG_PTS[4][1], PLAT_AUG_PTS[5][1], PLAT_AUG_PTS[0][1]]
        pz = [PLAT_AUG_PTS[0][2], PLAT_AUG_PTS[1][2], PLAT_AUG_PTS[2][2], PLAT_AUG_PTS[3][2],
              PLAT_AUG_PTS[4][2], PLAT_AUG_PTS[5][2], PLAT_AUG_PTS[0][2]]

        ax.plot(px, py, pz, 'ko')
        ax.plot(px, py, pz, 'k')

        plt.show()

    ### Private methods

    def getPlatformAugmentedPoints(self):
        C_ALPHA = cos(self.pitch)
        S_ALPHA = sin(self.pitch)

        C_BETA = cos(self.roll)
        S_BETA = sin(self.roll)

        C_GAMMA = cos(self.yaw)
        S_GAMMA = sin(self.yaw)

        ROT_POT_MTX = np.array([
            [C_GAMMA * C_BETA, C_GAMMA * S_BETA * S_ALPHA - S_GAMMA * C_ALPHA, C_GAMMA * S_BETA * C_ALPHA + S_GAMMA * S_ALPHA],
            [S_GAMMA * C_BETA, S_GAMMA * S_BETA * S_ALPHA + C_GAMMA * C_ALPHA, S_GAMMA * S_BETA * C_ALPHA - C_GAMMA * S_ALPHA],
            [-S_BETA, C_BETA * S_ALPHA, C_BETA * C_ALPHA]
        ])

        PLAT_PTS = np.array([self.p1, self.p2, self.p3, self.p4, self.p5, self.p6])

        A = np.array([self.x, self.y, self.z, self.pitch, self.roll, self.yaw]).transpose()

        X_NOT = A[0:3] - PLAT_PTS

        PLAT_AUG_PTS = np.zeros(PLAT_PTS.shape)

        # Compute rotation
        for i in xrange(6):
            PLAT_AUG_PTS[i, :] = np.dot(ROT_POT_MTX, PLAT_PTS[i, :])

        # Compute translation
        for i in xrange(6):
            PLAT_AUG_PTS[i, 0] = PLAT_AUG_PTS[i, 0] + self.x
            PLAT_AUG_PTS[i, 1] = PLAT_AUG_PTS[i, 1] + self.y
            PLAT_AUG_PTS[i, 2] = PLAT_AUG_PTS[i, 2] + self.z

        #print(PLAT_PTS)
        #print(PLAT_AUG_PTS)

        return PLAT_AUG_PTS