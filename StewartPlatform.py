# Numerical Analysis 1
# Group project 1 - Stewart Platform
#
# Authors:
#   Cody Smith
#   William Rooney
#   David Andrews
#   Yunzhou Li

class StewartPlatform:

    def __init__(self, l1, l2, l3, x1, x2, y1, y2):

        # Set class members
        self.setL(l1, l2, l3)
        self.setX1(x1)
        self.setX2(x2)
        self.setY1(y1)
        self.setY2(y2)

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

    # Set length 1
    def setL1(self, l1):
        self.l1 = l1

    # Set length 2
    def setL2(self, l2):
        self.l2 = l2

    # Set length 3
    def setL3(self, l3):
        self.l3 = l3

    # Set all lengths (L1, L2, L3)
    def setL(self, l1, l2, l3):
        self.setL1(l1)
        self.setL2(l2)
        self.setL3(l3)

    ### Solvers

    # Given L1, L2 and L3 as well as P1, P2 and P3, compute x, y and theta for each P1, P2, P3
    def forwardKinematics(self, p1, p2, p3):

        # Output x, y, theta for p1, p2, p3

        # NOTE: Nine values need computed

        # Multiple solutions?
        return

    # Given x, y and theta, compute p1, p2 and p3
    def inverseKinematics(self, x, y, theta):

        # Output p1, p2, p3
        return

    ### Private helper functions

    def __P1Squared(self, x, y):

        # p1 = x^2 + y^2

        import math

        return math.pow(x, 2) + math.pow(y, 2)

    def __P2Squared(self, x, y, theta):

        # p2^2 = (x + A2)^2 + (y + B2)^2

        # Used for solving x and y
        # p2^2 = x^2 + y^2 + 2 A2 x + 2 B2 y + A2^2 + B2^2
        # p2^2 = p1^2 + 2 A2 x + 2 B2 y + A2^2 + B2^2

        import math

        return math.pow(x + self.__A2(theta), 2)\
               + math.pow(y + self.__B2(), 2)

    def __P3Squared(self, x, y, theta, gamma):

        # p3 = (x + A3)^2 + (y + B3)^2

        # Used for solving x and y
        # p3^2 = x^2 + y^2 + 2 A3 x + 2 B3 y + A3^2 + B3^2
        # p3^2 = = p1^2 + 2 A3 x + 2 B3 y + A3^2 + B3^2

        import math

        return math.pow(x + self.__A3(theta, gamma), 2)\
               + math.pow(y + self.__B3(theta, gamma), 2)

    def __A2(self, theta):

        # A2 = L3 cos θ − x1

        import math

        return self.l3() * math.cos(theta) - self.x1

    def __A3(self, theta, gamma):

        # A3 = L2 cos(θ + γ) − x2
        # A3 = L2[cos θ cos γ − sin θ sin γ] − x2

        import math

        return self.l2 * math.cos(theta - gamma) - self.x2

    def __B2(self, theta):

        # B2 = L3 sin θ

        import math

        return self.l3 * math.sin(theta)

    def __B3(self, theta, gamma):

        # B3 = L2 sin(θ + γ) − y2
        # B3 = L2[cos θ sin γ + sin θ cos γ] − y2

        import math

        return self.l2 * math.sin(theta + gamma) - self.y2

    def __D(self, theta, gamma):

        # D = 2(A2 B3 - B2 A3)

        import math

        return 2 * (self.__A2(theta) * self.__B3(theta, gamma)
                    - self.__B2(theta) * self.__A3(theta, gamma))

    def __N1(self, x, y, theta, gamma):

        # N1 = B3(p2^2 − p1^2 − A2^2 − B2^2) − B2(p3^2 − p1^2 − A3^2 − B3^2)

        import math

        A2 = self.__A2(theta)
        A2S = math.pow(A2, 2)

        A3 = self.__A3(theta, gamma)
        A3S = math.pow(A3, 2)

        B2 = self.__B2(theta)
        B2S = math.pow(B2, 2)

        B3 = self.__B3(theta, gamma)
        B3S = math.pow(B3, 2)

        P1S = self.__P1Squared(x, y)
        P2S = self.__P2Squared(x, y, theta)
        P3S = self.__P3Squared(x, y, theta, gamma)

        return B3(P2S - P1S - A2S - A2S) - B2(P3S - P1S - A3S - B3S)

    def __N2(self, x, y, theta, gamma):

        # N2 = −A3(p2^2 − p1^2 − A2^2 − B2^2) + A2(p3^2 − p1^2 − A3^2 − B3^2)

        import math

        A2 = self.__A2(theta)
        A2S = math.pow(A2, 2)

        A3 = self.__A3(theta, gamma)
        A3S = math.pow(A3, 2)

        B2 = self.__B2(theta)
        B2S = math.pow(B2, 2)

        B3 = self.__B3(theta, gamma)
        B3S = math.pow(B3, 2)

        P1S = self.__P1Squared(x, y)
        P2S = self.__P2Squared(x, y, theta)
        P3S = self.__P3Squared(x, y, theta, gamma)

        return (-A3)(P2S - P1S - A2S - B2S) + A2(P3S - P1S - A3S - B3S)

    def __x(self, x, y, theta, gamma):

        # x = N1 / D
        # As long as D =/= 0

        return self.__N1(x, y, theta, gamma) / self.__D(theta, gamma)

    def __y(self, x, y, theta, gamma):

        # y = N2 / D
        # As long as D =/= 0

        return self.__N2(x, y, theta, gamma) / self.__D(theta, gamma)