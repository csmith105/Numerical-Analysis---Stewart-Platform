# Numerical Analysis 1
# Group project 1 - Stewart Platform
#
# Authors:
#   Cody Smith
#   William Rooney
#   David Andrews
#   Yunzhou Li

class StewartPlatform:

    def init(self, l1, l2, l3, gamma, x1, y1, x2, y2):

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

        import math

        # Run the squared methods
        P1S = self.P1S(x, y)
        P2S = self.P2S(x, y, theta)
        P3S = self.P3S(x, y, theta)

        # Take the square roots of those values
        P1 = math.sqrt(P1S)
        P2 = math.sqrt(P2S)
        P3 = math.sqrt(P3S)

        return P1, P2, P3

    ### Private methods

    def P1S(self, x, y):

        # p1 = x^2 + y^2

        import math

        return math.pow(x, 2) + math.pow(y, 2)

    def P2S(self, x, y, theta):

        # p2^2 = (x + A2)^2 + (y + B2)^2

        # Used for solving x and y
        # p2^2 = x^2 + y^2 + 2 A2 x + 2 B2 y + A2^2 + B2^2
        # p2^2 = p1^2 + 2 A2 x + 2 B2 y + A2^2 + B2^2

        import math

        return math.pow(x + self.A2(theta), 2) + math.pow(y + self.B2(theta), 2)

    def P3S(self, x, y, theta):

        # p3 = (x + A3)^2 + (y + B3)^2

        # Used for solving x and y
        # p3^2 = x^2 + y^2 + 2 A3 x + 2 B3 y + A3^2 + B3^2
        # p3^2 = = p1^2 + 2 A3 x + 2 B3 y + A3^2 + B3^2

        import math

        return math.pow(x + self.A3(theta), 2) + math.pow(y + self.B3(theta), 2)

    def A2(self, theta):

        # A2 = L3 cos θ − x1

        import math

        return self.l3() * math.cos(theta) - self.x1

    def A3(self, theta):

        # A3 = L2 cos(θ + γ) − x2
        # A3 = L2[cos θ cos γ − sin θ sin γ] − x2

        import math

        return self.l2 * math.cos(theta - self.gamma) - self.x2

    def B2(self, theta):

        # B2 = L3 sin θ

        import math

        return self.l3 * math.sin(theta)

    def B3(self, theta):

        # B3 = L2 sin(θ + γ) − y2
        # B3 = L2[cos θ sin γ + sin θ cos γ] − y2

        import math

        return self.l2 * math.sin(theta + self.gamma) - self.y2

    def D(self, theta):

        # D = 2(A2 * B3 - B2 * A3)

        import math

        return 2 * (self.A2(theta) * self.B3(theta) - self.B2(theta) * self.A3(theta))

    def N1(self, x, y, theta):

        # N1 = B3(p2^2 − p1^2 − A2^2 − B2^2) − B2(p3^2 − p1^2 − A3^2 − B3^2)

        import math

        # Store results for A
        A2 = self.A2(theta)
        A2S = math.pow(A2, 2)

        A3 = self.A3(theta)
        A3S = math.pow(A3, 2)

        # Store results for B
        B2 = self.B2(theta)
        B2S = math.pow(B2, 2)

        B3 = self.B3(theta)
        B3S = math.pow(B3, 2)

        # Store results for P
        P1S = self.P1S(x, y)
        P2S = self.P2S(x, y, theta)
        P3S = self.P3S(x, y, theta)

        return B3(P2S - P1S - A2S - B2S) - B2(P3S - P1S - A3S - B3S)

    def N2(self, x, y, theta):

        # N2 = −A3(p2^2 − p1^2 − A2^2 − B2^2) + A2(p3^2 − p1^2 − A3^2 − B3^2)

        import math

        # Store results for A
        A2 = self.A2(theta)
        A2S = math.pow(A2, 2)

        A3 = self.A3(theta)
        A3S = math.pow(A3, 2)

        # Store results for B
        B2 = self.B2(theta)
        B2S = math.pow(B2, 2)

        B3 = self.B3(theta)
        B3S = math.pow(B3, 2)

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

    # Activity #1
    def f(self, theta):

        import math

        # Store results for A
        A2 = self.A2(theta)
        A2S = math.pow(A2, 2)

        A3 = self.A3(theta)
        A3S = math.pow(A3, 2)

        # Store results for B
        B2 = self.B2(theta)
        B2S = math.pow(B2, 2)

        B3 = self.B3(theta)
        B3S = math.pow(B3, 2)

        # Store results for P
        P1 = self.p1
        P1S = math.pow(P1, 2)

        P2 = self.p2
        P2S = math.pow(P2, 2)

        P3 = self.p3
        P3S = math.pow(P3, 2)

        # Store results for N
        N1 = B3(P2S - P1S - A2S - B2S) - B2(P3S - P1S - A3S - B3S)
        N1S = math.pow(N1, 2)

        N2 = (-A3)(P2S - P1S - A2S - B2S) + A2(P3S - P1S - A3S - B3S)
        N2S = math.pow(N2, 2)

        # Store results for D
        D = self.D(theta)
        DS = math.pow(D, 2)

        # Compute f
        return N1S + N2S - P2S * DS;

    ### Tests

        # TODO: Implement activity 1 test