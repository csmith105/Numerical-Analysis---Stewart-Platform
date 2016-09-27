# Numerical Analysis 1
# Group project 1
#
# Authors:
#   Cody Smith
#   xxx
#   xxx
#   xxx

class StewartPlatform:

    def __init__(self, l1, l2, l3, x1, x2, y1):

        import math

        # Set class members
        self.setL(l1, l2, l3)
        self.setX1(x1)
        self.setX2(x2)
        self.setY1(y1)

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

    # Given x, y and theta, compute p1, p2 and p3
    def inverseKinematics(self, x, y, theta):

        # Output p1, p2, p3
    
    ### Private helper functions

    def __P1(self):

    def __P2(self):

    def __P3(self):

