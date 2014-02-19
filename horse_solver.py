#--------------------------------------------------------------
# This script is an attempt to solve the logic problem.
#               HORSE
#            + SADDLE
#            --------
#              GALLOP
#--------------------------------------------------------------

import itertools

def createNum(myTup):
    myJoin = ''.join(str(x) for x in myTup)
    return(int(myJoin))

if __name__ == '__main__':
    nums = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    for permute in itertools.permutations(nums):
        A, D, E, G, H, L, O, P, R, S = permute

        if S != 0 and H != 0:
            horse = (H,O,R,S,E)
            saddle = (S,A,D,D,L,E)
            gallop = (G,A,L,L,O,P)

            hNum = createNum(horse)
            sNum = createNum(saddle)
            gNum = createNum(gallop)

            if hNum + sNum == gNum:
                print "A %s, D %s, E %s, G %s, H %s, L %s, O %s, P %s, R %s, S %s" % (A, D, E, G, H, L, O, P, R, S)
                print "Saddle equals: %s" % sNum
                print "Horse equals: %s" % hNum
                print "Gallop equals: %s" % gNum

