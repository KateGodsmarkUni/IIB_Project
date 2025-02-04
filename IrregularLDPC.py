import random
import tempfile
import os
import subprocess

from libs.TannerGraph import *

'''
A class for the handling of Regular LDPC matrices in tanner graph form

The tanner graph is stored as a dictionary, row indices (check nodes) are mapped to lists of column indices (variable
nodes) to indicate bipartite connections. Although this class defines regular matrices, it is not a requirement that
row and column weightages be constant. This is attempted in the respective constructions, but is not always possible
given the following premise for completely regular codes:
h = w * (c/r) where h = height, w = width, c = column weightage, r = row weightage

args: input to enable construction. input follows the following construction patern: h = n * (c / r) where c / r is not
explicitly simplified

construction specifies the method by which the matrix corresponding to args will be created
'''


class IrregularLDPC(TannerGraph):

    # parameters:
    #   args: list of arguments needed for irregular ldpc construction
    #   construction: the type of construction to be used
    # return:
    #   a fully defined Irregular LDPC code
    def __init__(self, args, construction, verbose=False):
        TannerGraph.__init__(self, args, construction=construction)

        self.width = int(self.args[0])

        #
        # args provided [width (n_bits), height (n_checks), degree distribution]
        #
        # Because r is dependent on width, height, and c (assuming regularity), defining c results
        # in limiting the constructor to one possible r value.
        # The resulting n, c, r values are passed to the constructor.
        #
        if len(self.args) == 3:
            self.height = int(self.args[1])
            self.n = int(self.args[0])
            self.deg_dist = self.args[2]
            # self.r = int((self.width / self.height) * self.c)
        else:
            raise RuntimeError("invalid input provided")

        self.tanner_graph = \
            IrregularLDPC.get_parity_check_graph(self.n, self.height, self.deg_dist, self.construction)

    # parameters:
    #   n: int, the width of the LDPC code, the codeword length
    #   height: int, the height of the LDPC code, the number of check nodes
    #   r: int, the weight of each row of the LDPC code
    #   c: int, the weight of each column of the code
    #   method: String, the construction to be employed
    @staticmethod
    def get_parity_check_graph(n, m, dist, method):

        if method == "peg":
            # use PEG (progressive edge growth) algorithm
            # we will use the peg/ library for this purpose

            # first step is to create a degree distribution file
            # which is trivial in our case since we have regular codes
            if len(dist) % 2 != 0:
                raise RuntimeError(f"Degree distrbution must be list with even number of entries, got {dist}")
            with tempfile.NamedTemporaryFile(mode='w', delete=False, dir=os.getcwd()) as f:
                degFileName = f.name
                num_degrees = int(len(dist) / 2)
                degs = []
                weights = []
                i = 0
                while i < len(dist):
                    if dist[i] % 1 != 0:
                        raise RuntimeError("Degrees must be integers")
                    else:
                        dist[i] = int(dist[i])
                    degs.append(dist[i])
                    weights.append(round(dist[i+1], 5))
                    i += 2
                assert len(degs) == num_degrees
                assert len(weights) == num_degrees
                f.write(f'{num_degrees}\n')  # number of degrees
                for deg in degs:
                    f.write(str(deg) + '     ')  # degree of variable node
                f.write('\n')
                for w in weights:
                    f.write(str(w) + ' ')  # probability of the degree occurring
                f.write('\n')

            peg_library_path = os.path.join(os.path.dirname(__file__),
                                            os.path.join(os.path.pardir, os.path.pardir), 'peg')
            peg_exec_path = os.path.join(peg_library_path, 'MainPEG')

            # now create temporary file name to be used for output of peg
            with tempfile.NamedTemporaryFile(mode='w', delete=False, dir=os.getcwd()) as f:
                outFileName = f.name

            # run PEG
            # create seed first
            peg_seed = random.randrange(10000000)
            print(outFileName)
            subprocess.run(peg_exec_path + ' -numM ' + str(m) + ' -numN ' + str(n) +
                           ' -codeName ' + outFileName + ' -degFileName ' + degFileName +
                           ' -q -seed ' + str(peg_seed), shell=True)

            # os.remove(degFileName)
            # create the initial empty graph
            tanner_graph = {}
            for i in range(m):
                tanner_graph[i] = []

            # now read the graph generated by PEG
            with open(outFileName) as f:
                # verify first two lines are n and m respectively
                assert n == int(f.readline().rstrip('\n'))
                assert m == int(f.readline().rstrip('\n'))
                _ = f.readline()
                # ignore third line which contains number of columns for remaining lines
                check_num = 0
                for line in f:
                    vals = [int(val) for val in line.rstrip('\n').rstrip(' ').split(' ')]
                    for val in vals:
                        if val != 0:  # 0 is used to denote absence of variable node
                            tanner_graph.get(check_num).append(val - 1)
                    check_num += 1
                assert check_num == m
            os.remove(outFileName)
            return tanner_graph
        else:
            raise RuntimeError('Invalid construction method')
