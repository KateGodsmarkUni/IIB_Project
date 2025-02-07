import numpy as np
from decimal import *
import random
from collections import defaultdict
import re
from time import time
import sys
import os.path
import subprocess

import subprocess

def execute(cmd):
    popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
    # for line in popen.stdout: 
    #     print(line.decode(), end='')
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)

#
# Set the maximum recursion depth to 2000
sys.setrecursionlimit(20000)


class Watermark():
    def __init__(self, N, pipd=0.04):
        self.insmax = 100
        self.Ps = Decimal(0.00)
        self.Pi = Decimal(pipd) 
        self.Pd = Decimal(pipd)
        self.Pt = Decimal(1) - (self.Pi + self.Pd)
        self.limiting = True
        self.N = N
        self.s = None
        self.info = None
        self.w = self.make_watermark(0.5)
        self.r = None
        self.Nr = None
        self.F_store  = defaultdict(dict)
        self.B_store  = defaultdict(dict)
        self.priors = None
        self.using_edge = False # whether to use the rigt-most column of the trellis
        self.data = None

    def reset_store(self):
        # reset forward and backward probabilities 
        self.F_store = defaultdict(dict)
        self.B_store = defaultdict(dict)

    def over_limit(self, i, j):
        # check if location within 5 sd of diagonal
        if abs(i - j) < 10:
            return False
        # max_drift = 5 * i**0.5 / 3 # set for pi = pd = 0.1
        max_drift = 3 * i **0.5 / 2.1 # for pi=pd=0.04, found through simulations
        if abs(i - j) > max_drift:
            return True
        else:
            return False
        
    def make_watermark(self, d):
        w = []
        for i in range(self.N):
            r = random.random()
            w.append(1) if r < d else w.append(0)
        assert len(w) == self.N
        return w

    def sparse(self, d):
        # sparsify data to have length N
        s = []
        info = []
        f = len(d) / self.N
        for i in range(len(d)):
            while random.random() > f:
                s.append(0)
            s.append(d[i])
            info.append(len(s) - 1)
        while len(s) != self.N:
            if len(s) < self.N:
                s.append(0)
            else:
                s, info = self.sparse(d)
        return s, info
    
    def add2(self, u, v):
        assert len(u) == len(v)
        sum = []
        for i in range(len(v)):
            sum.append((u[i] + v[i]) % 2)
        return sum
    
    def channel(self, data):
        id_events = 0
        r = []
        p = random.random()
        i = 0
        ins = 0
        while i < len(data):
            if self.Pi < p < (self.Pi + self.Pd):
                i += 1
                id_events += 1
                ins = 0
            elif p > (self.Pi + self.Pd):
                if random.random() < self.Ps:
                    r.append(1 - data[i])
                else:
                    r.append(data[i])
                i += 1
                ins = 0
            else:
                if ins < self.insmax:
                    r.append(0) if random.random() < 0.5 else r.append(1)
                    ins += 1
                    id_events += 1
            p = random.random()
        print(f"Passed through channel with {id_events} id events")
        return r, id_events

    def send_data(self, data):
        # sparsify data, add to w and send through channel
        self.data = data
        assert len(data) < self.N
        self.s, self.info = self.sparse(data)
        t = self.add2(self.w, self.s)
        self.r, id_ev = self.channel(t)
        self.Nr = len(self.r)

    def compute_forward(self, i, j, logging=False):
        if i == 0 and j == 0:
            f = Decimal(1)
        elif i == 0:
            f = self.forward(i, j - 1)*self.gamma(i, j - 1, i, j)
        elif j == 0:
            f = self.forward(i - 1, j)*self.gamma(i - 1, j, i, j)
        elif i > self.N or j > self.Nr:
            f = 0
        elif i < 0 or j < 0:
            f = 0
        else:
            f = self.forward(i - 1, j)*self.gamma(i - 1, j, i, j) + self.forward(i, j - 1)*self.gamma(i, j - 1, i, j) + self.forward(i - 1, j - 1) * self.gamma(i - 1, j - 1, i, j)
        self.F_store[i][j] = f
        # print(f"Found forward of {i},{j} is {f}")
        if f == 0:
            f = Decimal(0)
        return f
    
    def compute_backward(self, i, j, logging=False):
        if i == self.N and j == self.Nr:
            b = Decimal(1)
        elif i == self.N:
            if self.using_edge:
                b = self.backward(i, j + 1, logging)*self.gamma(i, j, i, j + 1)
            # Control if we can end on insertions
            else:
                b = 0
        elif j == self.Nr:
            b = self.backward(i + 1, j)*self.gamma(i, j, i + 1, j)
        elif i > self.N or j > self.Nr:
            b = 0
        elif i < 0 or j < 0:
            b = 0
        else:
            b = self.backward(i + 1, j, logging)*self.gamma(i, j, i + 1, j) + self.backward(i, j + 1, logging)*self.gamma(i, j, i, j + 1) + self.backward(i + 1, j + 1, logging)*self.gamma(i , j , i + 1, j + 1)
        if b == 0:
            b = Decimal(0)
        self.B_store[i][j] = b
        return b
    
    def gamma(self, i1, j1, i2, j2):
        if i2 == i1 and j2 == j1 + 1:
            g = self.Pi / 2
        elif i2 == i1 + 1 and j2 == j1:
            g = self.Pd
        elif i2 == i1 + 1 and j2 == j1 + 1:
            if i1 not in self.info: # not a data bit, just watermark
            # print(f"Comparing rec{j1} and w{i1}")
                if self.r[j1] == self.w[i1]:
                    g = self.Pt*(1 - self.Ps)
                else:
                    g = self.Pt*self.Ps
            else:
                ind = self.info.index(i1) # Index in data string
                p = Decimal(self.priors[ind]) # prior probability that bit is 0
                
                if (self.r[j1] + self.w[i1]) % 2 == 1: #
                    g = p * self.Pt * self.Ps + (1 - p) * self.Pt * (1 - self.Ps)
                    # print(f"Probability that bit is zero is {p}, prob that became 1 is {g}")
                else:
                    g = p * self.Pt * (1 - self.Ps) + (1 - p) * self.Pt * self.Ps
                    # print(f"Probability that bit is zero is {p}, prob that became 0 is {g}")
        else:
            print("not a valid gamma")
        return Decimal(g)
    
    def forward(self, i, j, logging=False):
        if self.over_limit(i, j) and self.limiting:
            return 0
        if i in self.F_store:
            if j in self.F_store[i]:
                return self.F_store[i][j]
        return self.compute_forward(i, j, logging)
    
    def backward(self, i, j, logging=False):
        if self.over_limit(i, j) and self.limiting:
            return 0
        if i in self.B_store:
            if j in self.B_store[i]:
                return self.B_store[i][j]
        return self.compute_backward(i, j, logging)
    
    def likelihood(self, i, val):
        l = Decimal(0)
        for j in range(self.Nr):
            if self.over_limit(i, j) and self.limiting:
                continue
            if self.r[j] == (val + self.w[i]) % 2:
                diag = self.Pt*(1 - self.Ps)
            else:
                diag = self.Pt*self.Ps

            l2 = self.forward(i, j) * diag * self.backward(i + 1, j + 1)
            l += l2
    
            l3 = self.forward(i, j) * self.Pd * self.backward(i + 1, j)
            l += l3
  
        l += self.forward(i, self.Nr) * self.Pd * self.backward(i + 1, self.Nr) # top row of trellis
        return l
    
    def mix_strings(self, s1, s2, p):
        # Function to vary priors for testing
        assert len(s1) == len(s2)
        new = []
        for i in range(len(s1)):
            if random.random() < p:
                new.append(s1[i])
            else:
                new.append(s2[i])
        return new
    
    def bin_list_to_file(self, ls, file):
        with open(file, 'w') as f:
            for l in ls:
                f.write(str(l))
    
    def list_to_file(self, ls, file):
        with open(file, 'w') as f:
            for l in ls:
                f.write(str(l) + ' ')
    
    def c_decode(self, llr_file, prior_file=None, test_name="t1"):
        N_data = len(self.info)
        print(self.Nr)
        if not prior_file:
            self.list_to_file([0.5]*N_data, f"{test_name}_outp")
            prior_file = f"{test_name}_outp"
        self.list_to_file(self.r, f"{test_name}_r")
        self.list_to_file(self.info, f"{test_name}_info")
        self.list_to_file(self.w, f"{test_name}_w")
        self.list_to_file(self.data, f"{test_name}_d")
        water_path = os.path.join(os.path.dirname(__file__), 'WaterDecode')
        execute(water_path + ' ' + prior_file + ' ' + llr_file + ' ' + f"{test_name}_r" + ' ' + f"{test_name}_info" + \
                       ' ' + f"{test_name}_w" + ' ' + str(self.N) + ' ' + str(self.Nr) + ' ' + str(N_data) + ' ' + \
                          str(0.04) + ' ' + f"{test_name}_d")
    
    def decode(self, llr_file, prior_file=None):
        # Takes in priors file and writes llr file
        t1 = time()
        successes = 0
        llrs = []
        one_llrs = []
        zero_llrs = []
        if not prior_file:
            self.priors = [0.5] * len(self.data)
        else:
            with open(prior_file) as f:
                s = f.read()
                matches = re.findall('([0-9.]+)', s)
                assert len(matches) == len(self.data)
                self.priors = []
                for m in matches:
                    self.priors.append(1 - float(m))
        # print(f"PRIORS: {self.priors}")
        
        for j in range(len(self.data)):
            p_win = 0
            pr = [0, 0]
            for v in range(2):
                p_ans = self.likelihood(self.info[j], v)
                if p_ans == 0:
                # print(f"Got zero from info bit {j} at position {info[j]}")
                    pass
                pr[v] = p_ans
                if p_win <= p_ans:
                    p_win = p_ans
                    ans = v

            if pr[0] != 0 and pr[1] != 0:
                lr = pr[0] / pr[1]
                llr = float(lr.ln())
                llrs.append(llr)
                if self.data[j] == 0:
                    zero_llrs.append(llr)
                else:
                    one_llrs.append(llr)
            else:
                pass
            if ans == self.data[j]:
                successes += 1
        # print(f"LLRS: {llrs}")
        
        st = ''
        for cl in llrs:
            s = str(cl)
            st = st + s
            st = st + ' '
        with open(llr_file, 'w') as f:
            f.write(st)
            f.close()
        print(f"{successes}/{len(self.data)} successes, ran in {time() - t1} s")