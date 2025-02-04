from Watermark import Watermark
from LDPC import LDPC
import random

def generate_data(L, p=0.5):
    data = [0] * L
    for j in range(L):
        if random.random() <= p:
            data[j] = 1
        else:
            data[j] = 0
    return data

def list_to_file(ls, file):
    with open(file, 'w') as f:
        for l in ls:
            f.write(str(l))

def file_to_list(file):
    ls = []
    with open(file) as f:
        s = f.read()
    for bit in s:
        if bit != '\n':
            ls.append(int(bit))
    return(ls)

def main():
    test_name = "ft3"
    K = 1700
    N = 3000
    n_checks = N - K
    NW = 4000
    
    m_data = generate_data(K)
    l_code = LDPC("code24", [2, 0.5, 3, 0.5], N, n_checks)
    list_to_file(m_data, f"{test_name}_m_data")
    l_code.encode(f"{test_name}_m_data", f"{test_name}_d_data")
    d_data = file_to_list(f"{test_name}_d_data")
    
    w1 = Watermark(NW, pipd=0.04)
    w1.send_data(d_data)
    w1.decode(f"{test_name}_llrs")

    l_code.decode(f"{test_name}_llrs", f"{test_name}_out", max_iters=20)

    for i in range(20):
        w1.reset_store()
        w1.decode(f"{test_name}_llrs", prior_file=f"{test_name}_outp")
        l_code.decode(f"{test_name}_llrs", f"{test_name}_out", max_iters=20)
        if i > 0:
            if file_to_list(f"{test_name}_d_data") == file_to_list(f"{test_name}_out"):
                print("SUCCESS!!!")
                # l_code.extract(f"{test_name}_out", f"{test_name}_message")
                break





if __name__ == '__main__':
    main()
    