import os.path
import subprocess
from libs.IrregularLDPC import IrregularLDPC

class LDPC():
    def __init__(self, name, deg_dist, n_bits, n_checks):
        self.pchk_file = f"{name}.pchk"
        self.gen_file = f"{name}.gen"
        self.n_bits = n_bits
        self.n_checks = n_checks
        self.ldpc_library_path = os.path.join(os.path.dirname(__file__), os.path.pardir, 'LDPC-codes')
        self.deg_dist = deg_dist
        if not os.path.isfile(self.gen_file) or not os.path.isfile(self.pchk_file):
            # make pcheck
            irregular_dimension_args = [self.n_bits, self.n_checks, self.deg_dist]
            ldpc_code = IrregularLDPC(irregular_dimension_args, 'peg',
                                verbose=True)  # setting verbose to true prints code info during construction
            self.write_graph_to_file(ldpc_code, self.pchk_file)
            # make generator
            
            ldpc_make_gen_path = os.path.join(self.ldpc_library_path, 'make-gen')
            subprocess.run(ldpc_make_gen_path + ' ' + self.pchk_file + ' ' + self.gen_file + ' ' + \
                       'sparse', shell=True)
    
    def intio_write(self, file, value):
        for i in range(3):
            b = value & 0xff

            bAsBinary = int(str(bin(b)), 2)
            binaryBtoBytes = bAsBinary.to_bytes(1, 'little')

            file.write(binaryBtoBytes)

            value >>= 8

        if value > 0:
            b = value
        else:
            b = (value + 256) % 256

        bAsBinary = int(str(bin(b)), 2)
        binaryBtoBytes = bAsBinary.to_bytes(1, 'little')
        file.write(binaryBtoBytes)
    
    def write_graph_to_file(self, ldpc_code, filepath):

        with open(filepath, "wb") as f:

            self.intio_write(f, (ord('P') << 8) + 0x80)

            self.intio_write(f, ldpc_code.height)
            self.intio_write(f, ldpc_code.width)

            for key in ldpc_code.tanner_graph:
                self.intio_write(f, -(key + 1))
                for value in sorted(ldpc_code.tanner_graph.get(key)):
                    self.intio_write(f, (value + 1))

            self.intio_write(f, 0)
        
    def encode(self, src_file, out_path):
        ldpc_encode_path = os.path.join(self.ldpc_library_path, 'encode')
        subprocess.run(ldpc_encode_path + ' ' + self.pchk_file + ' ' + self.gen_file +
                   ' ' + src_file + ' ' + out_path, shell=True)

    def decode(self, rec_file, out_path, max_iters=10, channel_type='misc', channel_value='0.0'):
        assert max_iters > 0
        p_file = out_path + 'p'
        ldpc_decode_path = os.path.join(self.ldpc_library_path, 'decode')
        subprocess.run(ldpc_decode_path + ' ' + self.pchk_file + ' ' + rec_file + ' ' + \
                       out_path + ' ' + p_file + ' ' + channel_type + ' ' + str(channel_value) + \
                       ' prprp ' + str(max_iters), shell=True)
        
    def extract(self, codeword_file, message_file):
        ldpc_extract_path = os.path.join(self.ldpc_library_path, 'extract')
        subprocess.run(ldpc_extract_path + ' ' + self.gen_file + ' ' + codeword_file + ' ' + \
                       message_file, shell=True)