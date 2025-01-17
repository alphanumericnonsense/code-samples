import os
import sys
import hashlib

n_iters = int(sys.argv[1])
width = int(sys.argv[2])

with open("msg.txt", 'w') as xf:
    with open("digest.txt", 'w') as yf:
        for i in range(n_iters):
            # msg
            plain = os.urandom(width//8)
            xh = plain.hex()
            xf.write(xh)
            xf.write('\n')
            # cipher
            h = hashlib.new('sha3_256')
            h.update(plain)
            yf.write(h.hexdigest())
            yf.write('\n')
