import os
import sys
import hashlib

n_iters = sys.argv[1]
n_iters = int(n_iters)

with open("msg.txt", 'w') as xf:
    with open("digest.txt", 'w') as yf:
        for i in range(n_iters):
            # msg
            plain = os.urandom(32)
            xh = plain.hex()
            xf.write(xh)
            xf.write('\n')
            # cipher
            h = hashlib.new('sha256')
            h.update(plain)
            yf.write(h.hexdigest())
            yf.write('\n')
