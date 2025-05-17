import os
import sys
from Cryptodome.Cipher import AES

n_iters = sys.argv[1]
n_iters = int(n_iters)

with open("key-plain.txt", 'w') as xf:
    with open("cipher.txt", 'w') as yf:
        for i in range(n_iters):
            # key
            key = os.urandom(16)
            xh = key.hex()
            xf.write(xh)
            xf.write('\n')
            # plain
            plain = os.urandom(16)
            xh = plain.hex()
            xf.write(xh)
            xf.write('\n')
            # cipher
            cipher = AES.new(key, AES.MODE_ECB)
            y = cipher.encrypt(plain)
            yh = y.hex()
            yf.write(yh)
            yf.write('\n')
