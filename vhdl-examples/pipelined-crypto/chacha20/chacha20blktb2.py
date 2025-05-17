import os
import sys
from Cryptodome.Cipher import ChaCha20

n_iters = sys.argv[1]
n_iters = int(n_iters)

with open("key-nonce.txt", 'w') as xf:
    with open("state.txt", 'w') as yf:
        for i in range(n_iters):
            # key
            key = os.urandom(32)
            xh = key.hex()
            xf.write(xh)
            xf.write('\n')
            # nonce
            nonce = os.urandom(8)
            xh = nonce.hex()
            xf.write(xh)
            xf.write('\n')
            # cipher
            cipher = ChaCha20.new(key=key, nonce=nonce)
            y = cipher.encrypt(bytes([0 for i in range(512//8)]))
            yh = y.hex()
            yf.write(yh)
            yf.write('\n')
