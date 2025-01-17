from Cryptodome.Cipher import AES

key = bytes.fromhex("cafebabedeadbeefb01dfaceba5eba11")
plain = bytes(16)

refaes = AES.new(key, AES.MODE_ECB)
with open("gt.txt", 'w') as cipherfile:
    for i in range(100):
        cipher = refaes.encrypt(plain)
        cipherfile.write(cipher.hex())
        cipherfile.write('\n')
        plain = cipher
