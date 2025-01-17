aestxt: aestxt.c AES.h
	gcc -o aestxt AES.h aestxt.c

aestest: aestest.c AES.h
	gcc -o aestest AES.h aestest.c

aestxtexample: aestxt
	./aestxt key.txt plain.txt cipher e
	./aestxt key.txt decrypted.txt cipher d

clean:
	rm -f aestest aestxt
