import mceliece_full_scheme as mceliece
import saber_scheme as saber
import frodo

print("\nRunning some FRODO trials...")
frodo.random_trials()

print("\nRunning some SABER trials...")
saber.random_trials()

print("\nRunning some MCELIECE trials...")
mceliece.random_trials(num_trials=5)
