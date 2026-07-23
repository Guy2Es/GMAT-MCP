import load_gmat  # runs the bootstrap (Setup, etc.)
from load_gmat import gmat  # directly retrieves the real gmatpy module

print("GMAT loaded successfully")
sat = gmat.Construct("Spacecraft", "TestSat")
print("Object created:", sat.GetName())