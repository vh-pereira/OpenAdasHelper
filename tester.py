import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

from src.main import OpenADAS


adas = OpenADAS()

ne, Te = 1.19e19, 180.65

get_ExcPEC = adas.impact_excitation_pec('H',0,(3,2)) # ADAS API to ExcPEC
get_RecPEC = adas.recombination_pec('H',0,(3,2)) # ADAS API to RecPEC

print(get_ExcPEC.evaluate(ne,Te))
print(get_RecPEC.evaluate(ne, Te))