import apogee_tools as ap
import os

if __name__ == "__main__":

#psuedo work flow:
	# init = ap.initialize()
	# fit = ap.fitMCMC(init)

	init_par, step_par, fiber = ap.initialize()

	print(fiber)