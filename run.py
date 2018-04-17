import apogee_tools as ap
import os

if __name__ == "__main__":

	if 'config.yaml' not in os.listdir():
		print('\nError: config.yaml not found in the current working directory. \
			Using default file found inside apogee_tools.\n')

	init_par, step_par, fiber = ap.initialize()

	print(init_par)
	print(step_par)
	print(fiber)

	# fit = ap.fitMCMC(init_par, step_par, fiber)
