#imports 
import apogee_tools as ap
import numpy as np
import emcee
import matplotlib.pyplot as plt
import corner
import random

#toy_model, this should be replaced by apogee_tools model class
class ToyModel(object):

    def __init__(self, **kwargs):
        theta=[0.0, 0.0, 0.0, 0.0, 0.0]
        chi=0.0
        
    @property
    def theta(self):
        #this is a random model should be replaced by the actual labels
        #must be arranged as a list [apparently EMCEE likes working with lists]
        # as [teff, logg,  fe_h, rv, vsini]
        return [random.uniform(150, 4000.0),
               random.uniform(0.0,6.0), 
               random.uniform(-3.0, 3.0),
               random.uniform(-50.0,  50.0),
               random.uniform(-100, 100.0)]
    @property
    def chi(self):
        #chi-square
        #this is a random chi-square should be replaced by the fitting code
        return random.uniform(0.0, 10000.0)

def lnlike(params):
    """
    likelihood function for our MCMC
    theta: should be able to generate a model from parameters (notice I don't use this because
    I return a random model) 
    chi: chi-sqaure
    """
    model = ToyModel() #replace this by actual getmodel from apogee_tools as model=ap.get_model(params)
    lnp = lnprior(model)
    lnposterior = -0.5*model.chi
    
    return lnp+lnposterior
    
def lnprior(model):
    """
    This function uses a flat prior by constraining parameters between some range
    """
    teff, logg,  fe_h, rv, vsini = model.theta  #parameters should be arranged in this manner
    
    if  (prior["teff"][0] < teff < 4000) and (3.0 < logg  < 6.0) \
    and (-2.0 < fe_h < 2.0) and (-50.0< rv < 50.0) \
    and (0 < vsini < 20):
        return 0.0
    return -np.inf

def run_MCMC(**kwargs):
    """
    MCMC run using EMCEE
    params: model (must have a chi-square and a list of parameters)
    
    optional params: 
    first guess: guess
    nwalkers: number of walkers
    nsamples: number of runs default: 10000
    
    """
    #model to use for initial guess if not given
    guess_model=kwargs.get('guess_model', ToyModel())
    guess= kwargs.get('guess', guess_model.theta)
    params=np.array(guess)
    
    #noise=kwargs.get('noise', spectrum.noise.value)
    ndim=len(model.theta)
    nwalkers=kwargs.get('nwalkers', 10)
    nsamples=kwargs.get('nsamples', 1000)
    
    #make the first guess the same for all walkers with some random number added
    p0=[guess for n in range(0,nwalkers)]
    for n in range(1, nwalkers):
        p0[n]=[p+ np.random.random(1)[0]*0.0001*p for p in params]
    
    print ("walkers {} initial guess for each walker {} samples {}".format(nwalkers, p0, nsamples))

    
    #make this a numpy array
    p0= np.array(p0)
        
    #emcee sampler object
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnlike)
    
    #run the MCMC 
    ps, lnps, rstate= sampler.run_mcmc(p0,nsamples)
    
    #get the posteriors
    samples = sampler.chain.reshape((-1, ndim))
    
    #get mean parameters
    pmean=[np.nanmean((samples.T)[i]) for i in range(0, ndim)]
    #get standard deviations
    pstd=[np.nanstd((samples.T)[i]) for i in range(0, ndim)]
    
    #pf= [np.mean(((samples.T)[i])[-10:]) for i in range(0, ndim)]
    
    #visualization 
    param_names=['teff', 'logg',  'fe_h', 'rv', 'vsini']
    #std keys
    param_er_names=['teff_er', 'logg_er',  'fe_h_er', 'rv_er', 'vsini_er']
    
      
    if kwargs.get('show_corner', True):
        plt.figure()
        fig = corner.corner(samples, labels=param_names)
        plt.show()
        
        plt.figure()
        fig, ax = plt.subplots(ndim, sharex=True, figsize=(12, 6))
        for i in range(ndim):
            ax[i].plot(sampler.chain[:, :, i].T, '-k', alpha=0.2)
        fig.show()
    
    return samples, dict(zip(param_names, pmean)), dict(zip(param_er_names, pstd))
    
    
if __name__=='__main__':
	model=ToyModel()
	run_MCMC(guess_model=model)