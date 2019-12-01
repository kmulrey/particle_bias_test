import check_bias as bias
import numpy as np

path='/vol/astro3/lofar/sim/pipeline/events/'


'''
    for example:
    event number: 196034009
    zenith=90-63.97
    azimuth=307.92
    core x=-105.58
    core y=-46.23
    '''
event=203459514
#event=203


# information that comes from fit routine
coreas_path='/vol/astro3/lofar/sim/pipeline/events/'+str(event)+'/0/coreas/'
zenith=26.0*np.pi/180
azimuth=307.9*np.pi/180
xcore=-105.58
ycore=-46.23



ntrials=500
min_percentage_trigger=95.0



(bias_passed, min_chance_of_hit, chance_of_hit_all_sims) = bias.find_bias(event, coreas_path, zenith, azimuth, xcore, ycore, ntrials,min_percentage_trigger=min_percentage_trigger)


print bias_passed
print min_chance_of_hit
print chance_of_hit_all_sims

