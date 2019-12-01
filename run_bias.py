import check_bias_new19nov as bias
import numpy as np
import trigger

path='/vol/astro3/lofar/sim/pipeline/events/'


'''
for example:
event number: 196034009
zenith=90-63.97
azimuth=307.92
core x=-105.58
core y=-46.23
'''

coreas_path='/vol/astro3/lofar/sim/pipeline/events/196034009/0/coreas/'
zenith=26.0*np.pi/180
azimuth=307.9*np.pi/180
xcore=-105.58
ycore=-46.23

## selection of events with different trigger criteria

#event=196034009
event=93990550
#event=95749948
#event=87892283
#event=272981612
#event=174634099

working_detectors, global_trigger_condition, local_trigger_condition, trigger_type=trigger.find_trigger(event)

# working detectors= array of 0/1 for each scintillator based on daily lora counts
# trigger type = 'd' for detector (#/20) or 's' for station (#/5)
# global trigger condition= #/20 detectors or #/5 stations
# local trigger contition = lora station condition.  normally 3/4, except some instances where it was changed to 2/4 if one detector was broken

print '_____________________________'
print trigger_type
print working_detectors
print global_trigger_condition
print local_trigger_condition
print '_____________________________'


ntrials=500
min_percentage_trigger=95.0



(bias_passed, min_chance_of_hit, chance_of_hit_all_sims) = bias.find_bias(coreas_path, zenith, azimuth, xcore, ycore, ntrials, working_detectors,trigger_type,local_trigger_condition=local_trigger_condition,min_percentage_trigger=min_percentage_trigger, trigger_condition_nof_detectors=global_trigger_condition)


print bias_passed
print min_chance_of_hit
print chance_of_hit_all_sims
