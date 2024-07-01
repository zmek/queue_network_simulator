
import sys
sys.path.insert(0, 'functions')


import queue_network_simulator as qns
import plotting_functions as pf



analysisID = "clinic"
analysis_location = "data"

number_of_runs = 20
warm_up_duration = 100
sim_duration = 60
plot_length = 60


model_parameters = qns.read_input(analysisID,analysis_location)


raw_data = {}


for r in range(number_of_runs):
    print (f"Run {r+1} of {number_of_runs}")
    my_model = qns.PathwayModel(r, model_parameters, warm_up_duration, sim_duration)
    my_model.run(r)
    raw_data[r+1] = my_model.output
    print()
    

plot_from = warm_up_duration + 1
plot_to = warm_up_duration + plot_length

time_transformation = 1


pf.plot_time_in_system(raw_data, plot_from, plot_to, analysisID, analysis_location, time_transformation)

for node_iter in my_model.node_ids:
    pf.plot_waiting_time(raw_data, node_iter, plot_from, plot_to, analysisID, analysis_location, time_transformation)

pf.plot_wait_over_time(raw_data, plot_from, plot_to, analysisID, analysis_location, time_transformation)