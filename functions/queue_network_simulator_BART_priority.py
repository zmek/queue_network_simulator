


import pandas as pd
import os
import simpy as sp
import numpy as np
import random as rnd
import math
import matplotlib.pyplot as mpl


class Params:
    
    analysisID = "" #name of the folder containing input and output folders of the current analysis
    
    analysis_location = "" #relative path from this Python script to the folder containing [analysisID]'s folder
    
    number_of_runs = 0 #ideally 100 or above
    warm_up_duration = 0 #days - warm-up period, to be set to 0 if we allow a non-empty system as starting state
    sim_duration = 0 #days - duration of the actual simulation after the warm-up period
    

class Node:
    
    def __init__(self, env, node_id, resource_name, capacity):
        self.id = node_id
        self.resource_name = resource_name
        
        self.resource = sp.PriorityResource(env,capacity=capacity)
        

# Class representing the patients
class Patient:
    
    def __init__(self, pat_id, traj_name, starting_node_id, arrival_time, time_joined_current_queue):
        self.id = pat_id
        self.traj_name = traj_name
        
        self.current_node_id = starting_node_id
        self.traj_step_iterator = 1
        
        self.arrival_time = arrival_time
        
        self.total_waiting_time = 0
        
        self.priority = 1
        
 


# Class representing the pathway model
class PathwayModel:
    
    def __init__(self, run_number, data, warm_up_duration, sim_duration):
    
        self.env = sp.Environment()
        self.run_number = run_number
        self.warm_up_duration = warm_up_duration
        self.sim_duration = sim_duration
        
        self.patient_counter = 0
        
        self.node_ids = data["node_ids"]
        self.node_names = data["node_names"]
        self.capacities = data["capacities"]
        self.node_cost_per_time = data["node_cost_per_time"]
        
        self.nodes = {}
        for key in self.node_ids:
            self.nodes[key] = Node(self.env, key, self.node_names[key], self.capacities[key])
       
        
        self.arrival_distribution = data["arrival_distribution"] 
        self.arrival_distr_parameters = float(data["arrival_distr_parameter"])
        
        self.trajectories = data["trajectories"]
        self.traj_prop = data["traj_prop"]
        
        self.output = None
        
            
    #function generating arrivals for each trajectory
    #arrival_process_parameters are a dictionary as they might differ depending on the underlying probability distribution
    def generate_arrivals(self, arrival_distribution, arrival_distr_parameter):
        
        if arrival_distr_parameter > 0:  #if rate is 0, then this function does not do anything
            
            while True:
                
                #determine the next arrival time and freeze until then
                if arrival_distribution == "poisson":
                    time_to_next_arrival = rnd.expovariate(arrival_distr_parameter)
                else: #assuming the inter-arrival times are deterministic
                    time_to_next_arrival = 1 / arrival_distr_parameter
                
                yield self.env.timeout(time_to_next_arrival)
                
                #create the new patient
                self.patient_counter += 1
                
                #sample patient trajectory
                traj_list = list(self.traj_prop.keys())
                traj_weights = list(self.traj_prop.values())
                
                new_pat_traj_name = rnd.choices(traj_list, weights = traj_weights)[0]
                new_pat_start_node = self.trajectories[new_pat_traj_name]["step"][0]
                
                p = Patient(self.patient_counter, new_pat_traj_name, new_pat_start_node, self.env.now, self.env.now)

                #start patient p's journey through the network
                self.env.process(self.patient_journey(p))
                
                
            
    
    def patient_journey(self, patient):
        
        go_on = True
        
        while go_on:
            
            patient.time_joined_current_queue = self.env.now
            queue_start = patient.time_joined_current_queue
            
            with self.nodes[patient.current_node_id].resource.request(priority=patient.priority) as req:
                
                #Freeze the function until the request can be met
                yield req
                
                queue_finish = self.env.now
                
                # Determine the time the patient will spend at the current node
                serv_distr = self.trajectories[patient.traj_name]["service_time_distribution"][patient.traj_step_iterator-1]
                serv_distr_par = [float(x) for x in str(self.trajectories[patient.traj_name]["length_of_need"][patient.traj_step_iterator-1]).split(';')]
                
                if serv_distr == "exponential":
                    serv_time = rnd.expovariate(1 / serv_distr_par[0])
                elif serv_distr == "lognormal":
                    obs_mean = serv_distr_par[0]
                    obs_sd = serv_distr_par[1]
                    mu = math.log(obs_mean/math.sqrt(1+obs_sd**2/obs_mean**2))
                    sigma = math.sqrt(math.log(1+obs_sd**2/obs_mean**2))
                    serv_time = np.random.lognormal( mu, sigma )
                else: #assuming the service times are deterministic
                    serv_time = serv_distr_par[0]
                
                # Freeze this function until that time has elapsed
                yield self.env.timeout(serv_time)
                
                #record time at which the patient moves to the next node
                move_on_time = self.env.now
                
                patient.total_waiting_time = patient.total_waiting_time + queue_finish - queue_start
                
                #record waiting time and update waiting time output
                new_wait_time_record = pd.DataFrame({
                    "patID":[patient.id],
                    "trajectory":[patient.traj_name],
                    "node":[self.nodes[patient.current_node_id].id],
                    "time_joined_queue":[queue_start],
                    "time_started_service":[queue_finish],
                    "waiting_time":[queue_finish - queue_start],
                    "time_moved_on":[move_on_time]
                    })
                self.output.waiting_time = pd.concat([self.output.waiting_time, new_wait_time_record])
                
                
                #identify next node in the pathway or terminate patient's journey
                if patient.traj_step_iterator < len(self.trajectories[patient.traj_name]["step"]):
                    if patient.current_node_id == "Barts":
                        patient.priority = 0
                    else:
                        patient.priority = 1
                        
                    next_node_id = self.trajectories[patient.traj_name]["step"][patient.traj_step_iterator]
                    patient.traj_step_iterator += 1
                    patient.current_node_id = next_node_id
                else:
                    departure_time = self.env.now
                    go_on = False
                    
                    new_time_in_system_record = pd.DataFrame({
                        "patID":[patient.id],
                        "trajectory":[patient.traj_name],
                        "arrival_time":[patient.arrival_time],
                        "departure_time":[departure_time],
                        "time_in_system":[departure_time - patient.arrival_time],
                        "total_waiting_time":[patient.total_waiting_time]
                        })
                    self.output.time_in_system = pd.concat([self.output.time_in_system, new_time_in_system_record])
                    
                
        
    def run(self, run_number):
        
        self.output = Output(run_number, self.node_ids, self.warm_up_duration + self.sim_duration )
        
        for key in list(self.output.resource_utilisation.keys()):
            self.output.resource_utilisation[key]["capacity"] = [self.capacities[key] for i in range(int(np.ceil(self.warm_up_duration + self.sim_duration + 2)))]
        
        
        # Start entity generators
        self.env.process(self.generate_arrivals(self.arrival_distribution, self.arrival_distr_parameters))
        
       
        # Run simulation
        self.env.run( until = self.warm_up_duration + self.sim_duration )
        
        



# Class representing and manipulating the output of the simulation runs
class Output:
    
    def __init__(self, run_number, node_ids, time_frames):
        
        self.run_number = run_number
        
        self.time_in_system = pd.DataFrame()
        self.time_in_system["patID"] = []
        self.time_in_system["trajectory"] = []
        self.time_in_system["arrival_time"] = []
        self.time_in_system["departure_time"] = []
        self.time_in_system["time_in_system"] = []
        self.time_in_system["total_waiting_time"] = []
        
        self.waiting_time = pd.DataFrame()
        self.waiting_time["patID"] = []
        self.waiting_time["trajectory"] = []
        self.waiting_time["node"] = []
        self.waiting_time["time_joined_queue"] = []
        self.waiting_time["time_started_service"] = []
        self.waiting_time["waiting_time"] = []
        self.waiting_time["time_moved_on"] = []
        
        self.resource_utilisation = {}
        for node_key in node_ids:
            self.resource_utilisation[node_key] = pd.DataFrame()
            self.resource_utilisation[node_key]["time"] = [i for i in range(time_frames+2)]
            self.resource_utilisation[node_key]["capacity"] = [0 for i in range(time_frames+2)]
            self.resource_utilisation[node_key]["resource_utilisation"] = [0.0 for i in range(time_frames+2)]




# Function reading input from Excel spreadsheet
def read_input(analysisID,analysis_location):
    
    file_location = analysis_location + "/" + analysisID + "/input"
    
    #Read and store simulation parameters
    
    res = {}
        
    #Read and store information about trajectories
    trajectories = {}
    
    file_list = os.listdir(file_location + "/trajectories")
    
    for tr in file_list:
        file_type = os.path.splitext(tr)[1]
        if file_type == ".txt":
            traj_name = os.path.splitext(tr)[0]
            trajectories[traj_name] = pd.read_csv(file_location+"/trajectories/"+tr,delimiter="\t")
    res["trajectories"] = trajectories
    
    #Read and store capacity data
    capacity_data = pd.read_csv(file_location+"/capacities.txt",delimiter="\t")
    
    node_ids = capacity_data["step"]
    res["node_ids"] = node_ids
    
    node_names = {}
    capacities = {}
    node_cost_per_time = {}
    for i in range(len(node_ids)):
        node_names[node_ids[i]] = capacity_data["resource_name"][i]
        capacities[node_ids[i]] = capacity_data["number_of_resources"][i]
        node_cost_per_time[node_ids[i]] = capacity_data["resource_cost_per_time"][i]
    res["node_names"] = node_names
    res["capacities"] = capacities
    res["node_cost_per_time"] = node_cost_per_time
    
    
    #Read and store arrival data
    arrival_data = pd.read_csv(file_location+"/arrivals.txt",delimiter="\t")
    
    res["arrival_distribution"] = arrival_data["arrival_distribution"][0]
    res["arrival_distr_parameter"] = arrival_data["arrival_rate"][0]
    
    traj_prop = {}
    for traj in trajectories:
        traj_prop[traj] = arrival_data[traj][0]
    res["traj_prop"] = traj_prop

    return res




