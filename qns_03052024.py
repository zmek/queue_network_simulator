#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#Input: for each patient type, a matrix of probabilities of moving from one node
#to another node. Might use an initial dummy node for better generalisation.
#For each patient: 1) read the matrix of probabilities (and determine the first
#node if not using a dummy initial node); 2) select current node; 3) wait for
#patient to be served; 4) determine next node; 5) if next node is an "exit" node
#then stop, otherwise go to step 2.





import simpy
import numpy
import random
import pandas
import math
import matplotlib.pyplot
from matplotlib.backends.backend_pdf import PdfPages



class Params:
    
    analysisID = "test" #name of the folder containing input and output folders of the current analysis
    
    analysis_location = "data" #relative path from this Python script to the folder containing [analysisID]'s folder
    
    number_of_runs = 0 #ideally 100 or above
    warm_up_duration = 0 #days - warm-up period, to be set to 0 if we allow a non-empty system as starting state
    sim_duration = 0 #days - duration of the actual simulation after the warm-up period
    
    first_time_point = 0 #i-th week after warm-up period (between 0 and sim_duration) for which we want to have a snapshot of waiting times
    second_time_point = 0 #j-th day after warm-up period (between 0 and sim_duration) for which we want to have a snapshot of waiting times
    sojourn_time_threshold = 0 #target maximum waiting time (in days)
    
    use_initial_patients = False
    initial_state_source = None
    generate_initial_patients = False
    
    plot_overall_sojourn_time = False
    plot_comparison = False
    plot_each_node = False
    
    #Example of list of initial patients - technical parameters
    print_generated_patient_list = False
    patient_list_run = 1 #simulation run from which an example of initial patient list is taken
    initial_patient_list = pandas.DataFrame()
    initial_patient_list["patID"] = [] 
    initial_patient_list["pat_stream"] = []
    initial_patient_list["truth"] = []
    initial_patient_list["current_node"] = []
    initial_patient_list["time_waited_at_node"] = []
    initial_patient_list["additional_time_in_system"] = []
    
    
    


# Class representing the nodes of the network
# Nodes could either represent discrete servers or pools of resources
# Capacity is set up according to the type of resource, whereas service rates
# and replenishment levels/frequencies are defined based on patient types
class Node:
    
    def __init__(self, env, node_name, node_type, resource_type, capacity, service_distributions, service_parameters):
        self.id = node_name
        self.type = node_type  #allowed options: {"entry","internal","exit"}
        
        self.resource = None
        self.replenish_freq = float("inf")
        self.resource_type = resource_type  #allowed options: {"gradual","batch"}
        
        if self.type == "internal":
            if self.resource_type == "gradual":
                self.resource = simpy.PriorityResource(env,capacity=capacity[0])
            else: #i.e. self.resource_type = "batch"
                self.resource = simpy.PriorityResource(env,capacity=capacity[0])
                self.replenish_freq = capacity[1]
                
        if self.type == "exit":
            self.resource = simpy.Container(env,capacity=100,init=100)
            
        self.service_distributions = service_distributions
        self.service_parameters = service_parameters
        
        

# Class representing the patient types
class PatientType:
    
    def __init__(self, pat_type, network_structure, test_based_branching):
        self.id = pat_type
        self.network_structure = network_structure
        self.test_based_branching = test_based_branching
        
    #function that picks the next node for the patient given the current node
    #based on the transition probabilities and on possible test sensitivity/specificity
    def getNextNode(self,current_node,stream,truth):
        
        options = list(self.network_structure[current_node].keys())
        vals = list(self.network_structure[current_node].values())
        
        res = None
        
        if self.test_based_branching[current_node]["positive_outcome"] == None or self.test_based_branching[current_node]["negative_outcome"] == None:
        
            sel_vals = [x for x in vals if x>0]
            if sum(sel_vals) != 1:
                print("ERROR: elements of patient flow matrix must sum up to 1 for branching that does not involve sensitivity/specificity features")
                exit()
                
            sel_options = [options[i] for i in range(len(options)) if vals[i]>0]
            
            res = random.choices(sel_options, weights = sel_vals)[0]
        
        else:
            
            sens_spec_options = [ self.test_based_branching[current_node]["positive_outcome"], self.test_based_branching[current_node]["negative_outcome"] ]
            sens = self.test_based_branching[current_node]["sensitivity"]
            spec = self.test_based_branching[current_node]["specificity"]
        
            sel_vals = [x for x in vals if ( x>0 and x<1 ) ]
            if sum(sel_vals) >= 1:
                print("ERROR: elements of patient flow matrix must sum up to strictly less than 1 (excluding sensitivity/specificity nodes) for branching that involves sensitivity/specificity features")
                exit()
            
            sel_options = [options[i] for i in range(len(options)) if ( vals[i]>0 and vals[i]<1 ) ]
            
            sel_options.append("sens_spec")
            sel_vals.append( 1 - sum(sel_vals) )
            
            res = random.choices(sel_options, weights = sel_vals)[0]
            
            if res == "sens_spec":
                
                if truth == "positive":
                
                    res = random.choices(sens_spec_options, weights = [ sens, 1 - sens ] )[0]
                    
                else:
                    
                    res = random.choices(sens_spec_options, weights = [ 1 - spec, spec ] )[0]
        
        return res    


# Class representing the patients
class Patient:
    
    def __init__(self, pat_id, pat_type, arrival_time, category, time_joined_current_queue):
        self.id = pat_id
        self.type = str(pat_type).split('_')[0]
        self.current_node_id = ""
        
        self.arrival_time = arrival_time
        
        self.category = category
        
        self.stream = str(pat_type).split('_')[0]
        self.truth = str(pat_type).split('_')[1]
        
        self.time_joined_current_queue = time_joined_current_queue
        


# Class representing the pathway model
class PathwayModel:
    
    def __init__(self, run_number, data, warm_up_duration, sim_duration):
    
        self.env = simpy.Environment()
        self.run_number = run_number
        self.warm_up_duration = warm_up_duration
        self.sim_duration = sim_duration
        self.initial_conditions_offset = 0
        
        self.patient_counter = 0
        
        self.patient_types = {} #dictionary containing patient types (key is patient type's ID, value is PatientType object)
        for key in data["patient_types"]:
            self.patient_types[key] = PatientType(key, data["network_structure"][key], data["test_based_branching"][key])
        
        self.entry_node_ids = data["entry_nodes"]
        self.internal_node_ids = data["internal_nodes"]
        self.exit_node_ids = data["exit_nodes"]
        
        self.nodes = {} #dictionary containing network nodes (key is node's ID, value is Node object)
        for key in self.entry_node_ids:
            self.nodes[key] = Node(self.env, key, "entry", "", float("inf"), None, None)
        for key in self.internal_node_ids:
            self.nodes[key] = Node(self.env, key, "internal", data["node_types"][key], data["capacity_parameters"][key], data["service_distributions"][key], data["service_parameters"][key])
        for key in self.exit_node_ids:
            self.nodes[key] = Node(self.env, key, "exit", "", float("inf"), None, None)
        
        
        self.arrival_distributions = data["arrival_distributions"] 
        self.arrival_parameters = data["arrival_parameters"]
        
        if(Params.use_initial_patients):
            if Params.generate_initial_patients:
            
                self.initial_number_distributions = data["initial_number_distributions"]
                self.initial_number_parameters = data["initial_number_parameters"]
        
                self.initial_node_wait_distributions = data["initial_node_wait_distributions"]
                self.initial_node_wait_parameters = data["initial_node_wait_parameters"]
        
                self.initial_additional_overall_sojourn_time_distributions = data["initial_additional_overall_sojourn_time_distributions"]
                self.initial_additional_overall_sojourn_time_parameters = data["initial_additional_overall_sojourn_time_parameters"]
            
                self.proportion_positive_initial = data["proportion_positive_initial"]
        
            else:
            
                self.initial_number_distributions = None
                self.initial_number_parameters = None
        
                self.initial_node_wait_distributions = None
                self.initial_node_wait_parameters = None
        
                self.initial_additional_overall_sojourn_time_distributions = None
                self.initial_additional_overall_sojourn_time_parameters = None
            
                self.proportion_positive_initial = None
            
                self.initial_patient_list = data["initial_patients"]
        
        
        self.proportion_positive = data["proportion_positive"]
        
        self.dna_rates = data["dna_rates"]
        
            
    #function generating arrivals for the specified patient type and entry nodes
    #arrival_process_parameters are a dictionary as they might differ depending on the underlying probability distribution
    def generate_arrivals(self, pat_type_id, entry_node_id, arrival_process_distribution, arrival_process_parameters, prop_positive):
        
        if arrival_process_parameters[0] > 0:  #if rate is 0, then this function does not do anything
            
            while True:
                
                #determine the next arrival time and freeze until then
                if arrival_process_distribution == "poisson":
                    time_to_next_arrival = random.expovariate(arrival_process_parameters[0])
                else: #assuming the inter-arrival times are deterministic
                    time_to_next_arrival = 1 / arrival_process_parameters[0]
                
                yield self.env.timeout(time_to_next_arrival)
                
                #create the new patient
                if self.env.now >= self.initial_conditions_offset:
                    
                    self.patient_counter += 1
                    
                    truth_state = random.choices( [ "positive", "negative" ], weights = [ prop_positive, 1 - prop_positive ] )[0]
                    #p = Patient(self.patient_counter, pat_type_id, numpy.ceil(self.env.now))
                    p = Patient(self.patient_counter, pat_type_id + "_" + truth_state, self.env.now, "new_arrival", self.env.now)
                    
                    #start patient p's journey through the network
                    self.env.process(self.patient_journey(p,entry_node_id))
                
    
                
    
    def generate_initial_waiting_times(self, pat_type_id, current_node_id, initial_number_distribution, initial_number_parameters, initial_node_wait_distribution, initial_node_wait_parameters):
    
        if Params.generate_initial_patients:
            
            initial_patients = 0
            wait_times = []
        
            if initial_number_distribution == "poisson":
                initial_patients = int(numpy.random.poisson(initial_number_parameters[0]))
            else:
                initial_patients = int(initial_number_parameters[0])
    
    
            if initial_patients > 0:
                #generate node waiting times and sort them from the biggest to the smallest
                if initial_node_wait_distribution == "lognormal":
                    obs_mean = initial_node_wait_parameters[0]
                    obs_sd = initial_node_wait_parameters[1]
                    mu = math.log(obs_mean/math.sqrt(1+obs_sd**2/obs_mean**2))
                    sigma = math.sqrt(math.log(1+obs_sd**2/obs_mean**2))
                    #wait_times = [ numpy.ceil(numpy.random.lognormal( mu, sigma )) for i in range(initial_patients) ]
                    wait_times = [ numpy.random.lognormal( mu, sigma ) for i in range(initial_patients) ]
                else: #deterministic
                    #wait_times = [ numpy.ceil(initial_node_wait_parameters[0]) for i in range(initial_patients) ]
                    wait_times = [ initial_node_wait_parameters[0] for i in range(initial_patients) ]
               
                wait_times.sort(reverse=True)
                
        else:
            
            sel_patients = self.initial_patient_list[self.initial_patient_list["current_node"]==current_node_id]
            sel_patients = sel_patients.sort_values(by="time_waited_at_node",ascending=False)
            
            wait_times = list(sel_patients["time_waited_at_node"])
            
                   
        return wait_times
            
            
    
    def generate_initial_patients(self, pat_type_id, current_node_id, initial_wait_times, initial_additional_overall_sojourn_time_distribution, initial_additional_overall_sojourn_time_parameters, prop_positive):
        
        if len(initial_wait_times) > 0:
        
            if Params.generate_initial_patients:
            
                #generate additional overall sojourn time (to be added to the waiting time at the current node)
                if initial_additional_overall_sojourn_time_distribution == "lognormal":
                    obs_mean = initial_additional_overall_sojourn_time_parameters[0]
                    obs_sd = initial_additional_overall_sojourn_time_parameters[1]
                    mu = math.log(obs_mean/math.sqrt(1+obs_sd**2/obs_mean**2))
                    sigma = math.sqrt(math.log(1+obs_sd**2/obs_mean**2))
                    #additional_sojourn_times = [ numpy.ceil(numpy.random.lognormal( mu, sigma )) for i in range(initial_patients) ]
                    additional_sojourn_times = [ numpy.random.lognormal( mu, sigma ) for i in range(len(initial_wait_times)) ]
                else: #deterministic
                    #additional_sojourn_times = [ numpy.ceil(initial_additional_overall_sojourn_time_parameters[0]) for i in range(initial_patients) ]
                    additional_sojourn_times = [ initial_additional_overall_sojourn_time_parameters[0] for i in range(len(initial_wait_times)) ]
            
            else:
                
                sel_patients = self.initial_patient_list[self.initial_patient_list["current_node"]==current_node_id]
                sel_patients = sel_patients.sort_values(by="time_waited_at_node",ascending=False)
            
                additional_sojourn_times = list(sel_patients["additional_time_in_system"])
                
            
            #block resources for the initial offset time
            for i in range(int(self.nodes[current_node_id].resource.capacity)):
            
                yield self.env.timeout(0)
            
                self.env.process(self.offset_elapse(current_node_id,self.initial_conditions_offset))
            
            #determine sequence of inter-arrival times of patients of the current type at the current node
            inter_arr_times = [ self.initial_conditions_offset - initial_wait_times[0] ]
            for iter in range(len(initial_wait_times)-1):
                inter_arr_times.append( initial_wait_times[iter] - initial_wait_times[iter+1] )
            
            #generate patients and assign them a relative arrival time to the system
            for iter in range(len(initial_wait_times)):
                
                time_to_next_arrival = inter_arr_times[iter]
                
                yield self.env.timeout(time_to_next_arrival)
                
                if Params.generate_initial_patients:
                    truth_state = random.choices( [ "positive", "negative" ], weights = [ prop_positive, 1 - prop_positive ] )[0]
                else:
                    truth_state = list(sel_patients["truth"])[iter]
                
                #create the new patient ID
                self.patient_counter += 1
                #p = Patient(self.patient_counter, pat_type_id, numpy.ceil(self.env.now) - additional_sojourn_times[iter] - wait_times[iter] )
                p = Patient(self.patient_counter, pat_type_id + "_" + truth_state, self.env.now - additional_sojourn_times[iter], "initial", self.env.now )
                
                if self.run_number == Params.patient_list_run:
                    new_pat_list_record = pandas.DataFrame({
                        "patID":[p.id],
                        "pat_stream":[p.stream],
                        "truth":[p.truth],
                        "current_node":[current_node_id],
                        "time_waited_at_node":[initial_wait_times[iter]],
                        "additional_time_in_system":[additional_sojourn_times[iter]]
                        })
                    #Params.initial_patient_list = Params.initial_patient_list.append(new_pat_list_record)
                    Params.initial_patient_list = pandas.concat([Params.initial_patient_list, new_pat_list_record])
                
                #start patient p's journey through the network
                self.env.process(self.patient_journey(p,current_node_id))
    
                
    
    def offset_elapse(self, node_id, offset_time):
        
        with self.nodes[node_id].resource.request(priority=-1000) as req:
        
            yield req
            
            yield self.env.timeout(offset_time)
            
    
    def patient_journey(self, patient, entry_node_id):
        
        #Retrieve patient type
        pt = patient.type
        curr_patient_type = self.patient_types[pt]
        
        #Initialise current node's ID
        curr_node_id = entry_node_id
        patient.current_node_id = curr_node_id
        
        while True:
            
            if self.nodes[curr_node_id].type == "exit":
                yield self.nodes[curr_node_id].resource.get(1)
                yield self.env.timeout(0)
                yield self.nodes[curr_node_id].resource.put(1)
                
                #record sojourn time and update sojourn time output
                new_soj_time_record = pandas.DataFrame({
                    "patID":[patient.id],
                    "pat_truth":[patient.truth],
                    "category":[patient.category],
                    "pat_stream":[patient.stream],
                    "arrival_time":[patient.arrival_time],
                    "departure_time":[self.env.now],
                    "sojourn_time":[self.env.now - patient.arrival_time],
                    "outcome":[self.nodes[curr_node_id].id]
                    })
                #self.output.sojourn_time = self.output.sojourn_time.append(new_soj_time_record)
                self.output.sojourn_time = pandas.concat([self.output.sojourn_time, new_soj_time_record])
                break
                
                
            if self.nodes[curr_node_id].type == "entry":
                
                #identify and retrieve the next node in the pathway
                curr_node_id = curr_patient_type.getNextNode(curr_node_id,patient.stream,patient.truth)
                patient.time_joined_current_queue = self.env.now
                
                #Update patient's current node
                patient.current_node_id = curr_node_id
            
            #If current node is not an exit node, then request resource
            if self.nodes[curr_node_id].type == "internal":
                
                #keep track of arrival time at the current node
                #queue_start = numpy.ceil(self.env.now)
                queue_start = patient.time_joined_current_queue
                
                dna_priority = 0
                
                while True:
                
                    #Check whether the node processes demand gradually or in batches and manage resources accordingly
                    if self.nodes[curr_node_id].resource_type == "gradual":
                        
                        with self.nodes[curr_node_id].resource.request(priority=dna_priority) as req: 
                            #Freeze the function until the request can be met
                            yield req
                            
                            # Determine the time the patient will spend at the current node 
                            if self.nodes[curr_node_id].service_distributions == "exponential":
                                serv_time = random.expovariate(1 / self.nodes[curr_node_id].service_parameters[0])
                            elif self.nodes[curr_node_id].service_distributions == "lognormal":
                                obs_mean = self.nodes[curr_node_id].service_parameters[0]
                                obs_sd = self.nodes[curr_node_id].service_parameters[1]
                                mu = math.log(obs_mean/math.sqrt(1+obs_sd**2/obs_mean**2))
                                sigma = math.sqrt(math.log(1+obs_sd**2/obs_mean**2))
                                serv_time = numpy.random.lognormal( mu, sigma )
                            else: #assuming the service times are deterministic
                                serv_time = self.nodes[curr_node_id].service_parameters[0]
                        
                            #update output about resource utilisation
                            ind_start = int(numpy.floor(self.env.now))
                            ind_stop = min(int(numpy.floor(self.env.now + serv_time)),self.warm_up_duration+self.sim_duration+2)
                        
                            for t in range(ind_start,ind_stop):
                                
                                old_value = self.output.resource_utilisation[patient.current_node_id][patient.stream][t]
                                new_value = old_value + 1
                                self.output.resource_utilisation[patient.current_node_id].at[t,patient.stream] = new_value
                            
                                res_util_den = self.output.resource_utilisation[patient.current_node_id]["capacity"][t]
                                res_util_num = 0
                                for p in list(self.patient_types.keys()):
                                    res_util_num = res_util_num + self.output.resource_utilisation[patient.current_node_id][str(p).split('_')[0]][t]
                                    self.output.resource_utilisation[patient.current_node_id].at[t,"resource_utilisation"] = res_util_num / res_util_den
    
                            # Freeze this function until that time has elapsed
                            yield self.env.timeout(serv_time)
                            
                            #keep track of time when patient moves to the next node (note, here we define waiting time as waiting time in the queue plus service time)
                            queue_finish = self.env.now
                        
                    else: #assuming here that demand is processed in batches
                    
                        with self.nodes[curr_node_id].resource.request(priority=dna_priority) as req: 
                            # Freeze the function until the request can be met
                            yield req
                        
                            # Keep slots allocated for the appropriate time
                            curr_time = self.env.now
                            time_since_last_replenishment = curr_time % self.nodes[curr_node_id].replenish_freq
                            serv_time = self.nodes[curr_node_id].replenish_freq - time_since_last_replenishment
                    
                    
                            #update output about resource utilisation
                            ind_start = int(numpy.floor(self.env.now))
                            ind_stop = min(int(numpy.floor(self.env.now + serv_time)),self.warm_up_duration+self.sim_duration+2)
                        
                            for t in range(ind_start,ind_stop):
                                
                                old_value = self.output.resource_utilisation[patient.current_node_id][patient.stream][t]
                                new_value = old_value + 1
                                self.output.resource_utilisation[patient.current_node_id].at[t,patient.stream] = new_value
                            
                                res_util_den = self.output.resource_utilisation[patient.current_node_id]["capacity"][t]
                                res_util_num = 0
                                for p in list(self.patient_types.keys()):
                                    res_util_num = res_util_num + self.output.resource_utilisation[patient.current_node_id][str(p).split('_')[0]][t]
                                    self.output.resource_utilisation[patient.current_node_id].at[t,"resource_utilisation"] = res_util_num / res_util_den
    
                            # Freeze this function until that time has elapsed
                            yield self.env.timeout(serv_time)
                        
                            #keep track of time when patient moves to the next node (note, here we define waiting time as waiting time in the queue plus service time)
                            queue_finish = self.env.now
                            
                    
                    #Check if the patient actually showed up - if so exit the loop, otherwise repeat with higher priority
                    dna_prob = self.dna_rates[curr_node_id]
                    showed_up = random.choices(["no","yes"],weights=[ dna_prob, 1 - dna_prob ])[0]
                    
                    if showed_up == "yes":
                        break
                    else:
                        dna_priority = -100
                
                #record time at which the patient moves to the next node
                move_on_time = self.env.now
                
                #update current number of patients processed
                ind = int(numpy.floor(self.env.now))
                current_value = self.output.cumulative_patients[patient.current_node_id][patient.stream]["pat_gone"][ind]
                new_value = current_value + 1
                self.output.cumulative_patients[patient.current_node_id][patient.stream].at[ind,"pat_gone"] = new_value
                     
                #record waiting time and update waiting time output
                new_wait_time_record = pandas.DataFrame({
                    "patID":[patient.id],
                    "pat_truth":[patient.truth],
                    "category":[patient.category],
                    "pat_stream":[patient.stream],
                    "node":self.nodes[curr_node_id].id,
                    "time_joined_queue":[queue_start],
                    "time_started_service":[queue_finish],
                    "waiting_time":[queue_finish - queue_start],
                    "time_moved_on":[move_on_time]
                    })
                #self.output.waiting_time = self.output.waiting_time.append(new_wait_time_record)
                self.output.waiting_time = pandas.concat([self.output.waiting_time, new_wait_time_record])
                
                
                #identify and retrieve the next node in the pathway
                curr_node_id = curr_patient_type.getNextNode(curr_node_id,patient.stream,patient.truth)
                patient.time_joined_current_queue = self.env.now
                
                #Update patient's current node
                patient.current_node_id = curr_node_id
                
             
        
    def run(self, run_number):
                    
        if Params.use_initial_patients:
        
            # Generate initial set of patients and set up offset period
            initial_wait_times = {}
            for node_id in self.internal_node_ids:
                initial_wait_times[node_id] = {}
                for pt_id in list(my_model.patient_types.keys()):
                    if Params.generate_initial_patients:
                        initial_wait_times[node_id][pt_id] = self.generate_initial_waiting_times(pt_id, node_id, self.initial_number_distributions[node_id][pt_id],self.initial_number_parameters[node_id][pt_id],self.initial_node_wait_distributions[node_id][pt_id], self.initial_node_wait_parameters[node_id][pt_id])    
                    else:
                        initial_wait_times[node_id][pt_id] = self.generate_initial_waiting_times(pt_id, node_id, None, None, None, None)
                    self.initial_conditions_offset = max( self.initial_conditions_offset, max(initial_wait_times[node_id][pt_id]) )
                    
        self.output = Output(run_number, list(self.patient_types.keys()), self.internal_node_ids, int(numpy.ceil(self.warm_up_duration + self.sim_duration + self.initial_conditions_offset)))
        for key in list(self.output.resource_utilisation.keys()):
            self.output.resource_utilisation[key]["capacity"] = [self.nodes[key].resource.capacity for i in range(int(numpy.ceil(self.warm_up_duration + self.sim_duration + self.initial_conditions_offset + 2)))]
        
        
        if Params.use_initial_patients:
        
            for node_id in self.internal_node_ids:
                for pt_id in list(my_model.patient_types.keys()):
                    if Params.generate_initial_patients:
                        self.env.process(self.generate_initial_patients(pt_id, node_id, initial_wait_times[node_id][pt_id], self.initial_additional_overall_sojourn_time_distributions[node_id][pt_id], self.initial_additional_overall_sojourn_time_parameters[node_id][pt_id],self.proportion_positive_initial[node_id][pt_id]))
                    else:
                        self.env.process(self.generate_initial_patients(pt_id, node_id, initial_wait_times[node_id][pt_id], None, None, None))
                            
        
        # Start entity generators
        for node_id in self.arrival_distributions:
            for pt_id in self.arrival_distributions[node_id]:
                self.env.process(self.generate_arrivals(pt_id,node_id,self.arrival_distributions[node_id][pt_id],self.arrival_parameters[node_id][pt_id],self.proportion_positive[node_id][pt_id]))
        
       
        # Run simulation
        self.env.run(until=(self.warm_up_duration + self.sim_duration + self.initial_conditions_offset))
        
        #update cumulative number of patients processed
        for node_key in list(self.output.cumulative_patients.keys()):
            for pat_type_key in list(self.output.cumulative_patients[node_key].keys()):
                for t in self.output.cumulative_patients[node_key][pat_type_key]["time"][1:]:
                    self.output.cumulative_patients[node_key][pat_type_key]["cumul_pat_gone"][t] = sum(self.output.cumulative_patients[node_key][pat_type_key]["pat_gone"][:(t+1)])
        
        
        #remove records regarding the initial offset period
        t = self.initial_conditions_offset
        
        self.output.sojourn_time["arrival_time"] = [ ( x - t ) for x in self.output.sojourn_time["arrival_time"] ]
        self.output.sojourn_time["departure_time"] = [ ( x - t ) for x in self.output.sojourn_time["departure_time"] ]
                
        self.output.waiting_time["time_joined_queue"] = [ ( x - t ) for x in self.output.waiting_time["time_joined_queue"] ]
        self.output.waiting_time["time_started_service"] = [ ( x - t ) for x in self.output.waiting_time["time_started_service"] ]
        self.output.waiting_time["time_moved_on"] = [ ( x - t ) for x in self.output.waiting_time["time_moved_on"] ]
                              
        for node_key in list(self.output.cumulative_patients.keys()):
            
            self.output.resource_utilisation[node_key] = self.output.resource_utilisation[node_key][self.output.resource_utilisation[node_key]["time"]>t]
            self.output.resource_utilisation[node_key]["time"] = [ x for x in range(len(self.output.resource_utilisation[node_key])) ]
            self.output.resource_utilisation[node_key].index = self.output.resource_utilisation[node_key]["time"]
        
            for pat_type_key in list(self.output.cumulative_patients[node_key].keys()):
                           
                self.output.cumulative_patients[node_key][pat_type_key] = self.output.cumulative_patients[node_key][pat_type_key][self.output.cumulative_patients[node_key][pat_type_key]["time"]>t]
                self.output.cumulative_patients[node_key][pat_type_key]["time"] = [ x for x in range(len(self.output.cumulative_patients[node_key][pat_type_key])) ]
                self.output.cumulative_patients[node_key][pat_type_key].index = self.output.cumulative_patients[node_key][pat_type_key]["time"]



# Class representing and manipulating the output of the simulation runs
class Output:
    
    def __init__(self, run_number, patient_type_ids, node_ids, time_frames):
        
        self.run_number = run_number
        
        patient_streams = list(set([str(x).split('_')[0] for x in patient_type_ids]))
        
        self.cumulative_patients = {}
        for node_key in node_ids:
            self.cumulative_patients[node_key] = {}
            for pat_type_key in patient_streams:
                self.cumulative_patients[node_key][pat_type_key] = pandas.DataFrame()
                self.cumulative_patients[node_key][pat_type_key]["time"] = [i for i in range(time_frames+2)]
                self.cumulative_patients[node_key][pat_type_key]["pat_gone"] = [0 for i in range(time_frames+2)]
                self.cumulative_patients[node_key][pat_type_key]["cumul_pat_gone"] = [0 for i in range(time_frames+2)]
        
        self.sojourn_time = pandas.DataFrame()
        self.sojourn_time["patID"] = []
        self.sojourn_time["pat_truth"] = []
        self.sojourn_time["category"] = []
        self.sojourn_time["pat_stream"] = []
        self.sojourn_time["arrival_time"] = []
        self.sojourn_time["departure_time"] = []
        self.sojourn_time["sojourn_time"] = []
        self.sojourn_time["outcome"] = []
        
        self.waiting_time = pandas.DataFrame()
        self.waiting_time["patID"] = []
        self.waiting_time["pat_truth"] = []
        self.waiting_time["category"] = []
        self.waiting_time["pat_stream"] = []
        self.waiting_time["node"] = []
        self.waiting_time["time_joined_queue"] = []
        self.waiting_time["time_started_service"] = []
        self.waiting_time["waiting_time"] = []
        self.waiting_time["time_moved_on"] = []
        
        self.resource_utilisation = {}
        for node_key in node_ids:
            self.resource_utilisation[node_key] = pandas.DataFrame()
            self.resource_utilisation[node_key]["time"] = [i for i in range(time_frames+2)]
            for pat_type_key in patient_streams:
                self.resource_utilisation[node_key][pat_type_key] = [0 for i in range(time_frames+2)]
            self.resource_utilisation[node_key]["capacity"] = [0 for i in range(time_frames+2)]
            self.resource_utilisation[node_key]["resource_utilisation"] = [0.0 for i in range(time_frames+2)]

    
    def print(self, file_name, start, finish):
        
        pp = PdfPages(file_name)
        
        for node_key in list(self.cumulative_patients.keys()):
            for pat_type_key in list(self.cumulative_patients[node_key].keys()):
                plot = matplotlib.pyplot.figure()
                matplotlib.pyplot.plot(self.cumulative_patients[node_key][pat_type_key]["cumul_pat_gone"][start:finish])
                matplotlib.pyplot.title("Cumulative patients processed - node: " + node_key + " -  patient type: " + pat_type_key)
                matplotlib.pyplot.xlabel("time")
                matplotlib.pyplot.ylabel("Number of patients")
                pp.savefig(plot)
                
        
        hist_soj = matplotlib.pyplot.figure()
        matplotlib.pyplot.hist(self.sojourn_time["sojourn_time"][start:finish])
        matplotlib.pyplot.title("Sojourn time distribution")
        matplotlib.pyplot.xlabel("sojourn time")
        matplotlib.pyplot.ylabel("frequency")
        pp.savefig(hist_soj)
        
        hist_wait = matplotlib.pyplot.figure()
        matplotlib.pyplot.hist(self.waiting_time["waiting_time"][start:finish])
        matplotlib.pyplot.title("Waiting time distribution")
        matplotlib.pyplot.xlabel("waiting time")
        matplotlib.pyplot.ylabel("frequency")
        pp.savefig(hist_wait)
        
        for node_key in list(self.resource_utilisation.keys()):
            plot = matplotlib.pyplot.figure()
            matplotlib.pyplot.plot(self.resource_utilisation[node_key]["resource_utilisation"][start:finish])
            matplotlib.pyplot.title("Resource utilisation - node: " + node_key)
            matplotlib.pyplot.xlabel("time")
            matplotlib.pyplot.ylabel("proportion")
            matplotlib.pyplot.ylim([0, 1.05])
            pp.savefig(plot)
        
        pp.close()
        



def time_to_diagnosis_plot(data, finish, threshold):
    
    start = 0
    
    all_results = data[1].sojourn_time
    if len(data)>1:
        for i in range(2,len(data)):
            #all_results = all_results.append(data[i].sojourn_time)
            all_results = pandas.concat([all_results, data[i].sojourn_time])
    
    sel_records = all_results
    sel_records = sel_records[sel_records["arrival_time"] >= start + Params.warm_up_duration]
    sel_records = sel_records[sel_records["departure_time"] < finish + Params.warm_up_duration]
    
    soj_times = list(sel_records["sojourn_time"])
    prop_in_target = 0
    if len(soj_times) > 0:
        prop_in_target = int( round( len([x for x in soj_times if x <= threshold ]) / len(soj_times) * 100, 0) )
    
    sel_records_positive = sel_records[sel_records["outcome"] == "CancerDiagnosis"]
    soj_times_positive = list(sel_records_positive["sojourn_time"])
    prop_in_target_positive = 0
    if len(soj_times_positive) > 0:
        prop_in_target_positive = int( round( len([x for x in soj_times_positive if x <= threshold ]) / len(soj_times_positive) * 100, 0) )
    
    sel_records_negative = sel_records[sel_records["outcome"] != "CancerDiagnosis"]
    soj_times_negative = list(sel_records_negative["sojourn_time"])
    prop_in_target_negative = 0
    if len(soj_times_negative) > 0:
        prop_in_target_negative = int( round( len([x for x in soj_times_negative if x <= threshold ]) / len(soj_times_negative) * 100, 0) )
    
    if len(soj_times) > 0:
        x_max = max( threshold, max( soj_times ) )
        
        counts, bins = numpy.histogram(soj_times,numpy.arange(start=0,stop=math.ceil(x_max)+0.1,step=1),density=True)
        y_max = max(counts)
    else:
        x_max = threshold
        y_max = 1
        
    x_max = 50
    y_max = 0.5
    
    fig, ((ax1, ax2), (ax3, ax4)) = matplotlib.pyplot.subplots(2, 2, figsize=(10,10))
    fig.suptitle("Overall time on pathway", fontsize = 24)
    fig.tight_layout(pad=4)
    
    ax1.hist(soj_times,density=True,bins=numpy.arange(start=0,stop=math.ceil(x_max)+0.1,step=1))
    ax1.axvline(x = threshold, color = 'r')
    ax1.set_title("All patients: " + str(prop_in_target) + "% left within " + str(threshold) + " days")
    ax1.set_xlabel("Days on pathway")
    ax1.set_ylabel("Proportion of patients")
    ax1.set_xlim(0,x_max*1.1)
    ax1.set_ylim(0,y_max*1.1)
    
    #y_vals = ax1.get_yticks()
    #ax1.set_yticklabels([int(round(y/Params.number_of_runs,0)) for y in y_vals])
    
    tp = len(sel_records_positive[sel_records_positive["pat_truth"]=="positive"])
    fp = len(sel_records_positive[sel_records_positive["pat_truth"]=="negative"])
    fn = len(sel_records_negative[sel_records_negative["pat_truth"]=="positive"])
    tn = len(sel_records_negative[sel_records_negative["pat_truth"]=="negative"])
    
    df = pandas.DataFrame({
        "":["No-cancer patients","Cancer patients"],
        "No-cancer diagnoses":[tn,fn],
        "Cancer diagnoses":[fp,tp]
        })
    
    
    ax2.axis("off")
    ax2.table(cellText=df.values,
              colLabels=["","No-cancer diagnoses","Cancer diagnoses"],
              loc="center")
    
    ax3.hist(soj_times_positive,density=True,bins=numpy.arange(start=0,stop=math.ceil(x_max)+0.1,step=1),color="orange")
    ax3.axvline(x = threshold, color = 'r')
    ax3.set_title("Cancer diagnosis: " + str(prop_in_target_positive) + "% left within " + str(threshold) + " days")
    ax3.set_xlabel("Days on pathway")
    ax3.set_ylabel("Proportion of patients")
    ax3.set_xlim(0,x_max*1.1)
    ax3.set_ylim(0,y_max*1.1)
    
    ax4.hist(soj_times_negative,density=True,bins=numpy.arange(start=0,stop=math.ceil(x_max)+0.1,step=1),color="green")
    ax4.axvline(x = threshold, color = 'r')
    ax4.set_title("No-cancer diagnosis: " + str(prop_in_target_negative) + "% left within " + str(threshold) + " days")
    ax4.set_xlabel("Days on pathway")
    ax4.set_ylabel("Proportion of patients")
    ax4.set_xlim(0,x_max*1.1)
    ax4.set_ylim(0,y_max*1.1)

    return fig
    

def compare_time_points_overall(data, start, finish, threshold):
    
    all_results = data[1].sojourn_time
    if len(data)>1:
        for i in range(2,len(data)):
            #all_results = all_results.append(data[i].sojourn_time)
            all_results = pandas.concat([all_results, data[i].sojourn_time])
              
    sel_records = all_results
    sel_records = sel_records[sel_records["arrival_time"] < start + Params.warm_up_duration]
    sel_records = sel_records[sel_records["departure_time"] >= start + Params.warm_up_duration]
    
    soj_times_1 = list(sel_records["sojourn_time"])
     
    sel_records = all_results
    sel_records = sel_records[sel_records["arrival_time"] < finish + Params.warm_up_duration]
    sel_records = sel_records[sel_records["departure_time"] >= finish + Params.warm_up_duration]
    
    soj_times_2 = list(sel_records["sojourn_time"])
     
    if len(soj_times_1 + soj_times_2) > 0:
        x_max = max( threshold, max( soj_times_1 + soj_times_2 ) )
        
        counts, bins = numpy.histogram(soj_times_1+soj_times_2,numpy.arange(start=0,stop=math.ceil(x_max)+0.1,step=1),density=True)
        y_max = max(counts)
    else:
        x_max = threshold
        y_max = 1
        
    x_max = 50
    y_max = 0.5
    
    fig, ((ax1, ax2), (ax3, ax4)) = matplotlib.pyplot.subplots(2, 2, figsize=(10,10))
    fig.suptitle("Overall time on pathway - comparison", fontsize = 24)
    fig.tight_layout(pad=4)
    
    ax1.hist(soj_times_1,density=True,bins=numpy.arange(start=0,stop=math.ceil(x_max)+0.1,step=1))
    ax1.axvline(x = threshold, color = 'r')
    ax1.set_title("Day " + str(start))
    ax1.set_xlabel("Current total waiting time (days)")
    ax1.set_ylabel("Proportion of patients")
    ax1.set_xlim(0,x_max*1.1)
    ax1.set_ylim(0,y_max*1.1)
    
    #y_vals = ax1.get_yticks()
    #ax1.set_yticklabels([int(round(y/Params.number_of_runs,0)) for y in y_vals])
    
    ax2.hist(soj_times_2,density=True,bins=numpy.arange(start=0,stop=math.ceil(x_max)+0.1,step=1))
    ax2.axvline(x = threshold, color = 'r')
    ax2.set_title("Day " + str(finish))
    ax2.set_xlabel("Current total waiting time (days)")
    ax2.set_ylabel("Proportion of patients")
    ax2.set_xlim(0,x_max*1.1)
    ax2.set_ylim(0,y_max*1.1)
    
    #y_vals = ax2.get_yticks()
    #ax2.set_yticklabels([int(round(y/Params.number_of_runs,0)) for y in y_vals])
    
    ax3.axis("off")
    ax4.axis("off")
       
    return fig



def waiting_time_by_node(data, nodeID, finish, threshold):
    
    start = 0
    
    all_results = data[1].waiting_time
    if len(data)>1:
        for i in range(2,len(data)):
            #all_results = all_results.append(data[i].waiting_time)
            all_results = pandas.concat([all_results, data[i].waiting_time])
              
    sel_records = all_results[all_results["node"]==nodeID]
    sel_records = sel_records[sel_records["time_joined_queue"] >= start + Params.warm_up_duration]
    sel_records = sel_records[sel_records["time_moved_on"] < finish + Params.warm_up_duration]
    
    waiting_times = list(sel_records["waiting_time"])
     
    if len(waiting_times) > 0:
        x_max = max( threshold, max( waiting_times ) )
        
        counts, bins = numpy.histogram(waiting_times,numpy.arange(start=0,stop=math.ceil(x_max)+0.1,step=1),density=True)
        y_max = max(counts)
    else:
        x_max = threshold
        y_max = 1
        
    x_max = 50
    y_max = 1

#    if nodeID == "Triage":
#        y_max = 0.5
#    if nodeID == "Pathology":
#        y_max = 0.8
#    if nodeID == "MDT":
#        y_max = 0.6
    
    fig, ((ax1, ax2), (ax3, ax4)) = matplotlib.pyplot.subplots(2, 2, figsize=(10,10))
    fig.suptitle("Waiting time distribution at " + nodeID, fontsize = 24)
    fig.tight_layout(pad=4)
    
    ax1.hist(waiting_times,density=True,bins=numpy.arange(start=0,stop=math.ceil(x_max)+0.1,step=1))
    ax1.set_xlabel("Waiting time (days)")
    ax1.set_ylabel("Proportion of patients")
    ax1.set_xlim(0,x_max*1.1)
    ax1.set_ylim(0,y_max*1.1)
    
    #y_vals = ax1.get_yticks()
    #ax1.set_yticklabels([int(round(y/Params.number_of_runs,0)) for y in y_vals])
    
    
    time_steps = range( start, finish + 1 )
    
    sel_results = {}
    
    
    for pat_type in list(data[1].cumulative_patients[nodeID].keys()):
        
        sel_results[pat_type] = pandas.DataFrame()
        
        for r in range(1,len(data)+1):
            
            df = data[r].cumulative_patients[nodeID][pat_type]
            
            values = list(df["cumul_pat_gone"][df["time"]>=start+Params.warm_up_duration][df["time"]<=finish+Params.warm_up_duration])
            
            offset = float(df["cumul_pat_gone"][df["time"]==start+Params.warm_up_duration])
            
            sel_results[pat_type][r] = [ x - offset for x in values ]
        
        summary = pandas.DataFrame()
        summary["time"] = list(time_steps)
        summary["avg"] = list(sel_results[pat_type].mean(axis=1))
        summary["std_err"] = list(sel_results[pat_type].sem(axis=1))
        
        ax3.plot(time_steps, list(summary["avg"]), label=pat_type)
        #matplotlib.pyplot.step(time_steps, list(summary["avg"]), label=pat_type) #use this if a "step" drawstyle is preferred
        ax3.fill_between(time_steps, [summary["avg"][i]-summary["std_err"][i] for i in range(len(time_steps))], [summary["avg"][i]+summary["std_err"][i] for i in range(len(time_steps))], alpha=0.2)
        
    ax3.set_title("Cumulative patients processed at " + nodeID)
    ax3.set_xlabel("Day")
    ax3.set_ylabel("Number of patients")
    ax3.legend(loc=2)
    
    ax2.axis("off")
    ax4.axis("off")
    
    return fig



def compare_time_points_by_node(data, nodeID, start, finish, threshold):
    
    all_results = data[1].waiting_time
    if len(data)>1:
        for i in range(2,len(data)):
            #all_results = all_results.append(data[i].waiting_time)
            all_results = pandas.concat([all_results, data[i].waiting_time])
              
    sel_records = all_results[all_results["node"]==nodeID]
    sel_records = sel_records[sel_records["time_joined_queue"] < start + Params.warm_up_duration]
    sel_records = sel_records[sel_records["time_moved_on"] >= start + Params.warm_up_duration]
    
    wait_so_far_1 = [ ( start + Params.warm_up_duration - x ) for x in list(sel_records["time_joined_queue"]) ]
     
    sel_records = all_results[all_results["node"]==nodeID]
    sel_records = sel_records[sel_records["time_joined_queue"] < finish + Params.warm_up_duration]
    sel_records = sel_records[sel_records["time_moved_on"] >= finish + Params.warm_up_duration]
    
    wait_so_far_2 = [ ( finish + Params.warm_up_duration - x ) for x in list(sel_records["time_joined_queue"]) ]
    
    if len(wait_so_far_1 + wait_so_far_2) > 0:
        x_max = max( threshold, max( wait_so_far_1 + wait_so_far_2 ) )
        
        counts, bins = numpy.histogram(wait_so_far_1+wait_so_far_2,numpy.arange(start=0,stop=math.ceil(x_max)+0.1,step=1),density=True)
        y_max = max(counts)
    else:
        x_max = threshold + 1
        y_max = 1
        
    x_max = 50
    y_max = 1
    
#    if nodeID == "Triage":
#        y_max = 0.5
#    if nodeID == "Pathology":
#        y_max = 0.8
#    if nodeID == "MDT":
#        y_max = 0.6

    
    fig, ((ax3, ax4), (ax1, ax2)) = matplotlib.pyplot.subplots(2, 2, figsize=(10,10))
    fig.suptitle("Waiting time comparison at " + nodeID, fontsize = 24)
    fig.tight_layout(pad=4)
    
    ax1.hist(wait_so_far_1,density=True,bins=numpy.arange(start=0,stop=math.ceil(x_max)+0.1,step=1))
    ax1.set_title("Day " + str(start))
    ax1.set_xlabel("Waiting time (days)")
    ax1.set_ylabel("Proportion of patients")
    ax1.set_xlim(0,x_max*1.1)
    ax1.set_ylim(0,y_max*1.1)
    
    #y_vals = ax1.get_yticks()
    #ax1.set_yticklabels([int(round(y/Params.number_of_runs,0)) for y in y_vals])
    
    ax2.hist(wait_so_far_2,density=True,bins=numpy.arange(start=0,stop=math.ceil(x_max)+0.1,step=1))
    ax2.set_title("Day " + str(finish))
    ax2.set_xlabel("Waiting time (days)")
    ax2.set_ylabel("Proportion of patients")
    ax2.set_xlim(0,x_max*1.1)
    ax2.set_ylim(0,y_max*1.1)
    
    #y_vals = ax2.get_yticks()
    #ax2.set_yticklabels([int(round(y/Params.number_of_runs,0)) for y in y_vals])

    
    
    time_steps = range( start, finish + 1 )
    
    sel_results = {}
    
    
    for pat_type in list(data[1].cumulative_patients[nodeID].keys()):
        
        sel_results[pat_type] = pandas.DataFrame()
        
        for r in range(1,len(data)+1):
            
            df = data[r].cumulative_patients[nodeID][pat_type]
            
            values = list(df["cumul_pat_gone"][df["time"]>=start+Params.warm_up_duration][df["time"]<=finish+Params.warm_up_duration])
            
            offset = float(df["cumul_pat_gone"][df["time"]==start+Params.warm_up_duration])
            
            sel_results[pat_type][r] = [ x - offset for x in values ]
        
        summary = pandas.DataFrame()
        summary["time"] = list(time_steps)
        summary["avg"] = list(sel_results[pat_type].mean(axis=1))
        summary["std_err"] = list(sel_results[pat_type].sem(axis=1))
        
        ax3.plot(time_steps, list(summary["avg"]), label=pat_type)
        #matplotlib.pyplot.step(time_steps, list(summary["avg"]), label=pat_type) #use this if a "step" drawstyle is preferred
        ax3.fill_between(time_steps, [summary["avg"][i]-summary["std_err"][i] for i in range(len(time_steps))], [summary["avg"][i]+summary["std_err"][i] for i in range(len(time_steps))], alpha=0.2)
        
    ax3.set_title("Cumulative patients processed at " + nodeID)
    ax3.set_xlabel("Day")
    ax3.set_ylabel("Number of patients")
    ax3.legend(loc=2)
    
    ax4.axis("off")
    
    return fig






        
        
# Function reading input from Excel spreadsheet
def read_input(file_location):
    
    #Read and store simulation parameters
    simul_pars = pandas.read_excel(file_location+"/simulation_parameters.xlsx",header=None)
    
    names = [x for x in simul_pars[0]]
    values = [x for x in simul_pars[1]]
    
    Params.number_of_runs = values[names.index("Number of runs")]
    Params.sim_duration = values[names.index("Simulation length")]
    Params.second_time_point = values[names.index("Time at final state")]
    Params.sojourn_time_threshold = values[names.index("Target time-to-diagnosis")]
    
    if(values[names.index("Use initial state")]=="Yes"):
        Params.use_initial_patients = True
        Params.first_time_point = values[names.index("Time at initial state")]
        
        if(values[names.index("Initial state mode")]=="generate initial patients"):
            Params.generate_initial_patients = True
            Params.initial_state_source = file_location+"/initial_state_parameters.xlsx"
        else:
            Params.initial_state_source = "initial_patients.xlsx"
    else:
        Params.warm_up_duration = values[names.index("Warm up period")]
        
    if(values[names.index("Plot overall time to diagnosis")]=="Yes"):
        Params.plot_overall_sojourn_time = True
        
    if(values[names.index("Plot waiting time and activity by step")]=="Yes"):
        Params.plot_each_node = True
        
    if(values[names.index("Plot time-point comparisons")]=="Yes" and Params.use_initial_patients):
        Params.plot_comparison = True
    

    #Read and store information about patient arrivals and flow, as well as any branching based on test sensitivity/specificity
    flow_data = pandas.read_excel(file_location+"/patient_flows.xlsx", sheet_name=None)
    
    res = {}

    patient_types = list(flow_data.keys())
    res["patient_types"] = patient_types
    
    nodes_from = []
    for row in list(flow_data[patient_types[0]]["ToFrom"]):
        if(pandas.isna(row)):
            break
        else:
            nodes_from.append(row)
    res["nodes_from"] = nodes_from
    
    nodes_to = []
    for col in list(flow_data[patient_types[0]].columns.values)[1:]:
        if("Unnamed" in col):
            break
        else:
            nodes_to.append(col)
    res["nodes_to"] = nodes_to

    entry_nodes = [x for x in nodes_from if x not in nodes_to]
    res["entry_nodes"] = entry_nodes
    
    internal_nodes = [x for x in nodes_from if x in nodes_to]
    res["internal_nodes"] = internal_nodes
    
    exit_nodes = [x for x in nodes_to if x not in nodes_from]
    res["exit_nodes"] = exit_nodes


    arrival_distributions = {}
    for node in entry_nodes:
        arrival_distributions[node] = {}
        for pt in patient_types:
            tmp = list(flow_data[pt]["ArrivalDistribution"])
            ind = list(flow_data[pt]["ToFrom"]).index(node)
            arrival_distributions[node][pt] = tmp[ind]
    res["arrival_distributions"] = arrival_distributions

    arrival_parameters = {}
    for node in entry_nodes:
        arrival_parameters[node] = {}
        for pt in patient_types:
            tmp = list(flow_data[pt]["ArrivalParameters"])
            ind = list(flow_data[pt]["ToFrom"]).index(node)
            arrival_parameters[node][pt] = [float(x) for x in str(tmp[ind]).split(';')]
    res["arrival_parameters"] = arrival_parameters
    
    proportion_positive = {}
    for node in entry_nodes:
        proportion_positive[node] = {}
        for pt in patient_types:
            tmp = list(flow_data[pt]["ProportionPositive"])
            ind = list(flow_data[pt]["ToFrom"]).index(node)
            proportion_positive[node][pt] = float(tmp[ind])
    res["proportion_positive"] = proportion_positive
    
    network_structure = {}
    for pt in patient_types:
        network_structure[pt] = {}
        for node_from in (entry_nodes + internal_nodes):
            network_structure[pt][node_from] = {}
            for node_to in (internal_nodes + exit_nodes):
                tmp = list(flow_data[pt][node_to])
                ind = list(flow_data[pt]["ToFrom"]).index(node_from)
                network_structure[pt][node_from][node_to] = tmp[ind]
    res["network_structure"] = network_structure

    test_based_branching = {}
    for pt in patient_types:
        test_based_branching[pt] = {}
        for node_from in (entry_nodes + internal_nodes):
            test_based_branching[pt][node_from] = {}
            
            tmp = list(flow_data[pt]["TestBasedBranching"])
            ind = list(flow_data[pt]["ToFrom"]).index(node_from)
            
            nodes_to = [x for x in str(tmp[ind]).split(';')]
            if len(nodes_to) > 1:
                test_based_branching[pt][node_from]["positive_outcome"] = nodes_to[0]
                test_based_branching[pt][node_from]["negative_outcome"] = nodes_to[1]
                test_based_branching[pt][node_from]["sensitivity"] = float(list(flow_data[pt]["Sensitivity"])[ind])
                test_based_branching[pt][node_from]["specificity"] = float(list(flow_data[pt]["Specificity"])[ind])
            else:
                test_based_branching[pt][node_from]["positive_outcome"] = None
                test_based_branching[pt][node_from]["negative_outcome"] = None
                test_based_branching[pt][node_from]["sensitivity"] = None
                test_based_branching[pt][node_from]["specificity"] = None
    res["test_based_branching"] = test_based_branching


    #Read and store capacity data
    capacity_data = pandas.read_excel(file_location+"/capacities.xlsx")
    
    node_types = {}
    for node in internal_nodes:
        tmp = list(capacity_data["ServiceType"])
        ind = list(capacity_data["Step"]).index(node)
        node_types[node] = tmp[ind]
    res["node_types"] = node_types

    capacity_parameters = {}
    for node in internal_nodes:
        tmp = list(capacity_data["CapacityParameters"])
        ind = list(capacity_data["Step"]).index(node)
        capacity_parameters[node] = [float(x) for x in str(tmp[ind]).split(';')]
    res["capacity_parameters"] = capacity_parameters

    service_distributions = {}
    for node in internal_nodes:
        tmp = list(capacity_data["ServiceDistribution"])
        ind = list(capacity_data["Step"]).index(node)
        service_distributions[node] = tmp[ind]
    res["service_distributions"] = service_distributions
    
    service_parameters = {}
    for node in internal_nodes:
        tmp = list(capacity_data["ServiceParameters"])
        ind = list(capacity_data["Step"]).index(node)
        service_parameters[node] = [float(x) for x in str(tmp[ind]).split(';')]
    res["service_parameters"] = service_parameters

    dna_rates = {}
    for node in internal_nodes:
        tmp = list(capacity_data["DNArate"])
        ind = list(capacity_data["Step"]).index(node)
        dna_rates[node] = float(tmp[ind])
    res["dna_rates"] = dna_rates
    
    
    #Read and store information on initial patients  
    if(Params.use_initial_patients):
        
        if(Params.generate_initial_patients):
            
            initial_state_parameters = pandas.read_excel(file_location+"/initial_state_parameters.xlsx",sheet_name=None)
            
            initial_number_distributions = {}
            for node in internal_nodes:
                initial_number_distributions[node] = {}
                for pt in patient_types:
                    tmp = list(initial_state_parameters[pt]["InitialNumberDistribution"])
                    ind = list(initial_state_parameters[pt]["Step"]).index(node)
                    initial_number_distributions[node][pt] = tmp[ind]
            res["initial_number_distributions"] = initial_number_distributions
            
            initial_number_parameters = {}
            for node in internal_nodes:
                initial_number_parameters[node] = {}
                for pt in patient_types:
                    tmp = list(initial_state_parameters[pt]["InitialNumberParameters"])
                    ind = list(initial_state_parameters[pt]["Step"]).index(node)
                    initial_number_parameters[node][pt] = [float(x) for x in str(tmp[ind]).split(';')]
            res["initial_number_parameters"] = initial_number_parameters
            
            initial_node_wait_distributions = {}
            for node in internal_nodes:
                initial_node_wait_distributions[node] = {}
                for pt in patient_types:
                    tmp = list(initial_state_parameters[pt]["NodeWaitDistribution"])
                    ind = list(initial_state_parameters[pt]["Step"]).index(node)
                    initial_node_wait_distributions[node][pt] = tmp[ind]
            res["initial_node_wait_distributions"] = initial_node_wait_distributions
            
            initial_node_wait_parameters = {}
            for node in internal_nodes:
                initial_node_wait_parameters[node] = {}
                for pt in patient_types:
                    tmp = list(initial_state_parameters[pt]["NodeWaitParameters"])
                    ind = list(initial_state_parameters[pt]["Step"]).index(node)
                    initial_node_wait_parameters[node][pt] = [float(x) for x in str(tmp[ind]).split(';')]
            res["initial_node_wait_parameters"] = initial_node_wait_parameters

            initial_additional_overall_sojourn_time_distributions = {}
            for node in internal_nodes:
                initial_additional_overall_sojourn_time_distributions[node] = {}
                for pt in patient_types:
                    tmp = list(initial_state_parameters[pt]["AdditionalOverallWaitDistribution"])
                    ind = list(initial_state_parameters[pt]["Step"]).index(node)
                    initial_additional_overall_sojourn_time_distributions[node][pt] = tmp[ind]
            res["initial_additional_overall_sojourn_time_distributions"] = initial_additional_overall_sojourn_time_distributions
            
            initial_additional_overall_sojourn_time_parameters = {}
            for node in internal_nodes:
                initial_additional_overall_sojourn_time_parameters[node] = {}
                for pt in patient_types:
                    tmp = list(initial_state_parameters[pt]["AdditionalOverallWaitParameters"])
                    ind = list(initial_state_parameters[pt]["Step"]).index(node)
                    initial_additional_overall_sojourn_time_parameters[node][pt] = [float(x) for x in str(tmp[ind]).split(';')]
            res["initial_additional_overall_sojourn_time_parameters"] = initial_additional_overall_sojourn_time_parameters
            
            proportion_positive_initial = {}
            for node in internal_nodes:
                proportion_positive_initial[node] = {}
                for pt in patient_types:
                    tmp = list(initial_state_parameters[pt]["ProportionPositive"])
                    ind = list(initial_state_parameters[pt]["Step"]).index(node)
                    proportion_positive_initial[node][pt] = float(tmp[ind])
            res["proportion_positive_initial"] = proportion_positive_initial
            
        else:
            
            initial_patients = pandas.read_excel(file_location+"/initial_patients.xlsx")
            
            res["initial_patients"] = initial_patients
            

    return res



#Run simulation


input_location = Params.analysis_location + "/" + Params.analysisID + "/input"
output_location = Params.analysis_location + "/" + Params.analysisID + "/output"
model_parameters = read_input(input_location)

raw_data = {}

#print(model_parameters)

for r in range(Params.number_of_runs):
    print (f"Run {r+1} of {Params.number_of_runs}")
    my_model = PathwayModel(r, model_parameters, Params.warm_up_duration, Params.sim_duration)
    my_model.run(r)
    raw_data[r+1] = my_model.output
    #my_model.output.print(output_file, my_model.warm_up_duration, my_model.warm_up_duration + my_model.sim_duration - 1)
    print()



if Params.generate_initial_patients and Params.print_generated_patient_list:
    Params.initial_patient_list.to_excel(output_location + "/" + Params.analysisID + "_initial_patients_example.xlsx", index = False)


with PdfPages(output_location + "/" + Params.analysisID + "_report.pdf") as pdf:

    if Params.plot_overall_sojourn_time:
        pdf.savefig(time_to_diagnosis_plot(raw_data, Params.second_time_point, Params.sojourn_time_threshold))

    if Params.plot_each_node:
        for node_iter in my_model.internal_node_ids:
            pdf.savefig(waiting_time_by_node(raw_data, node_iter, Params.second_time_point, Params.sojourn_time_threshold))

    if Params.plot_comparison:
        #pdf.savefig(compare_time_points_overall(raw_data, Params.first_time_point, Params.second_time_point, Params.sojourn_time_threshold))
        for node_iter in my_model.internal_node_ids:
            pdf.savefig(compare_time_points_by_node(raw_data, node_iter, Params.first_time_point, Params.second_time_point, Params.sojourn_time_threshold))
    
    
if Params.plot_overall_sojourn_time:
    time_to_diagnosis_plot(raw_data, Params.second_time_point, Params.sojourn_time_threshold).savefig(output_location + "/" + Params.analysisID + "_time_on_pathway.png",dpi=1000)
    
