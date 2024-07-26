# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 16:58:39 2022

@author: shipatel
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-



import simpy
import random
import pandas as pd
import csv

# Class to store global parameter values.  
class g:
    utc_inter = 5
    mean_consult_one = 35
    mean_consult_two = 20
    mean_discharge = 15
    mean_diagnostic = 45
    prob_diag = 0.11
              
    #number_of_treatment_rooms= 7
    number_of_utc_clinician = 7
    
    sim_duration = 840
    warm_up_duration = 180
    number_of_runs = 200
    
# Class representing our patients coming in to the UTC where we store a
# patient ID and whether the patient will be sent for diagnostics
       
class UTC_Patient:
    def __init__(self, p_id, prob_diag):
        self.id = p_id
        self.prob_diag = prob_diag
        self.diag_patient = False
        # determines whether the patient requires diagnostics
        self.start_q_consult_one=0
        self.q_time_consult_one = 0
        self.duration_consult_one=0
        self.duration_diagnostics=0
        self.q_time_consult_two = 0
        self.duration_consult_two=0
        self.duration_discharge=0
        self.time_in_UTC=0
        self.flag_less4h=1
        
    # Method to determine whether or not this patient will require
    # diagnostics
    def determine_diag_destiny(self):
        if random.uniform(0, 1) < self.prob_diag:
            self.diag_patient  = True
       
# Class representing our model of the UTC
class UTC_Model:
    def __init__(self, run_number):
        self.env = simpy.Environment()
        self.patient_counter = 0
        self.utc_clinician = simpy.Resource(self.env,
                                           capacity=g.number_of_utc_clinician)
        #self.treatment_rooms = simpy.Resource(self.env, capacity=g.number_of_treatment_rooms)

        
        self.run_number = run_number
        self.mean_q_time_consult = 0
        
        self.results_df = pd.DataFrame()
        self.results_df["P_ID"] = []
        self.results_df["Diag_flag"]=[]
        self.results_df["Arrival_time"] = []
        self.results_df["Q_Time_consult_one"] = []
        self.results_df["Duration_Time_consult_one"] = []
        self.results_df["Duration_Time_Diagnostics"] = []
        self.results_df["Q_Time_consult_two"] = []
        self.results_df["Duration_Time_consult_two"] = []
        self.results_df["Duration_Time_Discharge"] = []
        self.results_df["Time_in_UTC"]=[]
        self.results_df["Flag_less4h"]=[]
        self.results_df.set_index("P_ID", inplace=True)

   # A method that generates patients arriving at the UTC
    def generate_UTC_arrivals(self):
       # Keep generating indefinitely whilst the simulation is running
       while True:
            # Increment patient counter by 1
            self.patient_counter += 1
            # Create a new patient
            p = UTC_Patient(self.patient_counter, g.prob_diag)
            
            # Determine the patient's UTC destiny by running the appropriate
            # method
            p.determine_diag_destiny()
            #p.default_priority()
            
            # Get the SimPy environment to run the UTC_patient_journey method 
            # with this patient
            self.env.process(self.UTC_patient_journey(p))
            
            # Randomly sample the time to the next patient arriving
            sampled_interarrival = random.expovariate(1.0 / g.utc_inter)
            
            # Freeze this function until that time has elapsed
            yield self.env.timeout(sampled_interarrival)
            
    def UTC_patient_journey(self, patient):
                #"""CONSULT ONE """
        # Record the time the patient started queuing for first consultation
        start_q_consult_one = self.env.now
        # Request a UTC clinician
        with self.utc_clinician.request()as req:
                        
            # Freeze the function until the request can be met
            yield req
            end_q_consult_one = self.env.now
            patient.q_time_consult_one = end_q_consult_one - start_q_consult_one
            
            #print (f"Patient {patient.id} is seen for UTC consult one after waiting {patient.q_time_consult_one:.0f} minutes")
            sampled_consult_one_duration = random.expovariate(1.0 / g.mean_consult_one)
            patient.duration_consult_one = sampled_consult_one_duration
            print (f"Patient {patient.id} is seen for UTC consult one after waiting {patient.q_time_consult_one:.0f} minutes and consult one takes {sampled_consult_one_duration:.2f} minutes")
            #print (f"Patient {patient.id} consult one takes {sampled_consult_one_duration:.2f} minutes.")
            #print (f"Patient {priority} ")
            # Freeze this function until that time has elapsed
            yield self.env.timeout(sampled_consult_one_duration)
            
             ##"""BRANCH - DIAGNOSTICS OR DISCHARGE ACTIVITY"""
         
        if patient.diag_patient == False:
            ##"" "DISHARGE ACTIVITY """
            sampled_discharge_duration = (
                random.expovariate(1.0 / g.mean_discharge))
            patient.duration_discharge = sampled_discharge_duration
            patient.time_in_UTC = (patient.q_time_consult_one+patient.duration_consult_one+ patient.duration_discharge)
            
            if patient.time_in_UTC <= 240:
                patient.flag_less4h=1
            else:
                patient.flag_less4h=0
            
            #Freeze this time until that time has elapsed
            yield self.env.timeout(sampled_discharge_duration)
            print (f"Patient {patient.id}  discharged without diagnostics after discharge ({sampled_discharge_duration:.2f} minutes) and the patient spent {patient.time_in_UTC:.2f} minutes in UTC in total.") 
            #print (f"Patient {patient.id} spent  {patient.time_in_UTC:.2f} minutes in UTC.")
        else:
            ##"""DIAGNOSTICS"""
            # the duration includes travel and queuing for diagnostics
            sampled_diagnostics_duration = (
                random.expovariate(1.0 / g.mean_diagnostic))
            patient.duration_diagnostics = sampled_diagnostics_duration
            yield self.env.timeout(sampled_diagnostics_duration)
            print (f"Patient {patient.id} returns from diagnostics after {sampled_diagnostics_duration:.2f} minutes.")
            
            
            ##"""CONSULT TWO"""
            start_q_consult_two =self.env.now
            
            # Request a UTC clinician
            with self.utc_clinician.request() as req:
            # Freeze the function until the request can be met
                yield req
                end_q_consult_two = self.env.now
                patient.q_time_consult_two = end_q_consult_two - start_q_consult_two
                #print (f"Patient {patient.id} is seen by UTC doctor for 2nd consult after waiting {patient.q_time_consult_two:.2f} minutes")
                sampled_consult_two_duration = random.expovariate(1.0 / g.mean_consult_two)
                patient.duration_consult_two = sampled_consult_two_duration
                #print (f"Patient {patient.id} is seen by UTC doctor and 2nd consult takes {sampled_consult_two_duration:.2f} minutes.")
                print (f"Patient {patient.id} is seen for UTC consult two after waiting {patient.q_time_consult_two:.0f} minutes and consult two takes {sampled_consult_two_duration:.2f} minutes")
                
                # Freeze this function until that time has elapsed
                yield self.env.timeout(sampled_consult_two_duration)
            
            ##"" "DISHARGE ACTIVITY """

            sampled_discharge_duration = (
                random.expovariate(1.0 / g.mean_discharge))
            patient.duration_discharge = sampled_discharge_duration
            patient.time_in_UTC = (patient.q_time_consult_one+patient.duration_consult_one
            +patient.q_time_consult_two+patient.duration_consult_two
            +patient.duration_diagnostics
            + patient.duration_discharge)
            
            if patient.time_in_UTC <= 240:
                patient.flag_less4h=1
            else:
                patient.flag_less4h=0
                
            #Freeze this time until that time has elapsed
            yield self.env.timeout(sampled_discharge_duration)
            print (f"Patient {patient.id} is discharged after 2nd consult and spent {patient.time_in_UTC:.2f} minutes in UTC in total.")
        # If the warm up time has passed, then call the store_patient_results 
        # method (this doesn't need to be processed by the environment, as it's
        # not a generator function)
        if self.env.now > g.warm_up_duration:
            self.store_patient_results(patient)
            
            
             ##########################################################
   # A method to store the patient's results (queuing times here) for this
    # run alongside their patient ID in the Pandas DataFrame of the ED_Model
    # class
    def store_patient_results(self, patient):        
        # First, because we have a branching path, this patient will have
        # either been discharged post-consult or gone for diagnostics 
        # and then second consult but not both.
        # Therefore, we need to check which happened, and insert NaNs
        # (Not A Number) in the entries for the other queue in the DataFrame.
        # NaNs are automatically ignored by Pandas when calculating the mean
        # etc.  We can create a nan by casting the string 'nan' as a float :
        # float("nan")
        if patient.diag_patient== False:
            patient.duration_diagnostics = float("nan")
            patient.q_time_consult_two = float("nan")
            patient.duration_consult_two = float("nan")
        else: None
        
            
        df_to_add = pd.DataFrame({"P_ID":[patient.id],"Diag_flag":[patient.diag_patient],
                                  "Arrival_time":[patient.start_q_consult_one],
                                  "Q_Time_Consult_One":[patient.q_time_consult_one],
                                  "Duration_Time_Consult_One":[patient.duration_consult_one],
                                  "Duration_Time_Diagnostics":[patient.duration_diagnostics],
                                  "Q_Time_Consult_Two":([patient.q_time_consult_two]),
                                  "Duration_Time_Consult_Two":([patient.duration_consult_two]),
                                  "Duration_Time_Discharge":([patient.duration_discharge]),
                                  "Time_in_UTC":([patient.time_in_UTC]),
                                  "Flag_less4h":([patient.flag_less4h])
                                  })   
        df_to_add.set_index("P_ID", inplace=True)
        self.results_df = pd.concat([self.results_df, df_to_add])
      
        
        
        
    # A method that calculates the average queuing times for each queue.  We
    # can call this at the end of each run
    def calculate_mean_q_times(self):
        self.patients_seen = self.results_df.shape[0]
        self.mean_q_time_consult_one = (
            self.results_df["Q_Time_Consult_One"].mean())
        self.mean_q_time_consult_two = (
            self.results_df["Q_Time_Consult_Two"].mean())
        self.mean_duration_consult_one = (
            self.results_df["Duration_Time_Consult_One"].mean())
        self.mean_duration_consult_two = (
            self.results_df["Duration_Time_Consult_Two"].mean())
        self.mean_duration_diagnostics = (
            self.results_df["Duration_Time_Diagnostics"].mean())
        self.mean_duration_discharge = (
            self.results_df["Duration_Time_Discharge"].mean())
        self.mean_time_in_UTC = (
            self.results_df["Time_in_UTC"].mean())
        self.mean_less4h = (
            self.results_df["Flag_less4h"].mean())
        
    # A method to write run results to file
    def write_run_results(self):
        with open("trial_UTC_results.csv", "a") as f:
            writer = csv.writer(f, delimiter=",")
            results_to_write = [self.run_number,
                                self.patients_seen,
                                self.mean_q_time_consult_one,
                                self.mean_q_time_consult_two,
                                self.mean_duration_consult_one,
                                self.mean_duration_consult_two,
                                self.mean_duration_diagnostics,
                                self.mean_duration_discharge,
                                self.mean_time_in_UTC,
                                self.mean_less4h]
            writer.writerow(results_to_write)
            
    # The run method starts up the entity generators, and tells SimPy to start
    # running the environment for the duration specified in the g class. After
    # the simulation has run, it calls the methods that calculate run
    # results, and the method that writes these results to file
    def run(self):
        # Start entity generators
        self.env.process(self.generate_UTC_arrivals())
                
        # Run simulation
        self.env.run(until=(g.sim_duration + g.warm_up_duration))
        
        # Calculate run results
        self.calculate_mean_q_times()
        
        # Write run results to file
        self.write_run_results()    
                
# Class to store, calculate and manipulate trial results
class Trial_Results_Calculator:
    def __init__(self):
        self.trial_results_df = pd.DataFrame()
        
    # A method to read in the trial results and print them for the user
    def print_trial_results(self):
        print ("TRIAL RESULTS")
        print ("-------------")
        
        # Read in results from each run
        self.trial_results_df = pd.read_csv("trial_UTC_results.csv")
        
        # Take average over runs
        trial_mean_patients_seen = (
            self.trial_results_df["Patients_Seen"].mean())
        trial_mean_q_time_consult_one = (
            self.trial_results_df["Mean_Q_Time_Consult_One"].mean())
        trial_mean_q_time_consult_two = (
            self.trial_results_df["Mean_Q_Time_Consult_Two"].mean())
        trial_mean_duration_consult_one = (
            self.trial_results_df["Mean_Duration_Consult_One"].mean())
        trial_mean_duration_consult_two = (
            self.trial_results_df["Mean_Duration_Consult_Two"].mean())
        trial_mean_duration_diagnostics = (
            self.trial_results_df["Mean_Duration_Diagnostics"].mean())
        trial_mean_duration_discharge = (
            self.trial_results_df["Mean_Duration_Discharge"].mean())
        trial_mean_time_in_UTC = (
            self.trial_results_df["Mean_Time_UTC"].mean())
        trial_less4h = 100*(self.trial_results_df["Mean_Less4h"].mean())
                
        #number of patients seen
        print ("Average patients discharged :",
               f"{trial_mean_patients_seen:.0f}")
        print ("Mean Queuing Time for Consult_One :",
               f"{trial_mean_q_time_consult_one:.2f}")
        print ("Mean Queuing Time for Consult_Two :",
               f"{trial_mean_q_time_consult_two:.2f}")
        print ("Mean Duration for Consult_One :",
               f"{trial_mean_duration_consult_one:.2f}")
        print ("Mean Duration for Consult_Two :",
               f"{trial_mean_duration_consult_two:.2f}")

        print ("Mean Duration for Diagnostics :",
               f"{trial_mean_duration_diagnostics:.2f}")
        print ("Mean Duration for Discharge :",
               f"{trial_mean_duration_discharge:.2f}")
        print ("Mean Time in UTC :",
              f"{trial_mean_time_in_UTC:.2f}")
        print ("Patients in UTC % < 4hours :",
              f"{trial_less4h:.2f}","%")
                
# Everything above is definition of classes and functions, but here's where
# the code will start actively doing things.        

# Create a file to store trial results
with open("trial_utc_results.csv", "w") as f:
    writer = csv.writer(f, delimiter=",")
    column_headers = ["Run",
                      "Patients_Seen",
                      "Mean_Q_Time_Consult_One",
                      "Mean_Q_Time_Consult_Two",
                      "Mean_Duration_Consult_One",
                      "Mean_Duration_Consult_Two",
                      "Mean_Duration_Diagnostics",
                      "Mean_Duration_Discharge",
                      "Mean_Time_UTC",
                      "Mean_Less4h"
                      ]
    writer.writerow(column_headers)

# For the number of runs specified in the g class, create an instance of the
# UTC_Model class, and call its run method
for run in range(g.number_of_runs):
    print (f"Run {run+1} of {g.number_of_runs}")
    my_utc_model = UTC_Model(run)
    my_utc_model.run()
    print ()

# Once the trial is complete, we'll create an instance of the
# Trial_Result_Calculator class and run the print_trial_results method
my_trial_results_calculator = Trial_Results_Calculator()
my_trial_results_calculator.print_trial_results()
            


           