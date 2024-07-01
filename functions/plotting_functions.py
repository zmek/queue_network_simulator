


import pandas as pd
import os
import simpy as sp
import numpy as np
import random as rnd
import math
import matplotlib.pyplot as mpl



#Plots
def plot_waiting_time(data, nodeID, start, finish, analysisID, analysis_location, time_transformation):
    
    file_name = analysis_location + "/" + analysisID + "/output" + "/" + analysisID + "_" + nodeID + "_waiting_times.png"
    
    all_results = data[1].waiting_time
    if len(data)>1:
        for i in range(2,len(data)):
            all_results = pd.concat([all_results, data[i].waiting_time])
            
    x_max = max(all_results["waiting_time"]) * time_transformation
           
    plot = mpl.figure()
    mpl.title("Waiting time distribution - node: " + nodeID)
    mpl.xlabel("Waiting time (hours)")
    mpl.ylabel("Proportion of patients")
    
    sel_records = all_results[all_results["node"]==nodeID]
    sel_records = sel_records[sel_records["time_moved_on"] >= start - 1 ]
    sel_records = sel_records[sel_records["time_moved_on"] < finish ]
    
    print(np.mean([x * time_transformation for x in sel_records["waiting_time"]]))
    
    mpl.hist([x * time_transformation for x in sel_records["waiting_time"]],density=True)
    #mpl.xlim(0,x_max)
    mpl.ylim(0,1)
    
    mpl.show()
    plot.savefig(file_name, dpi=1000)
    
    
def plot_wait_over_time(data, start, finish, analysisID, analysis_location, time_transformation):
    
    file_name = analysis_location + "/" + analysisID + "/output" + "/" + analysisID + "_wait_over_time.png"
    
    all_results = data[1].time_in_system
    if len(data)>1:
        for i in range(2,len(data)):
            all_results = pd.concat([all_results, data[i].time_in_system])
    
    y_max = max(all_results["time_in_system"]) * time_transformation
    
    plot = mpl.figure()
    mpl.title("Time in system over time")
    mpl.xlabel("Arrival time")
    mpl.ylabel("Time in system")
    
    sel_records = all_results
    sel_records = sel_records[sel_records["departure_time"] >= (start - 1) ]
    sel_records = sel_records[sel_records["departure_time"] < finish ]
    
    x_vals = [ x for x in range(start,finish) ]
    y_vals = [ 0 for y in range(start,finish) ]
    
    for x in range(len(x_vals)):
        sel = sel_records[sel_records["departure_time"] >= (x_vals[x]-1)]
        sel = sel[sel["departure_time"] < x_vals[x] ]
        #print(len(sel))
        if len(sel)>0:
            y_vals[x] = np.mean(sel["time_in_system"])
    
    mpl.scatter(x_vals * time_transformation, y_vals * time_transformation)
    mpl.xlim(start,finish)
    mpl.ylim(0,y_max)
    mpl.show()
    plot.savefig(file_name, dpi=1000)



def plot_time_in_system(data, start, finish, analysisID, analysis_location, time_transformation):
    
    file_name = analysis_location + "/" + analysisID + "/output" + "/" + analysisID + "_time_in_system.png"
    
    all_results = data[1].time_in_system
    if len(data)>1:
        for i in range(2,len(data)):
            all_results = pd.concat([all_results, data[i].time_in_system])
            
    x_max = max(all_results["time_in_system"]) * time_transformation
           
    plot = mpl.figure()
    mpl.title("Time in system distribution")
    mpl.xlabel("Time in system (days)")
    mpl.ylabel("Proportion of patients")
    
    sel_records = all_results
    sel_records = sel_records[sel_records["departure_time"] >= start - 1 ]
    sel_records = sel_records[sel_records["departure_time"] < finish ]
    
    print(np.mean([x * time_transformation for x in sel_records["time_in_system"]]))
    
    mpl.hist([x * time_transformation for x in sel_records["time_in_system"]],density=True)
    mpl.xlim(0,x_max)
    mpl.ylim(0,1)
    mpl.show()
    plot.savefig(file_name, dpi=1000)
    

def plot_total_waiting_time(data, start, finish, analysisID, analysis_location, time_transformation):
    
    file_name = analysis_location + "/" + analysisID + "/output" + "/" + analysisID + "_total_waiting_time.png"
    
    all_results = data[1].time_in_system
    if len(data)>1:
        for i in range(2,len(data)):
            all_results = pd.concat([all_results, data[i].time_in_system])
            
    x_max = max(all_results["time_in_system"]) * time_transformation
           
    plot = mpl.figure()
    mpl.title("Total waiting time distribution")
    mpl.xlabel("Total waiting time")
    mpl.ylabel("Proportion of patients")
    
    sel_records = all_results
    sel_records = sel_records[sel_records["departure_time"] >= start - 1 ]
    sel_records = sel_records[sel_records["departure_time"] < finish ]
    
    print(np.mean([x * time_transformation for x in sel_records["total_waiting_time"]]))
    
    mpl.hist([x * time_transformation for x in sel_records["total_waiting_time"]],density=True)
    mpl.xlim(0,x_max)
    mpl.ylim(0,1)
    mpl.show()
    plot.savefig(file_name, dpi=1000)

