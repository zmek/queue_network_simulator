


import pandas as pd
import os
import simpy as sp
import numpy as np
import random as rnd
import math
import matplotlib.pyplot as mpl



#Plots
def plot_waiting_time(data, nodeID, start, finish, analysisID, analysis_location, time_transformation, bins):
    
    file_name = analysis_location + "/" + analysisID + "/output" + "/" + analysisID + "_" + nodeID + "_waiting_time.png"
    
    all_results = data[1].waiting_time
    if len(data)>1:
        for i in range(2,len(data)):
            all_results = pd.concat([all_results, data[i].waiting_time])
                    
    sel_records = all_results[all_results["node"]==nodeID]
    sel_records = sel_records[sel_records["time_moved_on"] >= start - 1 ]
    sel_records = sel_records[sel_records["time_moved_on"] < finish ]
    
    x_max = max(sel_records["waiting_time"]) * time_transformation
    y_max = len(sel_records)
    
    avg = np.mean([x * time_transformation for x in sel_records["waiting_time"]])
    print("Expected waiting time at",nodeID,":",f"{avg:.2f}")
    
    plot = mpl.figure()
    mpl.title("Waiting time distribution at " + nodeID + " - Exp " + str(round(avg,2)) )
    mpl.xlabel("Waiting time")
    mpl.ylabel("Frequency")
    
    mpl.hist([x * time_transformation for x in sel_records["waiting_time"]],bins=bins,density=False)
    mpl.axvline(x = avg, color = 'b')
    mpl.xlim(0,x_max)
    #mpl.ylim(0,y_max)
    
    mpl.show()
    plot.savefig(file_name, dpi=1000)
    
    
def plot_wait_over_time(data, start, finish, analysisID, analysis_location, time_transformation):
    
    file_name = analysis_location + "/" + analysisID + "/output" + "/" + analysisID + "_overall_wait_over_time.png"
    
    all_results = data[1].time_in_system
    if len(data)>1:
        for i in range(2,len(data)):
            all_results = pd.concat([all_results, data[i].time_in_system])
    
    sel_records = all_results
    sel_records = sel_records[sel_records["departure_time"] >= (start - 1) ]
    sel_records = sel_records[sel_records["departure_time"] < finish ]
    
    y_max = max(sel_records["total_waiting_time"]) * time_transformation
    
    x_vals = [ x for x in range(start,finish) ]
    y_vals = [ 0 for y in range(start,finish) ]
    
    for x in range(len(x_vals)):
        sel = sel_records[sel_records["departure_time"] >= (x_vals[x]-1)]
        sel = sel[sel["departure_time"] < x_vals[x] ]
        if len(sel)>0:
            y_vals[x] = np.mean(sel["total_waiting_time"])
            
    plot = mpl.figure()
    mpl.title("Overall waiting time, depending on departure time")
    mpl.xlabel("Departure time")
    mpl.ylabel("Waiting time")
    
    mpl.scatter(x_vals * time_transformation, y_vals * time_transformation)
    mpl.xlim(start,finish)
    mpl.ylim(0,y_max)
    mpl.show()
    plot.savefig(file_name, dpi=1000)



def plot_overall_time_in_system(data, start, finish, analysisID, analysis_location, time_transformation, bins):
    
    file_name = analysis_location + "/" + analysisID + "/output" + "/" + analysisID + "_overall_time_in_system.png"
    
    all_results = data[1].time_in_system
    if len(data)>1:
        for i in range(2,len(data)):
            all_results = pd.concat([all_results, data[i].time_in_system])
            
    sel_records = all_results
    sel_records = sel_records[sel_records["departure_time"] >= start - 1 ]
    sel_records = sel_records[sel_records["departure_time"] < finish ]
    
    x_max = max(sel_records["time_in_system"]) * time_transformation
    y_max = len(sel_records)
    
    avg = np.mean([x * time_transformation for x in sel_records["time_in_system"]])
    print("Expected time spent in the system:",f"{avg:.2f}")
    
    plot = mpl.figure()
    mpl.title("Overall time in system distribution" + " - Exp " + str(round(avg,2)) )
    mpl.xlabel("Time in system")
    mpl.ylabel("Frequency")
    
    mpl.hist([x * time_transformation for x in sel_records["time_in_system"]],bins=bins,density=False)
    mpl.axvline(x = avg, color = 'b')
    mpl.xlim(0,x_max)
    #mpl.ylim(0,y_max)
    mpl.show()
    plot.savefig(file_name, dpi=1000)
    
    
def plot_performance(data, target, start, finish, analysisID, analysis_location, time_transformation, bins):
    
    file_name = analysis_location + "/" + analysisID + "/output" + "/" + analysisID + "_performance.png"
    
    all_results = data[1].time_in_system
    if len(data)>1:
        for i in range(2,len(data)):
            all_results = pd.concat([all_results, data[i].time_in_system])
            
    sel_records = all_results
    sel_records = sel_records[sel_records["departure_time"] >= start - 1 ]
    sel_records = sel_records[sel_records["departure_time"] < finish ]
    
    x_max = max(sel_records["time_in_system"]) * time_transformation
    y_max = len(sel_records)
    
    within_target = sel_records[sel_records["time_in_system"]<=target]
    
    perf = len(within_target) / len(sel_records)
    
    print("Proportion of patients leaving within target time:",f"{perf*100:.2f}","%")
    
    plot = mpl.figure()
    mpl.title("Proportion of patients leaving within target time:" + str(round(perf*100,2)) + "%")
    mpl.xlabel("Time in system")
    mpl.ylabel("Frequency")
    
    mpl.hist([x * time_transformation for x in sel_records["time_in_system"]],bins=bins,density=False)
    mpl.axvline(x = target, color = 'k')
    mpl.xlim(0,x_max)
    #mpl.ylim(0,y_max)
    mpl.show()
    plot.savefig(file_name, dpi=1000)
    

def plot_total_waiting_time(data, start, finish, analysisID, analysis_location, time_transformation, bins):
    
    file_name = analysis_location + "/" + analysisID + "/output" + "/" + analysisID + "_overall_waiting_time.png"
    
    all_results = data[1].time_in_system
    if len(data)>1:
        for i in range(2,len(data)):
            all_results = pd.concat([all_results, data[i].time_in_system])
    
    sel_records = all_results
    sel_records = sel_records[sel_records["departure_time"] >= start - 1 ]
    sel_records = sel_records[sel_records["departure_time"] < finish ]
    
    x_max = max(sel_records["total_waiting_time"]) * time_transformation
    y_max = len(sel_records)
    
    avg = np.mean([x * time_transformation for x in sel_records["total_waiting_time"]])
    print("Expected overall waiting time:",f"{avg:.2f}")
    
    plot = mpl.figure()
    mpl.title("Overall waiting time distribution" + " - Exp " + str(round(avg,2)) )
    mpl.xlabel("Waiting time")
    mpl.ylabel("Frequency")
    
    mpl.hist([x * time_transformation for x in sel_records["total_waiting_time"]],bins=bins,density=False)
    mpl.axvline(x = avg, color = 'b')
    mpl.xlim(0,x_max)
    #mpl.ylim(0,y_max)
    mpl.show()
    plot.savefig(file_name, dpi=1000)

