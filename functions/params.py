import pandas


class Params:

    analysisID = "test"  # name of the folder containing input and output folders of the current analysis

    analysis_location = "data"  # relative path from this Python script to the folder containing [analysisID]'s folder

    number_of_runs = 0  # ideally 100 or above
    warm_up_duration = 0  # days - warm-up period, to be set to 0 if we allow a non-empty system as starting state
    sim_duration = (
        0  # days - duration of the actual simulation after the warm-up period
    )

    first_time_point = 0  # i-th week after warm-up period (between 0 and sim_duration) for which we want to have a snapshot of waiting times
    second_time_point = 0  # j-th day after warm-up period (between 0 and sim_duration) for which we want to have a snapshot of waiting times
    sojourn_time_threshold = 0  # target maximum waiting time (in days)

    use_initial_patients = False
    initial_state_source = None
    generate_initial_patients = False

    plot_overall_sojourn_time = False
    plot_comparison = False
    plot_each_node = False

    # Example of list of initial patients - technical parameters
    print_generated_patient_list = False
    patient_list_run = (
        1  # simulation run from which an example of initial patient list is taken
    )
    initial_patient_list = pandas.DataFrame()
    initial_patient_list["patID"] = []
    initial_patient_list["pat_stream"] = []
    initial_patient_list["truth"] = []
    initial_patient_list["current_node"] = []
    initial_patient_list["time_waited_at_node"] = []
    initial_patient_list["additional_time_in_system"] = []
