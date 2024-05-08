import simpy


# Class representing the nodes of the network
# Nodes could either represent discrete servers or pools of resources
# Capacity is set up according to the type of resource, whereas service rates
# and replenishment levels/frequencies are defined based on patient types
class Node:

    def __init__(
        self,
        env,
        node_name,
        node_type,
        resource_type,
        capacity,
        service_distributions,
        service_parameters,
    ):
        self.id = node_name
        self.type = node_type  # allowed options: {"entry","internal","exit"}

        self.resource = None
        self.replenish_freq = float("inf")
        self.resource_type = resource_type  # allowed options: {"gradual","batch"}

        if self.type == "internal":
            if self.resource_type == "gradual":
                self.resource = simpy.PriorityResource(env, capacity=capacity[0])
            else:  # i.e. self.resource_type = "batch"
                self.resource = simpy.PriorityResource(env, capacity=capacity[0])
                self.replenish_freq = capacity[1]

        if self.type == "exit":
            self.resource = simpy.Container(env, capacity=100, init=100)

        self.service_distributions = service_distributions
        self.service_parameters = service_parameters
