# queue_network_simulator
Python code to demonstrate Discrete Event Simulation, developed by Luca Grieco, with Jupyter lab set-up created by Zella King, both of the Clinical Operational Research Unit at UCL

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/zmek/queue_network_simulator/HEAD?labpath=notebooks%2Ftest.ipynb)


# Multi-Instance Jupyter Lab Setup with Docker Compose

This repository includes configuration files to set up multiple instances of Jupyter Lab using Docker Compose. Each instance will run on a different port with its own environment settings.

## Prerequisites

- Docker: [Install Docker](https://docs.docker.com/get-docker/)
- Docker Compose: [Install Docker Compose](https://docs.docker.com/compose/install/)

## Setup Instructions

### Step 1: Clone the Repository

Clone this repository to your local machine:

```sh
git clone https://github.com/zmek/queue_network_simulator.git
cd your-repo
```

### Step 2: Create environment file in the root, calling it ```.env```. If creating multiple instances, save each one with a different name ```.env.instance1``` for instance1
```sh
JUPYTER_PASSWORD=instance1_password
JUPYTER_EXPOSED_PORT=8889
```

### Step 3: Can skip if not creating multiple instances. 

#### Option 1: Create Docker Compose Override Files

Create a Docker Compose files for each instance, by copying the existing one. Example content is provided below:


```yaml
version: '3.8'

services:
  jupyter:
    container_name: jupyter_instance1
    env_file:
      - .env.instance1

```

#### Option 2: Create one docker file with all services 

Step 4: Build and Start the Instances

Use Docker Compose to build and start each instance with the respective override file. -d starts the output as a daemon. Skip this if you want to see the output. 

For a single instance, either of these will suffice

``` sh
docker compose -f docker-compose.yml up --build
docker-compose up --build

```

For multiple instances option 1, replace the name of the docker-compose file as appropriate, and the env file. Including the project name ensures that the docker system process using a different namespace for each container. 

``` sh
docker compose -f docker-compose.instance1.yml --env-file .env.instance1 --project-name desworkshop_instance1 up --build -d
``` 


```sh
docker compose --env-file .env --project-name desworkshop-zella up --build -d jupyter-zella
```

### Step 5: Access Jupyter lab

Once the instances are running, you can access Jupyter Lab in your web browser using the following URLs:

Instance 1: http://servername:8889
Enter the token set by JUPYTER_PASSWORD in the respective .env file to log in.

Stopping the Instances

To stop the single instance
```sh
docker-compose down

```
To stop the running instances, use the following commands:

Stop Instance 1
```sh
docker-compose -f docker-compose.yml -f docker-compose.instance1.yml down
```
