FROM quay.io/jupyter/datascience-notebook:2024-03-26

SHELL ["/bin/bash", "-o", "pipefail", "-e", "-u", "-x", "-c"]

ARG WORK_VOLUME_MOUNT_PATH
ARG GITHUB_REPO_URL=https://github.com/zmek/queue_network_simulator.git
ARG GITHUB_REPO_BRANCH=main

USER root

# Extra Linux tools that are nice to have.
RUN apt-get update --yes && \
    apt-get install --yes --no-install-recommends sqlite3=3.37.2-2ubuntu0.3 git && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN pip install --no-cache-dir pip==24.0 && \
    pip install --no-cache-dir pre-commit==3.7.0 nbstripout==0.7.1

# Add the requirements.txt file
COPY requirements.txt /tmp/

# Install Python dependencies from requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt

# Clone the GitHub repository
RUN git clone --branch ${GITHUB_REPO_BRANCH} ${GITHUB_REPO_URL} ${WORK_VOLUME_MOUNT_PATH}

RUN mkdir -p ${WORK_VOLUME_MOUNT_PATH}
RUN chown -R ${NB_UID}:${NB_GID} ${WORK_VOLUME_MOUNT_PATH}

# Switch back to jovyan to avoid accidental container runs as root
USER ${NB_UID}



