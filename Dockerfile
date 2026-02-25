FROM opencfd/openfoam-dev-v2512

# Set working directory
WORKDIR /home/openfoam/project

# Copy your local repo into the container
COPY --chown=openfoam:openfoam . .

# Compile the CTT solver
RUN . /usr/lib/openfoam/openfoam2512/etc/bashrc && \
    cd applications/solvers/incompressible/cttNavierStokesFoam_CTT_Logic && \
    wmake
