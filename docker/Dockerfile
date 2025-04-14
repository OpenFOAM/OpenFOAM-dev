FROM ubuntu:24.04 AS openfoam

RUN apt-get update \
	&& apt-get install -y \
        build-essential \
        cmake \
        git \
        ca-certificates \
        flex \
        wget \
        software-properties-common \
        python-is-python3 \
        python3-pip \
        python3-venv

# Install OpenFOAM build dependencies
RUN wget -O - https://dl.openfoam.org/gpg.key > /etc/apt/trusted.gpg.d/openfoam.asc
RUN add-apt-repository http://dl.openfoam.org/ubuntu
RUN apt install -y openfoam-deps paraview-dev

WORKDIR /atomic/OpenFOAM
RUN git clone --single-branch --branch atomic https://github.com/Atomic-Industries/OpenFOAM-dev.git
RUN git clone https://github.com/OpenFOAM/ThirdParty-dev.git


# Set up environment and compile OpenFOAM
RUN echo 'source /atomic/OpenFOAM/OpenFOAM-dev/etc/bashrc' >> ~/.bashrc

COPY build.sh /atomic/OpenFOAM/
RUN chmod +x /atomic/OpenFOAM/build.sh
RUN /atomic/OpenFOAM/build.sh

# Set up Python
RUN python3 -m venv .venv
ENV PATH="/atomic/OpenFOAM/.venv/bin:${PATH}"
RUN python -m pip install numpy scipy jinja2 click
RUN cd $WM_PROJECT_DIR/python && python -m pip install --editable .