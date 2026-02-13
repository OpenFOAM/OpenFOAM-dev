FROM public.ecr.aws/ubuntu/ubuntu:24.04_stable AS atomic-openfoam

ARG USERNAME=ubuntu
ARG GROUPNAME=ubuntu

RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        git \
        ca-certificates \
        flex \
        libfl-dev \
        wget \
        software-properties-common \
	&& rm -rf /var/lib/apt/lists/*

# Install OpenFOAM build dependencies
RUN wget -O - https://dl.openfoam.org/gpg.key > /etc/apt/trusted.gpg.d/openfoam.asc \
	&& add-apt-repository http://dl.openfoam.org/ubuntu \
	&& apt-get update \
	&& apt-get install -y \
        openfoam-deps \
        paraview-dev \
        libscotch-dev \
        libptscotch-dev \
        libopenmpi-dev \
	&& rm -rf /var/lib/apt/lists/*

COPY . /atomic/OpenFOAM/OpenFOAM-dev/
WORKDIR /atomic/OpenFOAM
RUN git clone --depth 1 https://github.com/OpenFOAM/ThirdParty-dev.git

# Set up environment and compile OpenFOAM
RUN echo 'source /atomic/OpenFOAM/OpenFOAM-dev/etc/bashrc' >> ~/.bashrc

RUN cp /atomic/OpenFOAM/OpenFOAM-dev/build.sh .
RUN chmod +x /atomic/OpenFOAM/build.sh
RUN /atomic/OpenFOAM/build.sh
