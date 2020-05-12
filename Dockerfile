# creates an image containing vg and hal2vg

# build on compatible vg image
FROM quay.io/vgteam/vg:v1.24.0

# update system and install dependencies not present in vg image
RUN apt-get -qq update && apt-get -qq install -y libhdf5-dev build-essential python3-dev python3-pip

# copy current directory to docker
ADD . /hal2vg

# set working directory
WORKDIR /hal2vg

# build
RUN make

# add hal2vg to the PATH
ENV PATH /hal2vg:/hal2vg/deps/hal/bin:$PATH
