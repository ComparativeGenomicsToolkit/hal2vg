# creates an image containing vg and hal2vg

# build on compatible vg image
FROM quay.io/vgteam/vg:v1.11.0-74-gdab42acd-t242-build

# update system and install dependencies not present in vg image
RUN apt-get -qq update && apt-get -qq install -y libhdf5-serial-dev

# copy current directory to docker
ADD . /hal2vg

# set working directory
WORKDIR /hal2vg

# build
RUN make

# add hal2vg to the PATH
ENV PATH /hal2vg:$PATH
