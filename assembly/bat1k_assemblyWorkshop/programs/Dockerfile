# start from base
FROM condaforge/mambaforge:4.9.2-5
SHELL ["/bin/bash", "-c"]

MAINTAINER Martin Pippel <pippel@mpi-cbg.de>

ARG NCPUS=1
ARG PREFIX="/usr/local"
ENV BINDIR="$PREFIX/bin"
ENV PATH="$BINDIR:$PATH"

# Configure tzdata (timezone)
ENV TZ=Europe/London
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Fix error message about locales (at least sometimes)
ENV LANGUAGE=C \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

# Provide our convenient build script to reduce verbosity
COPY ./build-and-install.sh /opt/

# Override default shell and use bash
SHELL ["/bin/bash", "--login", "-c"]

COPY conda-linux-64.lock .
RUN mamba create -y --prefix /usr/local --file conda-linux-64.lock && conda clean -afy

# FastK
RUN BUILD_DEPS="build-essential git" && \
    apt-get update && apt-get install -y $BUILD_DEPS && \
    REPO=https://github.com/thegenemyers/FASTK.git \
    BRANCH=9be2fe404d294d05ccbe28cca4b1279251c84123 \
    /opt/build-and-install.sh FastK make && \
    apt-get remove -y $BUILD_DEPS && apt-get autoremove -y && apt-get clean

# GENESCOPE.FK
ENV REPO=https://github.com/thegenemyers/GENESCOPE.FK.git
ENV BRANCH=380815c420f50171f9234a0fd1ff426b39829b91
RUN BUILD_DEPS="build-essential git" && \
    apt-get update && apt-get install -y $BUILD_DEPS && \
    git clone "$REPO" /opt/GENESCOPE.FK && \ 
    cd /opt/GENESCOPE.FK && \
    git checkout $BRANCH && \
#    mkdir ~/R_libs && echo "R_LIBS=~/R_libs/" >> ~/.Renviron && \ 
#    Rscript install.R && cp GeneScopeFK.R "$BINDIR" && \  
    echo "install.packages('minpack.lm', repos = \"http://cran.us.r-project.org\")" > my_install.R && \
    echo "install.packages('argparse', repos = \"http://cran.us.r-project.org\")" >> my_install.R && \
    echo "install.packages('.', repos=NULL, type=\"source\")" >> my_install.R && \
    Rscript my_install.R && cp GeneScopeFK.R "$BINDIR" && \ 
    apt-get remove -y $BUILD_DEPS && apt-get autoremove -y && apt-get clean && rm -rf /opt/GENESCOPE.FK 

# MERQURY.FK
RUN BUILD_DEPS="build-essential git" && \
    apt-get update && apt-get install -y $BUILD_DEPS && \
    REPO=https://github.com/thegenemyers/MERQURY.FK.git \
    BRANCH=b66c9c8a2c4cd82e1c8a15635eaddfae901f11dd \
    /opt/build-and-install.sh MERQURY.FK make && \
    apt-get remove -y $BUILD_DEPS && apt-get autoremove -y && apt-get clean

# manualCurationHiC 
RUN BUILD_DEPS="build-essential git" && \
    apt-get update && apt-get install -y $BUILD_DEPS && \
    REPO=https://git.mpi-cbg.de/assembly/programs/manualcurationhic.git \
    BRANCH=6be7829b72a9808d80b45e0eecaef48ef4b45b7a \
    /opt/build-and-install.sh manualcurationhic make && \
    apt-get remove -y $BUILD_DEPS && apt-get autoremove -y && apt-get clean

# Override default shell and use bash
# SHELL ["conda", "run", "-n", "assembly_workshop", "/bin/bash", "-c"]

# Python program to run in the container
COPY assembly_workshop.sh .
ENTRYPOINT ["/bin/bash", "assembly_workshop.sh"]
