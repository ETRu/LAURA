FROM sarusso/metadesktop

# Build dependencies
RUN apt-get install -y  build-essential libxml2 libx11-dev libxml2-dev libglu1-mesa-dev freeglut3-dev man libpng-dev libx11-dev libjpeg-dev libxft-dev libxinerama-dev libtool automake bison byacc flex 

# Install igraph
RUN cd /tmp && git clone https://github.com/igraph/igraph
RUN cd /tmp/igraph && git pull && git checkout 2c61d45d9b5afde2fd78126a4f0d48f819c5ac94
RUN cd /tmp/igraph && ./bootstrap.sh && ./configure && make && make install

# Install fltk
RUN cd /tmp && git clone https://github.com/fltk/fltk
RUN cd /tmp/fltk && git pull && git checkout 1195fa189b1ef1c94cb9da247f5a5bb7b1923bc1
RUN cd /tmp/fltk && ./autogen.sh && ./configure -enable-gl && make && make install

# Add code
COPY / /opt/Rondine

# Compile code
RUN cd /opt/Rondine && make

