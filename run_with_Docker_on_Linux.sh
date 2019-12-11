#!/bin/bash

docker run -it -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix -e LD_LIBRARY_PATH=/usr/local/lib/ rondine /opt/Rondine/Project
 
