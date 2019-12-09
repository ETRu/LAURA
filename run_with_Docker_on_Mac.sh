#!/bin/bash

echo "Note: you need XQuartz and to tick XQuartz -> Preferences -> Security -> Allow connections from network clients"

xhost + 127.0.0.1

docker run -it -e DISPLAY=docker.for.mac.host.internal:0 -e LD_LIBRARY_PATH=/usr/local/lib/ rondine /opt/Rondine/Project


# p.s. https://github.com/docker/for-mac/issues/2965
 
