#!/bin/bash

DOCKER_IMG='beangoben/chemical_vae'
# port for jupyter
PORT=8892
SCRIPT=$1

echo "runing docker env in"
echo "port: " ${PORT} 
echo "with " ${SCRIPT}
docker run -p $PORT:8888 -v "$(pwd)":/home/jovyan/work -it $DOCKER_IMG start.sh /home/jovyan/work/${SCRIPT}
