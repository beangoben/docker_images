## Executing docker image with python script

1. Create a bash script for example **'runme.sh'** calling your script via python inside.
2. Use execute_script.sh as `./execute_script.sh runme.sh`and it will create a docker environment and launch your script.

## How to use

To run the software on any computer you need to install [docker](https://www.docker.com/).

Then you can either download or build the docker image.

To download running the following command in your favorite terminal:

```
docker pull beangoben/chemical_vae
```

of build it (good to change things) by moving to the git cloned repository :

```
docker build -t "beangoben/chemical_vae" .
```

And then move to whatever folder you want to work with and execute:

```
docker run -p 8888:8888 -v "$(pwd)":/home/jovyan/work -it beangoben/chemical_vae
```
