## Executing docker image with python script

1. Create a bash script for example 'train.sh' calling your script via python inside.
2. Use execute_script.sh as `./execute_script.sh train.sh`and it will create a docker environment and launch your script.

## How to use

To run the software on any computer you need to install [docker](https://www.docker.com/).

Then you can either download or build the docker image.

To download running the following command in your favorite terminal:

```
docker pull beangoben/mol_autoencoder
```

of build it (good to change things) by moving to the git cloned repository :

```
docker build -t "beangoben/mol_autoencoder" .
```

And then move to whatever folder you want to work with and execute:

```
docker run -p 8888:8888 -v "$(pwd)":/home/jovyan/work -it beangoben/mol_autoencoder
```
