## How to use

To run the software on any computer you need to install [docker](https://www.docker.com/).

Then you can either download or build the docker image.

To download running the following command in your favorite terminal:

```
docker pull beangoben/useful_chem
```

of build it (good to change things) by moving to the git cloned repository :

```
docker build -t "beangoben/useful_chem" .
```

And then move to whatever folder you want to work with and execute:

```
docker run -p 8888:8888 -v "$(pwd)":/home/jovyan/work -it "beangoben/useful_chem"
```