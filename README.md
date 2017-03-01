# Repo for docker images

## How to use along with Docker

To run the software on any computer you need to install [docker](https://www.docker.com/).

Then you can build the image by moving to the git cloned repository and running the command:

```
docker build -t "name" .
```

And then move to whatever folder you want to work with and execute:

```
docker run -p 8888:8888 -v "$(pwd)":/home/jovyan/work -it "name"
```