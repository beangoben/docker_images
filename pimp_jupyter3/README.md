# pimp my jupyter

Docker container for a pimped-up jupyter notebook.

## What's included?

- based on the [jupyter/scipy-notebook](https://github.com/jupyter/docker-stacks/tree/master/scipy-notebook) container, has joyvan user with two anaconda enviroments: python 3.5 (default) and python 2.7.
- [jupyter themes](https://github.com/merqurio/jupyter_themes), to modify syntax themes and code font.
- [RISE slideshow](https://github.com/damianavila/RISE), to create jupyter "Live" Reveal.js slideshows.
- [jupyter nbextensions](https://github.com/ipython-contrib/jupyter_contrib_nbextensions), many useful extensions for notebooks. Many usability and styling extensions are activated by default.
- **custom.css, custom.js** in jupyter's custom folder. They change the headers, default font and many more css changes. Default line numbers.

## How to use

To run the software on any computer you need to install [docker](https://www.docker.com/).

Then you can either download or build the docker image.

To download running the following command in your favorite terminal:

```
docker pull beangoben/pimp_jupyter
```

of build it (good to change things) by moving to the git cloned repository :

```
docker build -t "beangoben/pimp_jupyter" .
```

And then move to whatever folder you want to work with and execute:

```
docker run -p 8888:8888 -v "$(pwd)":/home/jovyan/work -it beangoben/pimp_jupyter
```
