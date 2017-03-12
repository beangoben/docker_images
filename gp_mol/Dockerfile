FROM  beangoben/pimp_jupyter

USER root

# chemio-informatics
RUN conda remove --name python2 -y -q --all  && conda clean --all
#RUN conda install -q -y anaconda && conda update --all && conda clean --all
RUN conda install -c rdkit rdkit -q -y && \
    conda clean --all
# Good GP from sheffieldML group
RUN pip install GPy && pip install gpyopt
# bayesopt
RUN pip install git+https://github.com/fmfn/BayesianOptimization.git




