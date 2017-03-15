FROM  beangoben/pimp_jupyter3

USER root

RUN conda install -c rdkit rdkit -q -y && \
    conda clean -tipsy
# Good GP from sheffieldML group
RUN pip install --no-cache GPy && pip install gpyopt
# bayesopt
RUN pip install --no-cache git+https://github.com/fmfn/BayesianOptimization.git




