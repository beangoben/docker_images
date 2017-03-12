FROM  beangoben/pimp_jupyter

USER jovyan

# chemio-informatics
RUN conda install -c rdkit rdkit -q -y -f && \
    conda clean --all
# Good GP from sheffieldML group
RUN pip install GPy && pip install gpyopt
# bayesopt
RUN pip install git+https://github.com/fmfn/BayesianOptimization.git




