FROM  beangoben/pimp_jupyter

USER jovyan

# chemio-informatics
RUN conda install -c rdkit rdkit -q -y -f && \
    conda clean --all
# bayesopt
RUN pip install git+https://github.com/fmfn/BayesianOptimization.git




