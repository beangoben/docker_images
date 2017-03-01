FROM  beangoben/pimp_jupyter

USER root
RUN apt-get update && \
    apt-get install -y --no-install-recommends imagemagick librsvg2-bin \
    libcairo2-dev libgd-dev libffi-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
USER jovyan
# update numpy scipy
RUN conda update --yes numpy scipy matplotlib seaborn
# chemio-informatics
RUN pip install cairosvg
RUN conda install -c rdkit rdkit --quiet --yes
RUN conda install -c openbabel openbabel --quiet --yes
RUN conda install -c richlewis scikit-chem --quiet --yes
RUN pip install imolecule cclib py3dmol cirpy pubchempy chemspipy


