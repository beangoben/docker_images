FROM  beangoben/pimp_jupyter

USER root
# update numpy scipy
RUN conda update --yes numpy scipy matplotlib seaborn && \
    conda clean --all
# chemio-informatics
#RUN pip install cairosvg
RUN conda install -c rdkit rdkit -q -y -f && \
    conda clean --all
RUN conda install -c openbabel openbabel --quiet --yes && \
    conda clean --all
RUN conda install -c richlewis scikit-chem --quiet --yes && \
    conda clean --all
RUN pip install imolecule cclib py3dmol cirpy pubchempy chemspipy


