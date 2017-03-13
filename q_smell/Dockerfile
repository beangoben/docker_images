FROM  beangoben/pimp_jupyter

USER root
RUN conda remove --name python2 -y -q --all  && conda clean --all
# svg2png + chemioinformatics
RUN pip install pubchempy cirpy py3Dmol decorator chemspipy
# update numpy scipy
RUN conda install -q -y -c conda-forge networkx flask cairocffi && conda clean --all
RUN conda update --yes pillow && conda clean --all
# rdkit
RUN conda install -c rdkit rdkit --quiet --yes && conda clean --all

EXPOSE 8000

USER jovyan

