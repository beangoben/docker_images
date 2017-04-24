FROM  beangoben/pimp_jupyter3

USER root
#ADD requirements.txt /home/jovyan/work
RUN pip install --no-cache cairosvg 
RUN conda install -c rdkit -q -y rdkit && \
    conda install -c conda-forge tensorflow keras && \
    conda clean --all

