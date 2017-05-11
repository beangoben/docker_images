FROM  beangoben/pimp_jupyter3

USER root
#ADD requirements.txt /home/jovyan/work
RUN pip install --no-cache cairosvg 
RUN conda install -c rdkit -q -y rdkit && \
    conda install pandas=0.19 && \
    conda install -c conda-forge tensorflow && \
    conda clean --all
RUN pip install keras
