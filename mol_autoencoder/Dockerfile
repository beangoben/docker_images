FROM  beangoben/pimp_jupyter

USER jovyan
ADD requirements.txt /home/jovyan/work
# chemio-informatics
conda install -n python2 -c rdkit rdkit
RUN source activate python2 && pip install --no-cache -r requirements.txt && rm requirements.txt
# change default enviroment to python2
RUN echo 'export PATH="/opt/conda/envs/python2/bin:$PATH"' >> ~/.bashrc


