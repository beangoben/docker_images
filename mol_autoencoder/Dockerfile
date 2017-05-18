FROM  beangoben/pimp_jupyter

USER jovyan
ADD requirements.txt /home/jovyan/work
# chemio-informatics
conda install -c -n python2 rdkit rdkit
RUN source activate python2 && pip install --no-cache -r requirements.txt && rm requirements.txt
# change default enviroment to python2
RUN echo 'export PATH="/opt/conda/envs/python2/bin:$PATH"' >> ~/.bashrc


