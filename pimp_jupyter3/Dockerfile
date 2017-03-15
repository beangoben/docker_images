FROM  beangoben/pimp_jupyter

# make bash default shell
USER root

RUN conda remove --name python2 -y -q --all  &&\
    conda uninstall -y -q boto scikit-image sqlachemy vincent \
    beautifulsoup4 sympy numexpr bokeh &&\
    conda update --all -y -q &&\
    conda clean --all

