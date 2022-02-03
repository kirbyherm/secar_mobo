FROM continuumio/anaconda3
#FROM python:3
#FROM mikebentley15/pagmo2

#FROM intel/oneapi-hpckit
MAINTAINER Kirby Hermansen <khermansen@protonmail.com>
RUN mkdir /cosy
COPY requirements.txt /cosy/

RUN pip install --no-cache-dir --upgrade pip &&\
    pip install --no-cache-dir -r /cosy/requirements.txt
#RUN conda config --add channels conda-forge &&\
#    conda config --set channel_priority strict &&\
#    conda install pygmo
RUN conda create -n pygmo_env -c conda-forge pygmo
SHELL ["conda","run","-n","pygmo_env","/bin/bash","-c"]
#RUN conda install conda-forge::mamba &&\
#    mamba install -c conda-forge pygmo
#CMD ["python3", "-m", "pip", "install", "numpy"]
#CMD ["python3", "-m", "pip", "install", "pandas"]
#RUN chown newuser /cosy
#USER newuser
ADD COSY10.0 /cosy/COSY10.0
ADD py /cosy/py
ADD fox /cosy/fox
RUN mkdir /cosy/output
ENTRYPOINT ["/cosy/COSY10.0/cosy","/cosy/fox/SECAR_an_Optics_DE.fox"]
##ENTRYPOINT ["cd", "/COSY10.0","&&","make","clean","&&","make"]
#WORKDIR /cosy/py
##RUN make clean
##RUN make
#ENTRYPOINT ["conda","run","-n","pygmo_env","python", "optimize.py", "1"]
