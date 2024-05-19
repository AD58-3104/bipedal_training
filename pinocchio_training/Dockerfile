FROM continuumio/anaconda3
RUN apt-get update
RUN pip install --upgrade pip && \
    pip install matplotlib
RUN conda install pinocchio -c conda-forge
CMD [ "/bin/bash" ]