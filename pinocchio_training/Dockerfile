FROM continuumio/anaconda3
RUN apt-get update ;\
    pip install --upgrade pip ;\
    pip install matplotlib ;\
    conda install pinocchio -c conda-forge 
RUN apt-get install curl wget;\
    wget https://github.com/stack-of-tasks/pinocchio/releases/download/v2.7.1/pinocchio-2.7.1.tar.gz ;\
    tar -xzf pinocchio-2.7.1.tar.gz;\
    pip install meshcat;
CMD [ "/bin/bash" ]