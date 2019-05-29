# Use an official Python runtime as a parent image
FROM ubuntu

RUN apt-get update
RUN apt-get upgrade -y 
RUN apt-get install git -y
#RUN apt-get install gcc -y
#RUN apt-get install gsl-bin -y
#RUN apt-get install libgsl-dev -y
#RUN git clone https://github.com/compomics/searchgui.git
#RUN git clone --branch consistent_search_params https://github.com/marcottelab/MSblender.git

RUN git clone https://github.com/marcottelab/diffrac.git

#WORKDIR /MSblender/src/
#RUN /MSblender/src/compile

RUN apt-get install python -y
RUN apt-get install python-dev -y
#RUN apt-get install default-jre -y
RUN apt-get install wget -y


RUN apt-get install python-pip -y
#RUN apt-get install python-mysqldb -y
#kdrew: still interactive
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get install r-base -y

RUN pip install numpy==1.13.1
RUN pip install pandas==0.20.3 
RUN pip install scipy==0.19.1
RUN pip install scikit-learn==0.19.0
RUN pip install pyemd

RUN pip install statsmodels
RUN pip install numexpr==2.4.6

#RUN pip install networkx
#RUN pip install AgglomCluster==1.0.2
#RUN pip install matplotlib
#RUN pip install seaborn
#RUN pip install rpy2==2.8.6
#RUN pip install requests
#RUN pip install biopython
#RUN pip install mpmath
#RUN pip install numexpr

# Set the working directory to /app
WORKDIR /diffrac_work/

#COPY protein_complex_maps/ /protein_complex_maps/

ENV PYTHONPATH "${PYTHONPATH}:/diffrac/"

#ENTRYPOINT ["bash", "/MSblender/runMS2.sh"]
