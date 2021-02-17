FROM  continuumio/anaconda3

# Copy project files into container
# (not globbing all files to prevent IntelliJ IDE files from breaking Docker caching):
COPY conda_spec.txt /intSiteUploader/conda_spec.txt
COPY ./intSiteUploader.R /intSiteUploader/intSiteUploader.R
COPY ./installPackages.R /intSiteUploader/installPackages.R

WORKDIR /intSiteUploader

RUN apt-get update -y && \
    apt-get install -y libcurl4-openssl-dev \
    libmariadb-dev \
    gcc \
    make  && \
    apt-get clean

RUN ln -s /usr/lib/x86_64-linux-gnu/libmariadb.so /usr/lib/x86_64-linux-gnu/libmysqlclient.so.18

RUN conda create --name intSiteUploader --file conda_spec.txt

# Make RUN commands use the new environment:
# https://pythonspeed.com/articles/activate-conda-dockerfile/
SHELL ["conda", "run", "-n", "intSiteUploader", "/bin/bash", "-c"]

# Install non-Conda dependencies:
RUN Rscript installPackages.R
WORKDIR /scratch
ENTRYPOINT ["conda", "run", "-n", "intSiteUploader", "Rscript", "/intSiteUploader/intSiteUploader.R"]
