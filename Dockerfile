FROM nfcore/base:1.14
LABEL authors="Barry Digby" \
      description="Docker container containing fastqc"

WORKDIR ./
COPY test.yml ./
RUN conda env create -f test.yml && conda clean -a
ENV PATH /opt/conda/envs/test_env/bin:$PATH
