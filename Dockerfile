FROM mambaorg/micromamba:1.5.8

ENV MAMBA_ROOT_PREFIX=/opt/conda
ARG MAMBA_DOCKERFILE_ACTIVATE=1

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba create -y -n illumeta -f /tmp/environment.yml && \
    micromamba clean -a -y

ENV PATH=/opt/conda/envs/illumeta/bin:/opt/conda/bin:$PATH \
    CONDA_PREFIX=/opt/conda/envs/illumeta \
    R_LIBS_USER=/opt/conda/envs/illumeta/illumeta-r-lib \
    ILLUMETA_RESPECT_R_LIBS_USER=1
WORKDIR /app
COPY . /app

# Pre-install R/Bioconductor dependencies outside /app so a user bind-mount of
# project files does not hide the installed R package library at runtime.
RUN mkdir -p "$R_LIBS_USER" && \
    ILLUMETA_FORCE_SETUP=1 micromamba run -n illumeta Rscript r_scripts/setup_env.R

ENTRYPOINT ["micromamba", "run", "-n", "illumeta", "python3", "illumeta.py"]
CMD ["--help"]
