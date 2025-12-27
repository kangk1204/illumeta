FROM mambaorg/micromamba:1.5.8

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba create -y -n illumeta -f /tmp/environment.yml && \
    micromamba clean -a -y

ENV PATH=/opt/conda/envs/illumeta/bin:/opt/conda/bin:$PATH
WORKDIR /app
COPY . /app

# Pre-install R/Bioconductor dependencies inside the container.
RUN ILLUMETA_FORCE_SETUP=1 Rscript r_scripts/setup_env.R

ENTRYPOINT ["python3", "illumeta.py"]
CMD ["--help"]
