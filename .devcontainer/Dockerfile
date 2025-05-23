FROM quay.io/jupyter/r-notebook:r-4.3.1

### Environment variables
ENV GITHUB_CLI_VERSION 2.30.0
ENV QUARTO_VERSION 1.4.549

###########################
### SYSTEM INSTALLATION ###
###########################
USER root

### System dependencies. Feel free to add packages as necessary.
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        # Basic system usage
        lmodern \
        file \
        curl \
        g++ \
        tmux \
        ###################################################
        ### Add your own system dependencies installed  ###
        ### with `apt-get` as needed below this comment ###
        ### Example (note the backslash after name):    ###
        ### neofetch \                                  ###
        ###################################################
        && \
    apt-get clean -y && \
    rm -rf /var/lib/apt/lists/* /tmp/library-scripts

### R compiler settings
RUN R -e "dotR <- file.path(Sys.getenv('HOME'), '.R'); if(!file.exists(dotR)){ dir.create(dotR) }; Makevars <- file.path(dotR, 'Makevars'); if (!file.exists(Makevars)){  file.create(Makevars) }; cat('\nCXXFLAGS=-O3 -fPIC -Wno-unused-variable -Wno-unused-function', 'CXX14 = g++ -std=c++1y -fPIC', 'CXX = g++', 'CXX11 = g++', file = Makevars, sep = '\n', append = TRUE)"
RUN chmod 666 ${HOME}/.R/Makevars

### CRAN mirror
RUN R -e "dotRprofile <- file.path(Sys.getenv('HOME'), '.Rprofile'); if(!file.exists(dotRprofile)){ file.create(dotRprofile) }; cat('local({r <- getOption(\"repos\")', 'r[\"CRAN\"] <- \"https://cloud.r-project.org\"', 'options(repos=r)', '})', file = dotRprofile, sep = '\n', append = TRUE)"
RUN chmod 666 ${HOME}/.Rprofile

### Quarto
# versions: https://quarto.org/docs/download/_download.json
# neat setup: https://github.com/jeremiahpslewis/reproducibility-with-quarto
RUN curl --silent -L --fail \
        https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb > /tmp/quarto.deb && \
    apt-get update && \
    apt-get install -y --no-install-recommends /tmp/quarto.deb && \
    rm -rf /tmp/quarto.deb /var/lib/apt/lists/* /tmp/library-script && \
    apt-get clean

#########################
### USER INSTALLATION ###
#########################
USER ${NB_USER}

### Anaconda (conda/mamba) packages
RUN mamba install -y -c conda-forge --freeze-installed \
        # Jupyter setup
        jupyter-server-proxy=4.1.0 \
        #######################################################
        ### Add your own conda dependencies installed with  ###
        ### `conda/mamba` as needed below this comment      ###
        ### Example (note the backslash after name):        ###
        ### scikit-learn \                                  ###
        #######################################################
        && \
    mamba clean --all

### PyPI (pip) packages
RUN pip install \
        nbgitpuller \
        jupyterlab-quarto==0.2.8 \
        radian==0.6.11 \
        ################################################
        ### Add your own PyPI dependencies installed ###
        ### with `pip` as needed below this comment  ###
        ### Example (note the backslash after name): ###
        ### scikit-ntk \                             ###
        ################################################
        && \
    jupyter labextension enable nbgitpuller

### R packages
    # Versioned
RUN R -q -e 'remotes::install_version("markdown", version="1.12", repos="cloud.r-project.org")' && \
    R -q -e 'remotes::install_version("languageserver", version="0.3.16", repos="cloud.r-project.org")' && \
    R -q -e 'remotes::install_version("httpgd", version="2.0.1", repos="cloud.r-project.org")' && \
    # Latest Dev Versions
    R -q -e 'remotes::install_github("ManuelHentschel/vscDebugger")' && \
    ##########################################################
    ### Add your own R dependencies installed as needed    ###
    ### below this comment but before `echo`.              ###
    ### Example (note the `&& \` after the command):       ###
    ### R -q -e 'install.packages("dplyr")' && \           ###
    ##########################################################
    echo

### GitHub CLI Installation
RUN wget https://github.com/cli/cli/releases/download/v${GITHUB_CLI_VERSION}/gh_${GITHUB_CLI_VERSION}_linux_amd64.tar.gz -O - | \
    tar xvzf - -C /opt/conda/bin gh_${GITHUB_CLI_VERSION}_linux_amd64/bin/gh --strip-components=2

### Prints Jupyter server token when terminal is opened
RUN echo "echo \"Jupyter server token: \$(jupyter server list 2>&1 | grep -oP '(?<=token=)[[:alnum:]]*')\"" > ${HOME}/.get-jupyter-url.sh && \
    echo "sh \${HOME}/.get-jupyter-url.sh" >> ${HOME}/.bashrc

RUN R -q -e 'install.packages(c("tidyverse", "patchwork", "parallel"))' && \
RUN R -q -e 'install.packages(c("rrpack", "MRCE", "CVglasso", "MultiVarSel", "glmnet"))' && \
    R -q -e 'remotes::install_github("ajmolstad/MSRL")' && \
    R -q -e 'remotes::install_github("roobnloo/dpglasso-fixed")' 
