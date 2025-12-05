FROM continuumio/miniconda3

WORKDIR /app

# Update system libraries to ensure modern libstdc++
RUN apt-get update && apt-get install -y libstdc++6 \
    && rm -rf /var/lib/apt/lists/*

# Use conda-forge for scientific packages
RUN conda config --add channels conda-forge \
    && conda config --set channel_priority strict

# Install all required packages into the base environment
RUN conda install -y \
    python=3.11 \
    numpy \
    pandas \
    scipy \
    biopython \
    openmm \
    mdanalysis \
    rdkit \
    requests \
    pdbfixer \
    openbabel \
    libstdcxx-ng \
    && conda clean -afy

# Install pip-only packages (note the hyphenated name on PyPI)
RUN pip install useful-rdkit-utils
RUN pip install openbabel-wheel
RUN pip install molscrub


# Copy your project files
COPY . .

CMD ["bash"]