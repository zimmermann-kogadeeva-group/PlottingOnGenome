# app/Dockerfile

FROM python:3.12-slim

WORKDIR /app

COPY . /app

# Installing deps
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    zip \
    software-properties-common

# Installing pog package
RUN pip install --no-cache ".[streamlit]"

# Adding BLAST
RUN curl https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.16.0/ncbi-blast-2.16.0+-x64-linux.tar.gz -o blast.tar.gz && \
    tar -C /opt/ -zxvf blast.tar.gz && \
    rm blast.tar.gz

# Updating path env var
ENV PATH="${PATH}:/opt/ncbi-blast-2.16.0+/bin"

# Cleaning up
RUN apt remove --purge curl zip && rm -rf /var/lib/apt/lists/*

EXPOSE 8501

HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health

ENTRYPOINT ["streamlit", "run", "streamlit_app.py", "--server.port=8501", "--server.address=0.0.0.0"]
