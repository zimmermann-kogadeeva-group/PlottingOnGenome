# app/Dockerfile

FROM python:3.11-slim

WORKDIR /app

COPY . /app

RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    software-properties-common \
    git \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install . && pip3 install streamlit

RUN curl https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.15.0+-x64-linux.tar.gz -o blast.tar.gz && \
    tar -C /opt/ -zxvf blast.tar.gz && \
    rm blast.tar.gz

ENV PATH="${PATH}:/opt/ncbi-blast-2.15.0+/bin"

EXPOSE 8501

HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health

ENTRYPOINT ["streamlit", "run", "streamlit_app.py", "--server.port=8501", "--server.address=0.0.0.0"]
