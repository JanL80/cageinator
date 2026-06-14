FROM python:3.11-slim

ARG APP_VERSION=unknown
LABEL version="${APP_VERSION}"
LABEL description="Metal-organic cage assembler"

### Install necessary system dependencies
### - wget and xz-utils for downloading and extracting xTB
### - libgomp1 for xTB's OpenMP multi-threading requirements
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    xz-utils \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

### Download and install xTB (v6.7.1 stable)
RUN mkdir -p /opt/xtb && \
    wget https://github.com/grimme-lab/xtb/releases/download/v6.7.1/xtb-6.7.1-linux-x86_64.tar.xz -O /tmp/xtb.tar.xz && \
    tar -xf /tmp/xtb.tar.xz -C /opt/xtb --strip-components=1 && \
    rm /tmp/xtb.tar.xz
### Set xTB environment variables
ENV PATH="/opt/xtb/bin:${PATH}"
ENV XTBPATH="/opt/xtb/share/xtb"
ENV OMP_NUM_THREADS=1

### Set working directory
WORKDIR /app
### Copy the entire project into the container
COPY . /app/

### Install the Cageinator package
RUN pip install --no-cache-dir .

### Expose default web server port
EXPOSE 5001

### Set the console script as the primary executable
ENTRYPOINT ["cageinator"]

### Default to the web interface if no arguments are passed
CMD []
