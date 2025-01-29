FROM ubuntu:latest  

# Install dependencies (curl, wget, gnupg2, lsb-release)
RUN apt-get update && apt-get install -y curl wget gnupg2 lsb-release && \
    # Download and run the Entrez Direct installation script
    curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -o install-edirect.sh && \
    sh install-edirect.sh

WORKDIR /data  
CMD ["sleep", "infinity"]
