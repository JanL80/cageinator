#!/bin/bash

set -e

IMAGE_NAME="cageinator"
IMAGE_TAG=$(grep -E '^__version__\s*=' cageinator.py | cut -d '"' -f 2)

echo "Building Docker image: ${IMAGE_NAME}:${IMAGE_TAG}..."

docker build -t ${IMAGE_NAME}:${IMAGE_TAG} .

echo ""
echo "======================================="
echo "            Build Complete!            "
echo "======================================="
echo ""
echo ""
echo ">>> TO RUN THE WEB UI:"
echo ""
echo "    docker run -d -p 5001:5001 --name cageinator_web ${IMAGE_NAME}:${IMAGE_TAG} --web"
echo "    (Access via http://localhost:5001 or your server's IP address)"
echo ""
echo ""
echo ">>> TO RUN THE CLI:"
echo ""
echo "    To process local files, a volume has to be mounted to the container."
echo "    Example:"
echo "    docker run --rm -v /path/to/local/data:/data ${IMAGE_NAME}:${IMAGE_TAG} --nodes /data/nodes --linkers /data/linkers --out /data/output"
echo ""
