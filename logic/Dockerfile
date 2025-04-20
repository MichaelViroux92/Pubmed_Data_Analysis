FROM jupyter/scipy-notebook:latest
WORKDIR /home/jovyan
COPY requirements.txt /home/jovyan/requirements.txt
RUN pip install --no-cache-dir  -r /home/jovyan/requirements.txt
COPY ./ /home/jovyan/
CMD ["uvicorn", "api_endpoints:app", "--host", "0.0.0.0", "--port", "8001"]
