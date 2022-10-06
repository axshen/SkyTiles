FROM python:3.8-slim
WORKDIR /app

RUN apt-get update && apt-get install -y procps

COPY requirements.txt .
RUN pip install -r requirements.txt

COPY pyCASATILE.py .
ENTRYPOINT ['python', 'pyCASATILE.py']