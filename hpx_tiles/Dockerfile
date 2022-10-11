FROM python:3.8-slim
WORKDIR /app

COPY requirements.txt .
RUN pip install -r requirements.txt

COPY hpx_tiling_logic.py .
ENTRYPOINT ["python", "hpx_tiling_logic.py"]
