# zxdcfd 多人网页版: OpenFOAM + gmsh + Streamlit 装在一个镜像里，计算在容器内直接跑（本机方式）。
FROM opencfd/openfoam-default:2412
USER root
# OpenFOAM 已装在基础镜像里；它的 apt 源有时同步不全会让 apt-get update 失败，这里用不到，先去掉。
RUN grep -rl openfoam /etc/apt/sources.list.d 2>/dev/null | xargs -r rm -f \
    && apt-get update && apt-get install -y --no-install-recommends \
        python3 python3-pip python3-venv libglu1-mesa libgl1 libxcursor1 libxft2 libxinerama1 \
        libxrender1 libfontconfig1 \
    && rm -rf /var/lib/apt/lists/*
RUN python3 -m venv /opt/venv && /opt/venv/bin/pip install --no-cache-dir gmsh 'streamlit>=1.50' plotly pandas
ENV PATH=/opt/venv/bin:$PATH \
    ZXDCFD_RUN_MODE=本机 \
    ZXDCFD_WORKDIR=/data \
    PYTHONUNBUFFERED=1
WORKDIR /app
COPY *.py /app/
COPY .streamlit /app/.streamlit
COPY tests/smoke_docker.py /app/tests/smoke_docker.py
RUN useradd -m -s /bin/bash zxdcfd && mkdir -p /data && chown zxdcfd:zxdcfd /data
USER zxdcfd
VOLUME /data
EXPOSE 8501
ENTRYPOINT []
CMD ["streamlit", "run", "app.py", "--server.address=0.0.0.0", "--server.port=8501", "--server.headless=true"]
