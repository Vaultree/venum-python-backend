FROM python:3.10-slim

ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

WORKDIR /workspace

COPY . .

RUN $VIRTUAL_ENV/bin/pip install build \
	&& $VIRTUAL_ENV/bin/pip install -e .

RUN $VIRTUAL_ENV/bin/python -m pytest tests/ \
	&& $VIRTUAL_ENV/bin/python -m build --wheel
