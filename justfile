package-name := "venum"
env := "env"
python := env + "/bin/python"
pip := env + "/bin/pip"
global-python := "python3"
container-engine := "podman"
image-name := package-name + "-build"
container-build-dir := "./dist_container"

test path='./':
    {{ python }} -m pytest tests/{{path}}

setup:
    @echo "Setting up virtual environment"
    @{{ global-python }} -m venv {{ env }}
    @{{ pip }} install -e .
    @{{ pip }} install build

build:
    {{ python }} -m build --wheel

build-container:
    #!/bin/sh
    {{ container-engine }} build . -f Containerfile -t {{ image-name }}
    temp_container=$({{ container-engine }} create {{ image-name }})
    mkdir -p {{ container-build-dir }}
    {{ container-engine }} cp $temp_container:/workspace/dist {{ container-build-dir }}
    {{ container-engine }} rm -v $temp_container
    mv {{ container-build-dir }}/dist/* {{ container-build-dir }}
    rmdir {{ container-build-dir }}/dist

clean:
    rm -rf dist/
    rm -rf {{ container-build-dir }}/
    rm -rf build/
    rm -rf *.egg-info

format:
    @{{ python }} -m black {{ package-name }}/ tests/

lint:
    @{{ python }} -m flake8 {{ package-name }}/ tests/

repl:
    @{{ python }}

run path:
    @{{ python }} {{ path }}
