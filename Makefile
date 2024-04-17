
SHELL = /bin/bash
ACTIVATE = source .venv/bin/activate

.PHONY: venv build upload check isort tags clean data

.venv:
	python3 -m venv $@

venv: .venv
	${ACTIVATE} && pip install -e .

requirements.txt:
	${ACTIVATE} && pip freeze > $@

build:
	${ACTIVATE} && python3 -m build

upload: build
	${ACTIVATE} && python3 -m twine upload -r pypi -u __token__ dist/*

check:
	flake8 src

isort:
	isort src

tags:
	ctags-universal --recurse src

clean:
	rm -rf src/*.egg-info/ **/__pycache__/ build/ dist/ report.xml

dist/plotting_on_genome:
	${ACTIVATE} && pyinstaller -F -w -n $(notdir $@) pyinstaller_entry.py

