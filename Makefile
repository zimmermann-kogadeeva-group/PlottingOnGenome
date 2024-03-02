
SHELL = /bin/bash
ACTIVATE = source .venv/bin/activate

.venv:
	python3 -m venv $@

.PHONY: install_packages
install_packages: .venv
	${ACTIVATE} && pip install -e .

.PHONY: requirements.txt
requirements.txt:
	${ACTIVATE} && pip freeze > $@

.PHONY: check
check:
	flake8 src

.PHONY: clean
clean:
	rm -rf src/*.egg-info build dist **/__pycache__/ plotting_on_genome.spec

dist/plotting_on_genome:
	${ACTIVATE} && pyinstaller -F -w -n $(notdir $@) pyinstaller_entry.py

