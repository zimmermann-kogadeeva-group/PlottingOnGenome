
SHELL = /bin/bash
ACTIVATE = source .venv/bin/activate

.venv:
	python3 -m venv $@

.PHONY: install_packages
install_packages: .venv requirements.txt
	${ACTIVATE} && pip install -r requirements.txt

requirements.txt:
	${ACTIVATE} && pip freeze > $@

.git/hooks/pre-commit:
	ln -s ../../.pre-commit $@

.PHONY: install_precommit
install_precommit: .git/hooks/pre-commit

clean:
	rm -rf .venv *.egg-info build dist

dist/plotting_on_genome:
	pyinstaller -F -w -n $(notdir $@) pyinstaller_entry.py
