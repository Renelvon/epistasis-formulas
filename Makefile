NAME=epystatic
PIP=pip3
PYTHON=python3
SETUP=setup.py

.PHONY: all bdist_wheel build check clean dist distclean install test uninstall

all: build

bdist_wheel:
	$(PYTHON) $(SETUP) bdist_wheel

build:
	$(PYTHON) $(SETUP) build

check: test
	black $(SETUP) $(NAME)
	check-manifest
	pylint $(SETUP) $(NAME) tests
	pyroma -n 10 .

clean:
	git clean -xfd

dist:
	$(PYTHON) $(SETUP) sdist

distclean: clean

install: build
	$(PYTHON) $(SETUP) install

test:
	$(PYTHON) -m unittest

uninstall:
	$(PIP) uninstall -y $(NAME)
