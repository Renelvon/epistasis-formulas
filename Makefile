NAME=epystatic
PIP=pip3
PYTHON=python3
SETUP=setup.py
TESTS=tests

.PHONY: all bdist_wheel build check clean dist distclean install lint test uninstall

all: build

bdist_wheel:
	$(PYTHON) $(SETUP) bdist_wheel

build: bdist_wheel

check: lint test

clean:
	git clean -xfd

dist:
	$(PYTHON) $(SETUP) sdist

distclean: clean

install:
	$(PIP) install --user .

lint:
	black $(SETUP) $(NAME)
	check-manifest
	pylint $(SETUP) $(NAME) $(TESTS)
	pyroma -n 10 .

test:
	$(PYTHON) -m unittest

uninstall:
	$(PIP) uninstall -y $(NAME)
