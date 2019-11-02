NAME=epystatic
PIP=pip3
PYTHON=python3
SETUP=setup.py

.PHONY: all bdist_wheel build check clean dist distclean install installcheck uninstall

all: build

bdist_wheel:
	$(PYTHON) $(SETUP) bdist_wheel

build:
	$(PYTHON) $(SETUP) build

check:
	black $(SETUP) $(NAME)
	check-manifest
	pylint $(SETUP) tests
	pyroma -n 10 .

clean:
	git clean -xfd

dist:
	$(PYTHON) $(SETUP) sdist

distclean: clean

install: build
	$(PYTHON) $(SETUP) install

installcheck:
	nose2 tests

uninstall:
	$(PIP) uninstall -y $(NAME)
