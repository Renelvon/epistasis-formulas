NAME=epystatic
PIP=pip3
PYTHON=python3
SETUP=setup.py

.PHONY: all build build_ext check install

all: build

bdist_wheel:
	$(PYTHON) $(SETUP) bdist_wheel

build:
	$(PYTHON) $(SETUP) build

check:

clean:
	git clean -xfd

dist:
	$(PYTHON) $(SETUP) sdist

distclean: clean

install: build
	$(PYTHON) $(SETUP) install --user

uninstall:
	$(PIP) uninstall -y $(NAME)
