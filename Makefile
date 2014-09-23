#!/usr/bin/env make

# We need to build our own virtual env, including all of numpy and
# scipy, because the existing python environment under
# ~matsengrp/local is too old to be seful.  In particular there are
# routines used when calculting sums of log probabilities that are not
# available under the old verion of scipy installed in ~matsengrp.



# use the virtual environment
export PATH :=venv/bin:${PATH}

# build a python virtual environment
venv:
	virtualenv -p /usr/bin/python venv
	easy_install pip
	#pip install -U distribute
	pip install --egg SCons
	pip install git+https://github.com/fhcrc/nestly.git
	pip install -r python/requirements.txt



.phony: packages

