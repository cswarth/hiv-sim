#!/usr/bin/env make

# We need to build our own virtual env, including all of numpy and
# scipy, because the existing python environment under
# ~matsengrp/local is too old to be useful.  In particular there are
# routines used when calculating sums of log probabilities that are not
# available under the old verion of scipy installed in ~matsengrp.



# use the virtual environment
export PATH :=venv/bin:${PATH}

# build a python virtual environment
venv:
	virtualenv -p /usr/bin/python venv
	easy_install pip
	@#pip install -U distribute
	@# Unfortunately cannot install SCons from requirments.txt b/c
	@# installation requires the --egg option.
	@# http://stackoverflow.com/a/19697682/1135316
	pip install --egg -q SCons
	pip install -r python/requirements.txt



.PHONY: packages venv

