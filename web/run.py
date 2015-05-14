#!/usr/bin/env python

import os
import os.path
from hivweb import app as application
import logging

if __name__ == '__main__':
    # add files to be watched for changes
    # http://stackoverflow.com/a/9511655/1135316
    extra_files = [os.path.join(os.getcwd(), "util.py"),
                   os.path.join(os.getcwd(), "process.py"),
                   os.path.join(os.getcwd(), "filters.py")]
    application.config['DEBUG'] = True
    application.run("0.0.0.0", port=5000, debug=True, extra_files=extra_files)
