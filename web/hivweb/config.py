import atexit
import os.path
import shutil
import tempfile

# Defaults
LOG_PATH = '/extra/www/flask-webapps/sidbweb/error.log'
CONTENT_DIR = os.path.join(os.path.dirname(__file__), 'static', 'content')

DEBUG = True
# Email recipients
ADMINS = ['cwarth@fredhutch.org']

# Cache settings
CACHE_TYPE = 'filesystem'
CACHE_DIR = tempfile.mkdtemp(prefix='hivsim')
#CACHE_OPTIONS = {'cache_dir': CACHE_DIR}

def clear_cache(d):
    shutil.rmtree(d)
atexit.register(clear_cache, CACHE_DIR)
