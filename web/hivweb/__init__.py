from flask import Flask, g
from flask.ext.cache import Cache
import logging

from . import filters

app = Flask(__name__, template_folder='templates/')
cache = Cache(app)

# Initiate engine before the first request
@app.before_first_request
def before_first_request():
    filters.register(app)	# register jinja filters in the app
    # Configure logging
    if not app.debug:
        from logging.handlers import SMTPHandler, RotatingFileHandler
        from logging import Formatter
        mail_handler = SMTPHandler('127.0.0.1',
                                   'cwarth@fredhutch.org',
                                   app.config['ADMINS'], 'hiv-sim website Failed')
        mail_handler.setLevel(logging.ERROR)
        log_path = app.config['LOG_PATH']
        file_handler = RotatingFileHandler(log_path, maxBytes=1 << 20, backupCount=5)
        file_handler.setFormatter(Formatter(
            '%(asctime)s %(levelname)s: %(message)s '
            '[in %(pathname)s:%(lineno)d]'
            ))
        app.logger.setLevel(logging.INFO)
        file_handler.setLevel(logging.INFO)
        app.logger.addHandler(file_handler)
        app.logger.addHandler(mail_handler)

from . import views
