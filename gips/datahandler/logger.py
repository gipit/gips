import logging, logging.handlers
import os
import sys

from gips import utils
from . import queue


class Logger(object):
    __instance = None
    def __new__ (cls):
        if Logger.__instance is None:
            Logger.__instance = object.__new__(cls)
            try:
                port = utils.settings().LOG_PORT
            except:
                port = logging.handlers.DEFAULT_TCP_LOGGING_PORT

            try:
                server = utils.settings().GEOKIT_SERVER
            except:
                server = 'localhost'

            try:
                level = utils.settings().LOG_LEVEL
            except:
                level = gippy.Options.Verbose()

            rootLogger = logging.getLogger('')
            rootLogger.setLevel(cls.log_level(level))
            socketHandler = logging.handlers.SocketHandler(server, port)
            rootLogger.addHandler(socketHandler)
        return Logger.__instance

    def __init__(self):
        self.extra = {}
        try:
            self.extra['jobid'] = queue.get_job_name()
        except queue.NoCurrentJobError:
            self.extra['jobid'] = os.path.split(sys.argv[0])[1]


    def log (self, message, level=1):
        utils.verbose_out(message, level)
        extra = self.extra.copy()
        extra['caller'] = sys._getframe(1).f_code.co_name
        logging.log(self.log_level(level), message, extra=extra)


    @staticmethod
    def log_level (loglevel=None):
        if loglevel is None:
            loglevel = gippy.Options.Verbose()
        # logging priority is reversed and 10x gippy verbosity
        level = {
            1: logging.CRITICAL,
            2: logging.ERROR,
            3: logging.WARNING,
            4: logging.INFO,
            5: logging.DEBUG,
        }
        if loglevel < 1:
            loglevel = 1
        if loglevel > 5:
            loglevel = 5
        return level.get(loglevel)

