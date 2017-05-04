#!/usr/bin/env python

import argparse
import logging
import logging.handlers
from multiprocessing import Process
import os
import pickle
from SimpleXMLRPCServer import SimpleXMLRPCServer
from SimpleXMLRPCServer import SimpleXMLRPCRequestHandler
import socket
import SocketServer
import struct

from gips import utils
from gips.datahandler import geokit_api


class ThreadedRPC(SocketServer.ThreadingMixIn, SimpleXMLRPCServer):
    pass

class RequestHandler (SimpleXMLRPCRequestHandler):
    rpc_paths = ('/RPC2',)


class LogRecordStreamHandler (SocketServer.StreamRequestHandler):
    """Handler for a streaming logging request.

    This basically logs the record using whatever logging policy is
    configured locally.
    """

    def handle (self):
        """
        Handle multiple requests - each expected to be a 4-byte length,
        followed by the LogRecord in pickle format. Logs the record
        according to whatever policy is configured locally.
        """
        while True:
            chunk = self.connection.recv(4)
            if len(chunk) < 4:
                break
            slen = struct.unpack('>L', chunk)[0]
            chunk = self.connection.recv(slen)
            while len(chunk) < slen:
                chunk = chunk + self.connection.recv(slen - len(chunk))
            obj = self.unPickle(chunk)
            record = logging.makeLogRecord(obj)
            self.handleLogRecord(record)

    def unPickle (self, data):
        return pickle.loads(data)

    def handleLogRecord (self, record):
        # if a name is specified, we use the named logger rather than the one
        # implied by the record.
        if self.server.logname is not None:
            name = self.server.logname
        else:
            name = record.name
        logger = logging.getLogger(name)
        # N.B. EVERY record gets logged. This is because Logger.handle
        # is normally called AFTER logger-level filtering. If you want
        # to do filtering, do it at the client end to save wasting
        # cycles and network bandwidth!
        logger.handle(record)


class LogRecordSocketReceiver (SocketServer.ThreadingTCPServer):
    """
    Simple TCP socket-based logging receiver suitable for testing.
    """

    allow_reuse_address = 1

    def __init__ (self, host='localhost',
                 port=logging.handlers.DEFAULT_TCP_LOGGING_PORT,
                 handler=LogRecordStreamHandler):
        SocketServer.ThreadingTCPServer.__init__(self, (host, port), handler)
        self.abort = 0
        self.timeout = 1
        self.logname = None

    def serve_until_stopped (self):
        import select
        abort = 0
        while not abort:
            rd, wr, ex = select.select([self.socket.fileno()],
                                       [], [],
                                       self.timeout)
            if rd:
                self.handle_request()
            abort = self.abort


def serve_log (host, port):
    logging.basicConfig(format = (
        '---- %(asctime)s %(jobid)s %(caller)s -----\n'
        '%(message)s\n'
        '--------------------------------------------------------'
    ))
    tcpserver = LogRecordSocketReceiver(host, port)
    tcpserver.serve_until_stopped()


def serve_xmlrpc (host, port):
    #server = SimpleXMLRPCServer(
    server = ThreadedRPC(
        (host, port),
        requestHandler=RequestHandler,
        allow_none=True,
    )
    server.register_introspection_functions()

    server.register_function(geokit_api.get_datacatalog)
    server.register_function(geokit_api.submit_job)
    server.register_function(geokit_api.job_status)
    server.register_function(geokit_api.stats_request_results)
    server.register_function(geokit_api.stats_request_results_filter)

    server.serve_forever()

    
def main ():
    parser = argparse.ArgumentParser(prog=os.path.split(__file__)[1],
                                     description='Geokit daemon')

    parser.add_argument('--host', help='hostname or address',
                        default=utils.settings().GEOKIT_SERVER)
    parser.add_argument('--logport', help='port for log server',
                        default=utils.settings().LOG_PORT)
    parser.add_argument('-l', '--logfile', help='logfile',
                        default='geokitd.log')
    parser.add_argument('--xmlrpcport', help='port for xmlrpc server',
                        default=utils.settings().XMLRPC_PORT)

    args = parser.parse_args()

    p = Process(target=serve_log, args=(args.host, args.logport))
    p.start()

    p = Process(target=serve_xmlrpc, args=(args.host, args.xmlrpcport))
    p.start()

    
if __name__ == '__main__':
    main()
    

