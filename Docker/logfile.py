#! /usr/bin/env python3
from time import gmtime, strftime
import logging

def startlog(name):
 # create a file handler
 logfn = name + '.' + strftime("%Y%m%d", gmtime()) + '.log'
 logger = run(logfn)
 return(logger)

def run(logfn):
 level = logging.INFO
 logger = logging.getLogger(__name__)
 logger.setLevel(level)
 handler = logging.FileHandler(logfn)
 handler.setLevel(level)

 # create a logging format
 formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
 handler.setFormatter(formatter)

 # add the handlers to the logger
 logger.addHandler(handler)
 logger.addHandler(logging.StreamHandler())
 return(logger)

