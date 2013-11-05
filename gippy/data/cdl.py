#!/usr/bin/env python

from gippy.data.core import Data, main as datamain

class CDLData(Data):
	""" A CDL (Crop Data Layer) object """
	sensors = {'cdl': 'CDL'}
	rootdir = '/titan/data/CDL'
	_tiles_vector = 'states'
	_tiles_attribute = 'abbr'

def main(): datamain(CDLData)
