#!/APSshare/epd/rh6-x86_64/bin/python

"""Monitor the ORNL Cluster server at APS/sector 34 status and sets an EPICS PV to the current statusGUI to run 
	modified from Jon Tischler's version (APS/XSD, Argonne National Laboratory)
"""
import xmlrpc.client
import os
import sys
import datetime
from .display34ide import *

HOSTNAME = 'hpcs34' # name of host with cluster server to monitor

SERVER_NAME = 'ORNL cluster server'
MIN_SERVER_VERSION = 3.0			# this prgram requires MIN_SERVER_VERSION of the server

class Server(object):

	def __init__(self,two=False):
		self.two = two
		if MIN_SERVER_VERSION>=3.0:
			self.portNum = 8885
			if two: self.portNum = 8883
		elif MIN_SERVER_VERSION>=2.0:	self.portNum = 8886
		else:							self.portNum = 8888
		self._setStateUnconnected()


	# CONNECT to the server
	def connect(self,hostname):
		"""Try to estableis a new connection to a server"""
		if len(hostname)<1: hostname = self.hostname
		self.hostname = hostname
		if len(hostname)<1:
			txt = 'FAILED to connect to server, no hostname'
			#print txt
			return (False, txt)
		self.proxy = xmlrpc.client.ServerProxy('http://'+hostname+':'+str(self.portNum),allow_none=True)

		self.getStatus()
		if not self.connected:
			txt = 'connection to "'+hostname+'"  FAILED'
			#print txt
			return (False,txt)
		txt = '	connected to "'+hostname+'" on port '+str(self.portNum)+' with , status = '+str(self.status)

		# check server version
		(name,version) = self._getVersion()
		if name != SERVER_NAME:
			self._setStateUnconnected()
			#print txt

			#print 'len(name) =',len(name),'     len(SERVER_NAME) =',len(SERVER_NAME)

			txt = 'FAILED, server name is "'+name+'", but this client needs "'+SERVER_NAME+'"'
			#print txt
			return (False,txt)

		elif version < MIN_SERVER_VERSION:
			self._setStateUnconnected()
			#print txt
			txt = 'FAILED, server version is',version,', but we need at least ',MIN_SERVER_VERSION
			#print txt
			return (False,txt)

		txt = 'connected to "'+hostname+':'+str(self.portNum)+'", status = '+str(self.status)+', server version = '+str(version)
		self.serverVersion = version
		#print txt
		return (True,txt)


	# GET current status info about the server
	def getStatus(self):
		"""Get the server status, should be one of 'unknown', 'idle', or 'running'"""
		if self.proxy is None:
			self._setStateUnconnected()
			return self.status

		status = ''
		try:	status = self.proxy.status()
		except:	self._setStateUnconnected()

		if status.find('status=')==0:
			self.status = status[7:]
			self.connected = True
		else: self.status = 'unknown'
		return self.status


	def _getVersion(self):
		"""Get the version of the server"""
		if self.proxy is None: return ('',-1)
		try:	value = self.proxy.version()
		except:
			return ('',-1)

		try:
			name = value[0]
			version = float(value[1])
		except:
			return ('',-1)
		return (name,version)


	def _setStateUnconnected(self):
		"""Set the state to un-connected, and all associated variables too."""
		self.status = 'unknown'
		self.connected = False
		self.hostname = ''
		self.proxy = None
		self.version = -1
        
	def fraction(self):
		"""Get the fraction of current computation that is finished from server"""
		badValue = {'fraction':-1.0, 'elapsedSec':-1, 'remainSec':-1, 'endAt':''}
		if self.proxy is None: return badValue
		try:	fraction = self.proxy.fraction()
		except:	fraction = badValue

		# This little dance is needed because you cannot pass a timedelta thru xmlrpc
		try:
			elapsedSec = fraction['elapsedSec']
			fraction['elapsedSec'] = datetime.timedelta(elapsedSec[0],elapsedSec[1],elapsedSec[2])
		except:	fraction['elapsedSec'] = None
		try:
			remainSec = fraction['remainSec']
			fraction['remainSec'] = datetime.timedelta(remainSec[0],remainSec[1],remainSec[2])
		except:	fraction['remainSec'] = None
		try:
			endAt = datetime.datetime.strptime(fraction['endAt'].value, "%Y%m%dT%H:%M:%S")
			fraction['endAt'] = endAt
		except:	fraction['endAt'] = None
		return fraction


def getClusterStatus():
    server1 = Server()
    server1.connect(HOSTNAME)
    server1msg = server1.getStatus()
    percent = server1.fraction()['fraction']*100
    server1pct = '%6.1f' % percent + '%'
    server2 = Server(two=True)
    server2.connect(HOSTNAME)
    server2msg = server2.getStatus()
    percent = server2.fraction()['fraction']*100
    server2pct = '%6.1f' % percent + '%'
    return [(server1msg,server1pct), (server2msg,server2pct)]

def clusterpage():
    """ 
    main program to build a dynamical html page based on real-time PV values.
    """
    
    #imports
    from datetime import datetime
    
    # start html body, write page title
    print("Content-type: text/html")
    print('\n')
    print("<html>")
    print(htmltag('34IDE Computer Cluster Status','title',{}))
    print("<body bgcolor=\"#D6D6D6\">")
    print(htmltag('Current Cluster Status','h3',{'style':"line-height:90%;"}))
    message = datetime.now().strftime("%I:%M:%S %p, %a %Y-%m-%d")
    message = htmlpg(('@ ' + message, smallfont+itfont+"line-height:50%;"))
    print(message)
    
    # Construct main msg body
    serverStatus = getClusterStatus()
    for id in range(2):
        server = serverStatus[id]
        if server[0] == 'running':
            serverStatus[id] = htmlstylespan(','.join(server),boldfont+itfont+bluecolor)
        elif server[0] == 'idle':
            serverStatus[id] = htmlstylespan('idle',boldfont+dgreencolor)
        else:
            serverStatus[id] = htmlstylespan(server[0],boldfont)
    msglist = ['Main XMD Cluster: ' + serverStatus[0]]
    msglist.append('Extra Cluster: ' + serverStatus[1])
    print(htmlulist(msglist))

    # provide links
    print(htmlpg(htmllink('Back to 34IDE beamline status','http://sector34.xray.aps.anl.gov/34ide/status.cgi')))
    print(htmlhr) # horiz. line
    print(extlinks())
    
    # finish html body & return
    print("</body>")
    print("</html>")
    return
    
if __name__ == '__main__':
    clusterpage()
