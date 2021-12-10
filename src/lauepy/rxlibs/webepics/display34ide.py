#!/APSshare/epd/rh6-x86_64/bin/python
#

"""
cgi script for 34ide real-time status page

- Created Mon Apr 23 16:55:35 2012 by XU Ruqing
"""

# imports
import os
from epics import caget, ca
from .RXhtmlgen import htmlulist, htmlpg, htmllink, htmlstylespan, htmltag
from . import PVdef34ide as pvs

# string variables for font style specifications
greencolor = 'color:green;'
dgreencolor = 'color:darkgreen;'
redcolor = 'color:red;'
dorangecolor = 'color:darkorange;'
magentacolor = 'color:magenta;'
bluecolor = 'color:blue;'
boldfont = 'font-weight:bold;'
itfont = 'font-style:italic;'
ultext = 'text-decoration:underline;'
bigfont = 'font-size:large;'
smallfont = 'font-size:small;'
couriernewfont = "font-family:'Courier New', Courier, monospace;"

htmlbr = '<br />'
htmlhr = '<hr />'

def _caget(pv,valueType = None):
    """ do caget without giving error messages
    """
    from epics import caget
    try:
        if valueType is None:
            result = caget(pv, timeout = 1.5)
        elif type(valueType) == 'str':  #<<< TODO <<<
            pass
    except:
        result = None
    return result

    
def _caput(pv,value,wait = True):
    """ do caput without giving error messages
    """
    from epics import caput
    try:
        a = caput(pv, value, wait, timeout = 1.0)
        return 0
    except:
        return 1
    

def ca_test(pv):
    """ test if a PV is accessible, return value is boolean
    """
    temp = ca.initialize_libca()
    connect_ok = ca.connect_channel(ca.create_channel(pv))
    ca.finalize_libca()
    return connect_ok
    
def connectfailmsg(server):
    if server == pvs.prefix34idEPS :
        servername = '34ID EPS'
    else:
        servername = server
    message = 'Error connecting to the IOC of ' + servername + '!'
    return htmlpg((message,redcolor))
    
def pvmsg():
    """ 
    function that grabs PV values one by one and generate the main html text
    """
    html = ''
    ### scan part ###
    text = ''
    # test scan: busy? paused?
    scanbusy = 0.0
    for i in range(4):
        busy = _caget(pvs.ide34scanbusys[i])
        if busy is not None:
            scanbusy += busy 
    # if not busy, print last scan finished time, percentage, total time
    if scanbusy < 0.1:
        text += htmlpg('- Scan is idle')
        # make list of scan message
        msglist = []
        msg = "%5.1f" % _caget(pvs.ide34scanPercent)
        msglist.append("Last scan completion: " + msg +'%')
        msg = _caget(pvs.ide34scanEndTimeStr)
        if not msg:  # if empty string
            msglist.append("Last scan finished time: unknown")
        else:        # if not empty string
            msglist.append("Last scan finished at "+msg)
        # finish list
        text += htmlulist(msglist)
    
    # else, print more detailed scan info 
    else:
        # print scan running + pause_info
        paused = _caget(pvs.ide34scanPaused)
        if paused < 0.5:
            text += htmlpg(("- Scan is running...",dgreencolor+itfont+boldfont))
        else:
            text += htmlpg(("- Scan is paused!",dorangecolor+itfont+boldfont))
        # scan progress infos: percent, end time, remain time, paused time & times
        # start list
        msglist = []
        msg = "%6.2f" % _caget(pvs.ide34scanPercent)
        if paused <0.5:
            msglist.append("Scan progress: "+ htmlstylespan(msg + '%',dgreencolor+boldfont))
        else:
            msglist.append("Scan progress: "+ htmlstylespan(msg + '%',dorangecolor+boldfont))
        msg = _caget(pvs.ide34scanStartTimeStr)
        msglist.append('Scan started: ' + msg)
        msg = _caget(pvs.ide34scanRemainingTimeStr)
        msglist.append("Estimated remaining time: "+msg)
        msg = _caget(pvs.ide34scanEndTimeStr)
        msglist.append("Estimated end time: "+msg)
        msg = _caget(pvs.ide34scanElapsedTimeStr)
        msglist.append("Total elapsed time: "+msg)
        msg = _caget(pvs.ide34scanPauseTimeStr)
        msglist.append("Paused time: "+msg)
        # is there a number of paused times???
        # finish list
        text += htmlulist(msglist)
        
        # scan dimension + points in each dimension
        # start new list
        msglist = []
        scandim = int(round(_caget(pvs.ide34scanDimension)))
        msg = "%d" % scandim
        msglist.append("Scan dimensions: "+msg)
        msg = 'Number of points (outer loop first): '
        for i in range(scandim):
            if i>0 :
                msg += ' x '
            msg += "%d" % round(_caget(pvs.ide34scanPtsNos[scandim-i-1]))
        msglist.append(msg)
        
        # scan1 mode (fly or step)
        scan1mod = _caget(pvs.ide34scan1Mode)
        if scan1mod == 2 :
            stat = 'Fly Scan'; style = itfont + magentacolor
        elif scan1mod in [0,1]:
            stat = 'Step Scan'; style = ''
        else:
            stat = 'unknown (connection timeout)'; style = ''
        msg = 'Scan mode: ' + htmlstylespan(stat,style)
        msglist.append(msg)
        
        if scan1mod in [0,1]:
            # for step scans, display scan1 det triggers
            detlist = []
            for i in range(4):
                triggerpv = pvs.ide34scan1TriggerPVs[i]
                trigger = _caget(triggerpv)
                if trigger: # if not empty string
                    scalerfound = trigger.find(pvs.scalerTrigger)
                    PEfound = trigger.find(pvs.PEcommonStr)
                    # write expclicitly if one of PE detector
                    if PEfound >= 0:
                        PEid = int(trigger[7])  # '1' or '2' or '3'
                        PEname = 'PE detector #' + trigger[7]
                        # grab image number from PE ioc
                        # test connection
                        connect_ok = ca_test(pvs.PEfileNs[PEid-1])
                        if connect_ok :
                            PEfileNStr = '%d' % int(_caget(pvs.PEfileNs[PEid-1]))
                            detlist.append(PEname + ' (next file No. = ' + PEfileNStr + ')') 
                        else :
                            detlist.append(PEname + htmlstylespan(' (Error connecting to detector IOC)', dorangecolor))
                    # write explicitly if scaler
                    elif scalerfound >= 0:
                        detlist.append('Scalers')
                    # echo PV if unrecognized
                    else:
                        detlist.append(htmlstylespan(trigger,couriernewfont))

            msglist.append("Detectors triggered: " + str(len(detlist)))
            # give a message if detector list is empty
            if not detlist:
                detlist = [('NO detectors are triggered by scan1!', dorangecolor)]
            msglist.append(detlist)
        elif scan1mod == 2:
            # do nothing for fly scan.
            pass
        text += htmlulist(msglist)
    # end of 'scanbusy' condition
    
    
    # finish scan part
    html += htmltag('Scan Status: ', 'h3', {}, True) + text
    
    ### beam part ###
    text = ''
    # beam bad?
    beambad = _caget(pvs.ide34beambad)
    if beambad < 0.5: # good
        text += htmlpg(('- Beam is OK',dgreencolor + itfont))
    else:             # really bad
        text += htmlpg(('- Bad or no beam',redcolor))
    
    msglist = []
    
    # D-shutter
    DShutter = _caget(pvs.DShutterClosed)
    if DShutter > 0.5:
        stat = 'closed'; style = redcolor
    else:
        stat = 'open'; style = dgreencolor
    msg = 'D-shutter '+htmlstylespan(stat,style)
    msglist.append(msg)
    
    # FE shutter
    FEShutter = _caget(pvs.FEShutterClosed)
    if FEShutter > 0.5:
        stat = 'closed'; style = redcolor
    else:
        stat = 'open'; style = dgreencolor
    msg = 'Front End shutter ' + htmlstylespan(stat,style)
    msglist.append(msg)
    
    # mono/white?
    monoMode = _caget(pvs.ide34monoMode)
    if monoMode == 1:
        stat = 'white slitted'; style = boldfont
    elif monoMode ==2:
        stat = 'monochromatic'; style = itfont
    else:
        stat = 'pass-thru'; style = ''
    msg = 'Beam mode: ' + htmlstylespan(stat, style)
    msglist.append(msg)
    
    # mono-mode-specific info sublist, only displayed on mono mode
    if monoMode ==2:
        monolist = []
        # mono Energy
        monoE = _caget(pvs.ide34monoEnergy)
        msg = 'Beam energy = %.3f keV' % monoE
        monolist.append(msg)
        # ID synced with monoenergy?
        idsync = _caget(pvs.ide34monoIDSync)
        if idsync == 1 :
            stat = 'On'; style = ''
        else:
            stat = 'Off'; style = dorangecolor
        msg = 'Undulator sync: ' + htmlstylespan(stat, style)
        monolist.append(msg)
        msglist.append(monolist)
        
    # Undulator info sublist
    msglist.append('Undulator info:')
    idlist = []
    # Undulator control mode
    idmode = _caget(pvs.ID34dsControlMode)
    if int(idmode) == 0 :
        stat = 'User'; style = ''
    else:
        stat = 'Control Room'; style = dorangecolor
    msg = 'Control mode: ' + htmlstylespan(stat, style)
    idlist.append(msg)
    
    # id gap, taper
    idgap = _caget(pvs.ID34dsGapmm)
    msg = 'Gap = %.2f mm' % idgap
    idlist.append(msg)
    idtaper = _caget(pvs.ID34dsTapermm)
    msg = 'Taper = %.2f mm' % idtaper #<< conditional colorize?
    idlist.append(msg)
    # ID energy
    idenergy = _caget(pvs.ID34dsEnergy)
    idharm = _caget(pvs.ID34dsHarmonic)
    idenergy = idenergy/idharm
    msg = 'Energy = %.3f / %.3f  keV (1st / 3rd harmonic)' % (idenergy, idenergy*3)
    idlist.append(msg)
    # finish undulator info sublist
    msglist.append(idlist)
    
    # storage ring current
    srcurrent = _caget(pvs.ringCurrent)
    msg = 'Storage ring current = %.2f mA' % srcurrent
    msglist.append(msg)
    
    text += htmlulist(msglist)
    
    ## EPS status ##
    if beambad > 0.5: # do only if beam is bad
    
        # test if eps ioc is accessible #
        connet_eps_ok = ca_test(pvs.status34idEPS)
        if not connet_eps_ok:
            msg = '-- Beamline EPS unavailable, PLC or soft IOC could be down!'
            text += htmlpg((msg,redcolor))
        else:
            eps_idd = _caget(pvs.eps34iddOK)
            eps_ida = _caget(pvs.eps34idaOK)
            eps_vac = _caget(pvs.eps34idVacOK)
            if eps_idd*eps_ida*eps_vac:
                text += htmlpg('- Beamline EPS ' + htmlstylespan('OK', dgreencolor) + '.')
            else:
                text += htmlpg(('- Beamline EPS Fault:',redcolor))
                epslist = []
                if not eps_vac:
                    epslist.append(('Beamline vacuum bad',redcolor))
                if not eps_ida:
                    epslist.append(('Front End Shutter disabled',redcolor))
                if not eps_idd:
                    epslist.append(('D Shutter disabled',redcolor))
                    eps_ideHe = _caget(pvs.eps34ideHeOK)
                    if eps_ideHe:
                        msg = 'E hutch He pressure OK'
                    else:
                        msg = ('E hutch He pressure low',redcolor)
                    epslist.append([msg])
                text += htmlulist(epslist)
        
    # finish beam part
    html += htmltag('Beam Status: ', 'h3', {}, True) + text


    ### user custom part ### <<<<<< unfinished <<<<<<<
    # send PROC command to the EPICS Sequence
    msg = _caput(pvs.userPVrefresh,1,True)
    # see if there are valid PVs to display
    has_good_pvs = False
    for i in range(3):
        #pvenabled = _caget(pvs.userPVenables[i])
        pvvalid = _caget(pvs.userPVvalids[i])
        pvgood = pvvalid > 0.5 and pvvalid < 2.5 #pvenabled > 0.5 and  pvvalid < 0.5
        #if pvgood:
        #    msg = _caget(pvs.userPVstrings[1])
        has_good_pvs = pvgood or has_good_pvs
    
    if has_good_pvs :
        text = htmlpg(('Additional PVs:',boldfont))
        # read pv string & value
        msglist = []
        for i in range(3):
            #pvenabled = _caget(pvs.userPVenables[i])
            pvvalid = _caget(pvs.userPVvalids[i])
            if pvvalid > 0.5 and  pvvalid < 2.5 :
                pvstring = _caget(pvs.userPVstrings[i]).split(' ',1)[0]
                pvvalue = _caget(pvs.userPVvalues[i])
                msg = htmlstylespan(pvstring, couriernewfont) + '  :  ' + htmlstylespan(str(pvvalue), boldfont)
                msglist.append(msg)
        text += htmlulist(msglist)
    else:
        text = ''
    # end if has_good_pvs
    html += text
    
    return html
    
def extlinks():
    """ print some external links
    """
    html = htmlpg(('External Links',boldfont))
    linklist = []
    # 34ide home page
    linklist.append(htmllink('34IDE Home Page','http://www.aps.anl.gov/Sectors/33_34/microdiff/index.html'))
    # storage ring history
    linklist.append(htmllink('Storage ring history','http://www.aps.anl.gov/Facility/Status/'))
    # 34ide live cameras (internal access)
    linklist.append(htmllink('34IDE live cameras (internal access)','http://s34video.xray.aps.anl.gov'))
    # cluster status
    linklist.append(htmllink('HPC cluster status','http://shepherd.xray.aps.anl.gov/'))
    # ssm group page
    linklist.append(htmllink('SSM Group Home Page','http://www.aps.anl.gov/Sectors/33_34/'))
    
    html += htmlulist(linklist)
    
    return html

def main():
    """ 
    main program to build a dynamical html page based on real-time PV values.
    """
    
    #imports
    from datetime import datetime
    
    # start html body, write page title
    print("Content-type: text/html")
    print('\n')
    print("<html>")
    print(htmltag('34IDE Status','title',{}))
    print("<body bgcolor=\"#D6D6D6\">")
    print(htmltag('Current Status of 34ID-E','h2',{'style':"line-height:90%;"}))
    start_time = datetime.now()
    message = start_time.strftime("%I:%M:%S %p, %a %Y-%m-%d")
    message = htmlpg(('@ ' + message, smallfont+itfont+"line-height:50%;"))
    print(message)
    print(htmlhr) # add a horizontal line
    
    # setenv, try connect to 34ide
    os.environ["EPICS_CA_ADDR_LIST"]=' '.join([pvs.ip34ide,pvs.ip34ida, \
                         pvs.ip34idEPS, pvs.ipPE1, pvs.ipPE2, pvs.ipPE3])

    connect_ok = ca_test(pvs.status34ide)
        
    # on connection error give error message
    if not connect_ok:
        print(connectfailmsg(pvs.prefix34ide))
        print(htmlpg("-- Either the network or the VME crate is down."))
    # else construct main msg body
    else:
        print(pvmsg())

    # add cluster status link
    print(htmlhr) # horizontal line
    print(htmlpg(htmllink('Check Cluster Queue Status','http://sector34.xray.aps.anl.gov/34ide/serverstatus.cgi')))

    
    # finally: provide useful links
    print(htmlhr) # horiz. line
    print(extlinks())
    
    # print load time
    end_time = datetime.now()
    duration = (end_time-start_time).total_seconds()
    print(htmlpg(('Execution time: '+str(duration)+' seconds.', smallfont+itfont+"line-height:50%;")))

    # finish html body & return
    print("</body>")
    print("</html>")
    return
    
if __name__ == '__main__':
    #msglist = [('item1','color:grey'),'item2',(['item31','item32'],'color:yellow')]
    #message = htmlulist(msglist)
    #print message
    #message = 'a paragraph'
    #print htmlpg((message,'color:orange'))
    #main()
    pass
    
