#!/APSshare/epd/SunOS_5.10-x86/bin/python
#

"""
cgi script for 34ide real-time status page
- Created Mon Apr 23 16:55:35 2012 by XU Ruqing
"""

# imports
import os
from epics import caget, ca
from RXhtmlgen import htmlulist, htmlpg, htmllink, htmlstylespan, htmltag
import PVlistdef as pvs

# define some strings for font style specifications
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

htmlbr = '<br />'

def ca_test(pv):
    temp = ca.initialize_libca()
    connect_ok = ca.connect_channel(ca.create_channel(pv))
    ca.finalize_libca()
    return connect_ok
    
def connectfailmsg(server):
    if server == pvs.prefix34idEPS:
        server = '34ID EPS'
    message = 'Error connecting to the IOC of ' + server + '!'
    return htmlpg((message,redcolor))
    
def pvmsg():
    html = ''
    ### scan part ###
    text = ''
    # test scan: busy? paused?
    scanbusy = 0.0
    for i in range(4):
        scanbusy += caget(pvs.ide34scanbusys[i])
    # if not busy, print last scan finished time, percentage, total time
    if scanbusy < 0.1:
        text += htmlpg('Scan is idle.')
        # make list of scan message
        msglist = []
        msg = "%5.1f" % caget(pvs.ide34scanPercent)
        msglist.append("Last scan done " + msg +'%')
        msg = caget(pvs.ide34scanEndTimeStr)
        msglist.append("Last scan finished at "+msg)
        # finish list
        text += htmlulist(msglist)
    
    # else, print more detailed scan info 
    else:
        # print scan running + pause_info
        paused = caget(pvs.ide34scanPaused)
        if paused < 0.5:
            text += htmlpg(("Scan is running...",dgreencolor+itfont+boldfont))
        else:
            text += htmlpg(("Scan is paused!",dorangecolor+itfont+boldfont))
        # scan progress infos: percent, end time, remain time, paused time & times
        # start list
        msglist = []
        msg = "%6.2f" % caget(pvs.ide34scanPercent)
        msglist.append("Scan progress: "+ msg + '%')
        msg = caget(pvs.ide34scanRemainingTimeStr)
        msglist.append("Estimated remaining time: "+msg)
        msg = caget(pvs.ide34scanEndTimeStr)
        msglist.append("Estimated end time: "+msg)
        msg = caget(pvs.ide34scanElapsedTimeStr)
        msglist.append("Total elapsed time: "+msg)
        msg = caget(pvs.ide34scanPauseTimeStr)
        msglist.append("Paused time: "+msg)
        # is there a number of paused times???
        # finish list
        text += htmlulist(msglist)
        
        # scan dimension + points in each dimension
        # start new list
        #text += htmlpg("Scan Details:")
        msglist = []
        scandim = int(round(caget(pvs.ide34scanDimension)))
        msg = "%d" % scandim
        msglist.append("Scan dimensions: "+msg)
        msg = 'No. of points (outer loop first): '
        for i in range(scandim):
            if i>0 :
                msg = msg + ' x '
            msg = msg + "%d" % round(caget(pvs.ide34scanPtsNos[scandim-i]))
        msglist.append(msg)
        
        # scan1 mode (fly or step)
        scan1mod = caget(pvs.ide34scan1Mode)
        if scan1mod == 2 :
            stat = 'Fly Scan'; style = itfont
        else:
            stat = 'Step Scan'; style = ''
        msg = 'Scan mode: ' + htmlstylespan(stat,style)
        msglist.append(msg)
        
        # scan1 det triggers
        detlist = []
        for i in range(4):
            triggerpv = pvs.ide34scan1TriggerPVs[i]
            trigger = caget(triggerpv)
            if trigger: # if not empty string
                scalerfound = trigger.find(pvs.scalerTrigger)
                PEfound = trigger.find(pvs.PEcommonStr)
                # write expclicitly if one of PE detector
                if PEfound >= 0:
                    PEid = int(trigger[7])  # '1' or '2' or '3'
                    PEname = 'PE detector #' + trigger[7]
                    # grab image number from PE ioc
                    PEfileNStr = '%d' % int(round(caget(pvs.PEfileNs[PEid-1])) - 1)
                    detlist.append(PEname + ' (current file No. = ' + PEfileNStr + ')') 
                # write explicitly if scaler
                elif scalerfound >= 0:
                    detlist.append('Scalers')
                # echo PV if unrecognized
                else:
                    detlist.append(trigger)

        msglist.append("Detectors triggered: " + str(len(detlist)))
        # give a message if detector list is empty
        if not detlist:
            detlist = [('NO detectors are triggered by scan1!', dorangecolor)]
        msglist.append(detlist)
        text += htmlulist(msglist)
    # end of 'scanbusy' condition
    
    
    ## finish up the scan part by making it a form
    #text = htmltag("Scan Status",'legend', {'style':boldfont+bigfont}, True) + text
    #html += htmltag(text,'form',{}, True)
    
    #"form" doesn't really work
    html += htmltag('Scan Status: ', 'h3', {}, True) + text
    
    ### beam part ###
    text = ''
    # beam bad?
    beambad = caget(pvs.ide34beambad)
    if beambad < 0.5: # good
        text += htmlpg(('Beam is OK.',dgreencolor))
    else:             # really bad
        text += htmlpg(('Beam is Bad!',redcolor+boldfont))
    
    msglist = []
    
    # D-shutter
    DShutter = caget(pvs.DShutterClosed)
    if DShutter > 0.5:
        stat = 'closed'; style = redcolor
    else:
        stat = 'open'; style = dgreencolor
    msg = 'D-shutter '+htmlstylespan(stat,style)
    msglist.append(msg)
    
    # FE shutter
    FEShutter = caget(pvs.FEShutterClosed)
    if FEShutter > 0.5:
        stat = 'closed'; style = redcolor
    else:
        stat = 'open'; style = dgreencolor
    msg = 'Front End shutter ' + htmlstylespan(stat,style)
    msglist.append(msg)
    
    # mono/white
    monoMode = caget(pvs.ide34monoMode)
    if monoMode == 1:
        stat = 'white beam'; style = boldfont
    elif monoMode ==2:
        stat = 'monochromatic'; style = itfont
    else:
        stat = 'pass-thru'; style = ''
    msg = 'Monochromator mode: ' + htmlstylespan(stat, style)
    msglist.append(msg)
    
    # Undulator control mode
    idmode = caget(pvs.ID34dsControlMode)
    if int(idmode) == 0 :
        stat = 'User'; style = dgreencolor #<< colorize
    else:
        stat = 'Control Room'; style = redcolor #<< colorize
    msg = 'Undulator control mode: ' + htmlstylespan(stat, style)
    msglist.append(msg)
    
    # id gap (kev in ?th harmonic), taper
    idgap = caget(pvs.ID34dsGapmm)
    msg = 'Undulator gap = %.2f mm.' % idgap
    msglist.append(msg)
    idtaper = caget(pvs.ID34dsTapermm)
    msg = 'Undulator taper = %.2f mm.' % idtaper #<< conditional colorize?
    msglist.append(msg)
    
    # id synced with monoenergy?
    if monoMode ==2:
        idsync = caget(pvs.ide34monoIDSync)
        if idsync == 1 :
            msg = 'Undulator in sync with mono energy.' #<< colorize?
        else:
            msg = 'Undulator NOT in sync with mono energy.' #<< colorize?
        msglist.append(msg)
        
    # storage ring current
    srcurrent = caget(pvs.ringCurrent)
    msg = 'Storage ring current : %.2f mA.' % srcurrent
    msglist.append(msg)
    
    text += htmlulist(msglist)
    
    ## EPS status ##
    eps_idd = caget(pvs.eps34iddOK)
    eps_ida = caget(pvs.eps34idaOK)
    eps_vac = caget(pvs.eps34idVacOK)
    if eps_idd*eps_ida*eps_vac:
        text += htmlpg('Beamline EPS OK.')
    else:
        text += htmlpg(('Beamline EPS Fault:',redcolor))
        epslist = []
        if not eps_vac:
            epslist.append(('Beamline vacuum bad.',redcolor))
        if not eps_ida:
            epslist.append(('Front End Shutter disabled!',redcolor))
        if not eps_idd:
            epslist.append(('D Shutter disabled!',redcolor))
            eps_ideHe = caget(pvs.eps34ideHeOK)
            if eps_ideHe:
                msg = 'E hutch He pressure OK'
            else:
                msg = ('E hutch He pressure low!',redcolor)
            epslist.append([msg])
        text += htmlulist(epslist)
        
    ## finish up the beam part by making it a form
    #text = htmltag("Beam Status",'legend', {'style':boldfont+bigfont}, True) + text
    #html += htmltag(text,'form',{}, True)

    html += htmltag('Beam Status: ', 'h3', {}, True) + text

    return html
        
def extlinks():
    html = htmlpg(('External Links',boldfont))
    linklist = []
    # storage ring history
    linklist.append(htmllink('Storage ring history','http://www.aps.anl.gov/Facility/Status/'))
    # 34ide live cameras (internal access)
    linklist.append(htmllink('34IDE live cameras (internal access)','http://s34video.xor.aps.anl.gov'))
    # 34ide home page
    linklist.append(htmllink('34IDE Home Page','http://www.aps.anl.gov/Sectors/33_34/microdiff/index.html'))
    # ssm group page
    linklist.append(htmllink('SSM Group Home Page','http://www.aps.anl.gov/Sectors/33_34/'))
    
    html += htmlulist(linklist)
    
    return html

def main():
    #imports
    import os
    from datetime import datetime
    
    # start html body, write page title
    print("Content-type: text/html")
    print('\n')
    print("<html>")
    print(htmltag('34IDE Status','title',{}))
    print("<body bgcolor=\"#D6D6D6\">")
    print(htmltag('Current Status of 34IDE','h2',{}))
    message = datetime.now().strftime("%I:%M:%S %p, %a %Y-%m-%d")
    message = htmlpg(("Time of display: " + message, smallfont+itfont))
    print(message)
    
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

    # finally: provide useful links
    print("<hr />\n")
    print(extlinks())
    
    # finish html body & return
    print("</body>")
    print("</html>")
    return
    
#main()

if __name__ == '__main__':
    #msglist = [('item1','color:grey'),'item2',(['item31','item32'],'color:yellow')]
    #message = htmlulist(msglist)
    #print message
    #message = 'a para\ngraph'
    #print htmlpg((message,'color:orange'))
    #main()
    pass
    