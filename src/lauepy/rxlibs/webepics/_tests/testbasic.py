#!/APSshare/epd/SunOS_5.10-x86/bin/python
#

for i in range(4):
    triggerpv = '34ide:scan1.T' + ('%d' % (i+1)) + 'PV'
    print(triggerpv)
    