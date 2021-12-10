"""
functions that takes various forms of string inputs and generate html texts

-- Ruqing Xu, Apr. 2012
"""

def htmltag(content, tag, attributes, add_newline = False):
    """ Convert a string content into a tagged html text.
        content, tag - strings
        attributes - dict of {'attrname':stringvalue,...}
        add_newline - optional, if set true, will add a "\n" at the end of the text
    """
    attrstring = ''
    for attr in list(attributes.keys()):
        val = attributes[attr]
        attrstring += ' ' + str(attr)+'=\"' + str(val) + '\"'
    
    html = '<' + str(tag) + attrstring + '>' + content + '</' + str(tag.strip()) + '>'
    if add_newline:
        html += '\n'
    return html
    

def htmlulist(msgin):
    """ Convert a list of strings into a string for html unordered list display; 
        input list may be nested.

        Either the entire input or each individual list item may also be a 2-tuple,
        in which case the 1st item will be the display content and the 2nd item will 
        be a string giving the style specification.
        
        Example: to generate the following html list:
            <ul style = "color:red;font-weight:bold">
              <li>aa</li>
              <li style = "font-style:italic">bb</li>
              <ul>
                <li>cc</li>
              </ul>
            </ul>
        you should call the function as:
            htmllist((['aa',('bb','font-style:italic'),['cc']],'color:red;font-weight:bold'))
        
        TODO: type checking & raising error message for each element
    """
    # if input is a list
    if type(msgin) == list :
        html = '<ul>\n'
        msg = msgin
    # if input is (list, style_string)
    elif (type(msgin) == tuple) and (len(msgin) ==2) and (type(msgin[0]) == list):
        msg = msgin[0]
        style = msgin[1]
        html = '<ul style=\"' + style + '\">\n'
        
    for ss in msg: # loop over list items
        # if a sub-list
        if isinstance(ss,list): 
            html += htmlulist(ss)
        # if a sub-list with style specification
        elif isinstance(ss,tuple) and isinstance(ss[0],list) and (len(ss) ==2) :
            html += htmlulist(ss)
        # if just a string
        elif isinstance(ss,str) :
            html += htmltag(ss,'li',{},True) #'<li>' + ss + '</li>\n'
        # if a string with style specification
        elif isinstance(ss,tuple) and isinstance(ss[0],str) and (len(ss) ==2) :
            (txt, style) = ss
            html += htmltag(txt,'li',{'style':style},True) #'<li style=\"' + style + '\">' + txt + '</li>\n'
    html += '</ul>\n'

    return html

def htmlpg(msg):
    """ Convert a string into an html pragraph.
    
        The input message may also be a tuple of 2 strings, in which case the 1st string
        will be the display content and the 2nd will give the style specification.
        
        See the docstring of htmlulist for example of style specifications.
    """
    if type(msg) == str :
        txt = msg
        attr = {}
    elif (type(msg) == tuple) and (len(msg) ==2):
        (txt, style) = msg
        attr = {'style':style}
    else:
        return 'Unknown nput type'

    #if br_replace :
    #    txt = txt.replace('\r\n','<br />')
    #    txt = txt.replace('\n','<br />')
    
    return htmltag(txt, 'p', attr, True)
    
def htmlstylespan(msg,style=''):
    return htmltag(msg, 'span', {'style':style})

def htmllink(txt,link):
    """ give html format of a link represented by txt.
    """
    return htmltag(txt,'a',{'href':link},True)