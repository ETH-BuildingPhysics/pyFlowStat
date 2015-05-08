# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 16:28:00 2015

@author: mimmer
"""

import h5py

def saveDict(filename,dictionary,keyList=[],mode='w',dictName='dict'):
    '''
    Save dict as H5
    mode: 'w' for overwrite, 'w-' for new, 'a' for append
    '''
    if len(keyList)==0:
        keyList = dictionary.keys()
    keys = dictionary.keys()
    fwm = h5py.File(filename, mode)
    if dictName in fwm:
        gDict=fwm[dictName]
    else:
        gDict = fwm.create_group(dictName)
    for k in [k for k in keys if k in keyList]:
        gDict.create_dataset(k,data=dictionary[k])
    fwm.close()
    
def loadDict(filename,keyList=[],dictName='dict'):
    '''
    Load dict as H5
    '''
    dictionary=dict()
    fwm = h5py.File(filename, 'r')
    
    keys = fwm[dictName].keys()
    
    if len(keyList)==0:
        keyList=keys
    
    for k in [k for k in keys if k in keyList]:
        dictionary[k]=fwm[dictName][k].value
        
    fwm.close()
    return dictionary