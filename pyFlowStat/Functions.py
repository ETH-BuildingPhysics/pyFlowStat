import numpy as np

def sortNumStrList(numStrList,minVal=None,maxVal=None,step=1):
    '''
    Sort a list of number stored as a list of string. StrNumList is the list of
    number as string.
    
    Arguments:
        *numStrList*
    '''
    #tsListNum=[]
    if len(numStrList)==1:
        return numStrList
    numFltList = []
    for nb in numStrList:
        if is_number(nb)==True:
            numFltList.append(float(nb))
    numFltList = np.array(numFltList)
    numStrList_idx = np.argsort(numFltList)
    #numFltList_sorted = np.sort(numFltList)
    
    if minVal and maxVal:
        numStrList_idx_filt=[]
        for idx in numStrList_idx:
            f=numFltList[idx]
            if f>=minVal and f<=maxVal:
                numStrList_idx_filt.append(idx)
    else:
        numStrList_idx_filt=numStrList_idx
    
    if not step:
        step=1.0
    
    numStrList_sort = [numStrList[idx] for idx in numStrList_idx_filt][::step]
    return numStrList_sort


def is_number(string):
    '''
    
    '''
    try:
        float(string)
        return True
    except ValueError:
        return False