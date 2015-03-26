def sortNumStrList(numStrList):
    '''
    Sort a list of number stored as a list of string. StrNumList is the list of
    number as string.
    
    Arguments:
        *numStrList*
    '''
    #tsListNum=[]
    numFltList = []
    for nb in numStrList:
        if is_number(nb)==True:
            numFltList.append(float(nb))
    numFltList = np.array(numFltList)
    numStrList_idx = np.argsort(numFltList)
    numStrList_sort = [numStrList[idx] for idx in numStrList_idx]
    return numStrList_sort


def is_number(string):
    '''
    
    '''
    try:
        float(string)
        return True
    except ValueError:
        return False