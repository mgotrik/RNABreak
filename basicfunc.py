##### Called by Puzzle Definition: (in ~order)
from itertools import chain, count, groupby
import sys, os
import more_itertools as mit #allows to get contiguous items in a list

### Basic functions:

#   Inputs iterable, outputs iterable grouped by consecutive bases (required)
def find_ranges(iterable):
    """Yield range of consecutive numbers."""
    for group in mit.consecutive_groups(iterable):
        group = list(group)
        if len(group) == 1:
            yield group[0]
        else:
            yield group[0], group[-1]

#   Inputs a list, outputs the number of nested lists (required)
def depth(seq):
    for level in count():
        if not seq:
            return level
        seq = list(chain.from_iterable(s for s in seq if isinstance(s, Sequence)))

def getPairsFromList(List):
    """
    Inputs a flat list, groups it into pairs. Original version did no checks;
    this one asserts on list length
    """
    assert(len(list) % 2 == 0)
    return zip(List[::2], List[1::2])

#   Breaks all nested lists, in order, to give a single list
def FlattenList(List):
    flat_list = []
    flag=0
    try:
        len(List)
    except:
        List = [List]

    while flag==0:
        flag=1
        for i in range(len(List)):
            if type(List[i]) != list:
                flat_list.append(List[i])
            else:
                flag=0
                for j in List[i]:
                    flat_list.append(j)
        if flag==0:
            List=list(flat_list)
            flat_list=[]
    return flat_list

#   Removes duplicates, orders and flattens a list of HPs
def CleanList(ListOfHPs):
    return list(set(FlattenList(ListOfHPs)))

# Disable printing
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

# Restore
def enablePrint():
    sys.stdout = sys.__stdout__


##### HPGroups ++

def getAttributeOfFeatFromMemberID(Feat,MemberIDSearchAttr,MemberIDList,ReturnAttribute,anyorall='any',FeatIDList=None):
    ReturnList = []
    if FeatIDList==None:
        FeatIDList = Feat._IDlist
    MemberIDList = list(MemberIDList)
    for i in FeatIDList:
        if anyorall=='any':
            # print("Looking for", MemberIDList, "inside of", getattr(Feat[i],MemberIDSearchAttr))
            # if MemberIDList in getattr(Feat[i], MemberIDSearchAttr):
            #     print('yes')
            # print("Hint:", any(elem in getattr(Feat[i],MemberIDSearchAttr) for elem in MemberIDList))
            if any(elem in MemberIDList for elem in list(getattr(Feat[i],MemberIDSearchAttr))) or any(elem in [MemberIDList] for elem in list(getattr(Feat[i],MemberIDSearchAttr))):
                ReturnList.append(getattr(Feat[i],ReturnAttribute))
        else:
            if all(elem in MemberIDList for elem in getattr(Feat[i],MemberIDSearchAttr)):
                ReturnList.append(getattr(Feat[i],ReturnAttribute))
            # SearchList = [x for x in getattr(Feat[i], MemberIDSearchAttr) if x in MemberIDList]
            # if all(elem in MemberIDList for elem in SearchList):
            #     ReturnList.append(getattr(Feat[i],ReturnAttribute))
    return ReturnList

def getPropertyOfFeatsWithProperty(Feat,attribute,attributeout, value=None,IDList=None):
    if IDList==None:
        IDList = Feat._IDlist
    if value==None:
        return [getattr(Feat[x],attributeout) for x in IDList]
    AttributesIn = [getattr(Feat[x],attribute) for x in IDList]
    AttributesOut = [getattr(Feat[IDList[x]],attributeout) for x in range(0,len(AttributesIn)) if AttributesIn[x]==value]
    return AttributesOut

def groupNearbyFeaturesBasedOnAttribute(Feat,FeatID, attribute, value, nearorglobal, IDList=None):
#Searches a feature, starting at FeatID, for attribute of value, nearby, or globally, choosing from the IDList
    try:
        nearorglobal
    except:
        nearorglobal='near'
    HPList = [FeatID]
    ParentAttrValue = getattr(Feat[FeatID],attribute)
    if value is not None:
        if ParentAttrValue!=value:
            return [], []

    if nearorglobal=='near': #Only looks at features that are nearby, not globally
        NearbyFeats = Feat[FeatID].neighborID
        #NearbyFeats = Feat[FeatID].parentID + Feat[FeatID].childrenID
        RemainingHPs = list(NearbyFeats)
    else:
        RemainingHPs = list(Feat._IDlist)
        if IDList is not None:
            RemainingHPs = bf.FlattenList(IDList)
    if IDList is None:
        IDList = list(Feat._IDlist)

    while RemainingHPs!=[]:
        CurrHP = RemainingHPs.pop(0)
        #print("Currhp",CurrHP,RemainingHPs)
        CurrHPAttr = getattr(Feat[CurrHP],attribute)
        if CurrHPAttr==ParentAttrValue and nearorglobal=='near' and CurrHP not in HPList and CurrHP in IDList:
            HPList.append(CurrHP)
            NewHP = [x for x in Feat[CurrHP].neighborID if getattr(Feat[x],attribute)==ParentAttrValue and x not in HPList]
            #print("found new HP", NewHP, CurrHP)
            RemainingHPs.insert(0,NewHP)
            RemainingHPs = FlattenList(RemainingHPs)
        elif CurrHPAttr==ParentAttrValue and CurrHP not in HPList:
            HPList.append(CurrHP)
    return HPList,ParentAttrValue

def groupIDListBasedOnAttribute(Feat, attribute, value, nearorglobal=None, IDList=None):

#Searches a feature, starting at FeatID, for attribute of value, nearby, or globally, choosing from the IDList
    if nearorglobal!=None:
        nearorglobal='near'
    if nearorglobal=='global':
        if IDList is not None:
            AttrValues = [getattr(Feat[x],attribute) for x in IDList]
            UniqueValues = list(set(AttrValues))
            return [[[IDList[x] for x in range(0,len(IDList)) if AttrValues[x]==y] for y in UniqueValues], UniqueValues]
        else:
            IDList = Feat._IDlist
            AttrValues = [getattr(Feat[x],attribute) for x in IDList]
            UniqueValues = list(set(AttrValues))
            return [[[x for x in IDList if AttrValues[x]==y] for y in UniqueValues], UniqueValues]

    if IDList is not None:
        RemainingFeats = FlattenList(IDList)
    else:
        IDList = list(Feat._IDlist)
        RemainingFeats = list(Feat._IDlist)

    FeatureMap = []
    AttrValues = []
    for i in IDList:
        if i not in RemainingFeats:
            continue
        ParentFeatAttr = getattr(Feat[i], attribute)
        NearbyFeats = groupNearbyFeaturesBasedOnAttribute(Feat, i, attribute, ParentFeatAttr, 'near', IDList)[0]
        NearbyFeats = sorted([x for x in NearbyFeats if x in IDList and x in RemainingFeats])
        if NearbyFeats!=[]:
            RemainingFeats = [x for x in RemainingFeats if x not in NearbyFeats]
            CurrFeatAttr = getattr(Feat[NearbyFeats[0]], attribute)
            FeatureMap.append(NearbyFeats)
            AttrValues.append(CurrFeatAttr)
        else:
            FeatureMap.append([i])
            AttrValues.append(ParentFeatAttr)

    if value is None:
        return FeatureMap,AttrValues
    else:
        FeatureMap = [x if AttrValues[FeatureMap.index(x)]==value else [] for x in FeatureMap]
        AttrValues = [x for x in AttrValues if FeatureMap[AttrValues.index(x)]!=[]]
        FeatureMap = [x for x in FeatureMap if x!=[]]
        return FeatureMap, AttrValues

def groupAllFeaturesBasedOnAttribute(Feat,attribute,value=None,nearorglobal='near',IDList=None):
    if IDList is None:
        RemainingFeats = list(Feat._IDlist)
    else:
        RemainingFeats = FlattenList(IDList)
    FeatureMap = []
    AttrValues = []
    while RemainingFeats!=[]:
        CurrFeat = RemainingFeats.pop(0)
        NearbyFeats = groupNearbyFeaturesBasedOnAttribute(Feat,CurrFeat,attribute, value, nearorglobal, IDList)
        AttrValues.append(NearbyFeats[1])
        NearbyFeats = NearbyFeats[0]
        if NearbyFeats!=[]:
            RemainingFeats = [x for x in RemainingFeats if x not in NearbyFeats]
            FeatureMap.append(NearbyFeats)
        else:
            FeatureMap.append(CurrFeat)

    if value is not None:
        FeatureMap = [x for x in FeatureMap if AttrValues[FeatureMap.index(x)]!=[]]
        AttrValues = [x for x in AttrValues if x!=[]]

    return FeatureMap,AttrValues

