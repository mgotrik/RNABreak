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

#   Inputs a flat list, groups it into pairs
def getPairsFromList(List):
    NewList = [[List[0],List[1]]]
    NumPairs = int(len(List) / 2)
    if NumPairs > 1:
        for i in range(1, NumPairs):
            Idx1 = List[i * 2]
            Idx2 = List[i * 2 + 1]
            NewList.append([Idx1,Idx2])
    return NewList

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
    ReturnList = []
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

##### HPs + Bulges
#CRAWLING AROUND SEQUENCE


# def iterStep(Sequence,StartIdx,StepIter):
#     maxlength = len(Sequence)-1
#     NextIdx = StartIdx + StepIter
#     if NextIdx<0:
#         NextIdx = maxlength
#     elif NextIdx>=maxlength:
#         NextIdx = 0
#     return NextIdx

# def getContinuousBases(Sequence,PartnerIdxs,CurrIndex,StepIter):
#
#     if CurrIndex<0:
#         CurrIndex = 0
#     elif CurrIndex >=len(Sequence):
#         CurrIndex = len(Sequence)-1
#
#     if CurrIndex!=0 and CurrIndex!=len(Sequence)-1:
#         PrevIndex = iterStep(Sequence,CurrIndex,-StepIter)
#         while PartnerIdxs[PrevIndex]==-1 and CurrIndex>0 and CurrIndex<len(Sequence):
#             CurrIndex = PrevIndex
#             PrevIndex = iterStep(Sequence,CurrIndex,-StepIter)
#     IdxList = []
#     NextIdx = iterStep(Sequence,CurrIndex,StepIter)
#     #print("Checking", PartnerIdxs[CurrIndex],PartnerIdxs[NextIdx])
#
#     while PartnerIdxs[CurrIndex] == -1 and PartnerIdxs[NextIdx] == -1:
#         IdxList.append([CurrIndex, iterStep(Sequence, CurrIndex, StepIter)])
#         CurrIndex = NextIdx
#         NextIdx = iterStep(Sequence,CurrIndex,StepIter)
#
#     IdxList = sorted(CleanList(IdxList))
#
#     return IdxList, NextIdx

# def getAllSurroundingBPs(Sequence, StartIdx, Direction_Flag):
#     #print("Start Pair", StartPair)
#     PartnerIdxs = getPartnerIdxs(Sequence[2])[0]
#
#     if Direction_Flag==0:
#         StepIter = -1
#     else:
#         StepIter = 1
#
#     try:
#         CurrIndex = StartIdx[0]
#     except:
#         CurrIndex = StartIdx
#
#
#     AllIndex = []
#     #AllIndex.append(StartPair)
#     if PartnerIdxs[CurrIndex]!=-1:
#         CurrIndex = min([CurrIndex,PartnerIdxs[CurrIndex]])
#         AllIndex.append([min([CurrIndex,PartnerIdxs[CurrIndex]]),max([CurrIndex,PartnerIdxs[CurrIndex]])])
#
#
#     #print(StartIdx,"CI:",CurrIndex, PartnerIdxs[CurrIndex],CheckBulgeBoundaries(Sequence,CurrIndex,Direction_Flag),Direction_Flag)
#     #print("min",min([CurrIndex,PartnerIdxs[CurrIndex]]),AllIndex)
#     CurrIndex = iterStep(Sequence,CurrIndex,StepIter)
#     # if PartnerIdxs[CurrIndex]!=-1:
#     #     NewPair = getPairFromBaseID(Sequence[2], CurrIndex)
#     #     # print("found pair,", NewPair, CurrIndex)
#     #     CurrIndex = [x for x in NewPair if x != CurrIndex][0]
#     #     # print("found pair,", NewPair, CurrIndex)
#     #     CurrIndex = iterStep(Sequence, CurrIndex, StepIter)
#     #     if NewPair in AllIndex:
#     #         return sorted(AllIndex)
#     #     else:
#     #         AllIndex.append(NewPair)
#
#     while CheckBulgeBoundaries(Sequence,CurrIndex,Direction_Flag):
#         #print("Entering once,", CurrIndex)
#         if PartnerIdxs[CurrIndex] == -1:
#             CurrIndex = iterStep(Sequence,CurrIndex,StepIter)
#         else:
#             NewPair = getPairFromBaseID(Sequence[2], CurrIndex)
#             #print("found pair,", NewPair, CurrIndex)
#             #CurrIndex = [x for x in NewPair if x!=CurrIndex][0]
#             CurrIndex = PartnerIdxs[CurrIndex]
#             #print("found pair,", NewPair, CurrIndex)
#             #CurrIndex = iterStep(Sequence, CurrIndex, StepIter)
#             CurrIndex = iterStep(Sequence,CurrIndex,StepIter)
#             if NewPair in AllIndex:
#                 break
#             else:
#                 AllIndex.append(NewPair)
#     #print("Found all for", StartIdx,AllIndex)
#
#     return sorted(AllIndex)
#
#
# def getAllSurroundingBPs2(Sequence, StartIdx, Direction_Flag):
#     #print("Start Pair", StartPair)
#     PartnerIdxs = getPartnerIdxs(Sequence[2])[0]
#
#     if Direction_Flag==0:
#         StepIter = -1
#     else:
#         StepIter = 1
#
#     try:
#         CurrIndex = StartIdx[0]
#     except:
#         CurrIndex = StartIdx
#     #print("Currinedex:", CurrIndex, type(CurrIndex))
#     AllIndex = []
#     #AllIndex.append(StartPair)
#     CurrIndex = iterStep(Sequence,CurrIndex,StepIter)
#     # if PartnerIdxs[CurrIndex]!=-1:
#     #     NewPair = getPairFromBaseID(Sequence[2], CurrIndex)
#     #     # print("found pair,", NewPair, CurrIndex)
#     #     CurrIndex = [x for x in NewPair if x != CurrIndex][0]
#     #     # print("found pair,", NewPair, CurrIndex)
#     #     CurrIndex = iterStep(Sequence, CurrIndex, StepIter)
#     #     if NewPair in AllIndex:
#     #         return sorted(AllIndex)
#     #     else:
#     #         AllIndex.append(NewPair)
#
#     while CheckBulgeBoundaries(Sequence,CurrIndex,Direction_Flag):
#
#         #print("Entering once,", CurrIndex)
#         if PartnerIdxs[CurrIndex] == -1:
#             CurrIndex = iterStep(Sequence,CurrIndex,StepIter)
#         else:
#             NewPair = getPairFromBaseID(Sequence[2], CurrIndex)
#             #print("found pair,", NewPair, CurrIndex)
#             CurrIndex = [x for x in NewPair if x!=CurrIndex][0]
#             #print("found pair,", NewPair, CurrIndex)
#             CurrIndex = iterStep(Sequence, CurrIndex, StepIter)
#             if NewPair in AllIndex:
#                 return sorted(AllIndex)
#             else:
#                 AllIndex.append(NewPair)
#     return sorted(AllIndex)
# #This is the old one that I was using but has problems when the last base is paired:
# def getAllSurroundingBPs4(Sequence, StartIdx, Direction_Flag):
#     #print("Start Pair", StartPair)
#     PartnerIdxs = getPartnerIdxs(Sequence[2])[0]
#
#     try:
#         CurrIndex = StartIdx[0]
#     except:
#         CurrIndex = StartIdx
#     #print("Currinedex:", CurrIndex, type(CurrIndex))
#     AllIndex = []
#     #AllIndex.append(StartPair)
#     if Direction_Flag==1:
#         CurrIndex += 1
#         #print("CID1",CurrIndex,PartnerIdxs[CurrIndex])
#         while CurrIndex<len(Sequence):
#             if PartnerIdxs[CurrIndex]==-1:
#                 CurrIndex+=1
#             else:
#                 NewPair = getPairFromBaseID(Sequence[2],CurrIndex)
#                 CurrIndex = PartnerIdxs[CurrIndex]+1
#                 if NewPair in AllIndex:
#                     return sorted(AllIndex)
#                 else:
#                     AllIndex.append(NewPair)
#             if CurrIndex>=len(Sequence):
#                 CurrIndex = 0
#     else:
#         if CurrIndex==0:
#             if PartnerIdxs[CurrIndex]!=-1:
#                 return getPairFromBaseID(Sequence[2],CurrIndex)
#             else:
#                 CurrIndex = iterStep(Sequence,CurrIndex,0)
#
#         #CurrIndex-=1
#         while CurrIndex>=0:
#             #print("CID2",CurrIndex,PartnerIdxs[CurrIndex])
#             if PartnerIdxs[CurrIndex]==-1:
#                 CurrIndex-=1
#             else:
#                 NewPair = getPairFromBaseID(Sequence[2],CurrIndex)
#                 CurrIndex = PartnerIdxs[CurrIndex]-1
#                 #print("found pair",NewPair,AllIndex,CurrIndex)
#                 if NewPair in AllIndex:
#                     #print("sorting", AllIndex, sorted(AllIndex))
#                     return sorted(AllIndex)
#                 else:
#                     AllIndex.append(NewPair)
#             if CurrIndex<=0:
#                 CurrIndex = len(Sequence)-1

# def CheckBulgeBoundaries(Sequence,Index,DirectionFlag):
#     if DirectionFlag==0:
#         return Index < len(Sequence)
#     else:
#         return Index >= 0




# def getHPFromBaseIDList(BaseID,Sequence):
#     PairMap = getHPMap(Sequence[2])
#     IncludedHPs = []
#     for base in BaseID:
#         #IncludedHPs += [PairMap.index(x) for x in PairMap if x[0][0]<=base<=x[0][1] or x[1][0]<=base<=x[1][1]]
#         IncludedHPs += [PairMap.index(x) for x in PairMap if x[0][0]<=base<=x[-1][0] or x[0][1]<=base<=x[-1][1]]
#     IncludedHPs = CleanList(IncludedHPs)
#     return IncludedHPs

# def getPairFromBaseID(SS,BaseID):
#     List_Of_Partner_Idx, List_Of_Pairs = getPartnerIdxs(SS)
#     if type(BaseID)==list:
#         PairMap = []
#         for i in BaseID:
#             PairMap.append([[x,y] for [x,y] in List_Of_Pairs if i in [x,y]][0])
#         return PairMap
#     else:
#         return [Pair for Pair in List_Of_Pairs if BaseID in Pair][0]
#
# def getAllSurroundingHPs(Sequence,StartPair,Direction_Flag=1):
#
#     SurroundingPairs = getAllSurroundingBPs(Sequence, StartPair, Direction_Flag)
#     if SurroundingPairs==None:
#         SurroundingPairs = getAllSurroundingBPs(Sequence, StartPair, 1 - Direction_Flag)
#         if SurroundingPairs==None:
#             print("Error SurroundingPairs!", SurroundingPairs , Sequence,StartPair,Direction_Flag)
#             return None
#     SurroundingHPs = getHPFromBaseIDList(CleanList(SurroundingPairs), Sequence)
#
#     return SurroundingHPs
#
# def getLeafHPsinHPList(FG,HPList):
#     HPList = FlattenList(HPList)
#     LeafHPs = [x for x in HPList if any(elem in FG.HP[x].childrenID for elem in HPList)==False]
#     return LeafHPs













####    Weighting
# def _initWeight(Feat):
#     for i in Feat._IDlist:
#         #Feat._weightmap[i] = 0
#         Feat[i].Weight=0
#     print("initialized weight map", Feat.name)




# def getBulgeMap2(HP):
#
#     # We are starting at the first HP, so the first bulge has children 0
#     Sequence = HP.parentseq
#     ParentIDList = [[]]
#     NextHP = [0]
#     ChildrenIDList = [[0]]
#
#     if len(HP[0].sibling)==0:
#         # This assigns the outermost bulge, as the end of the puzzle to the start of the first HP
#         BulgeList = [list(range(0, HP[0].idx[0][0])), list(range(HP[0].idx[0][1] + 1, len(Sequence)))]
#         IDlist = list(HP._IDlist)
#         BulgeListAll = [BulgeList]  # List containing all bulges separately
#
#         while NextHP != []:  # While we don't find any more to add:
#             CurrHP = NextHP[0]
#             ChildrenHP = HP[CurrHP].childrenID
#             if ChildrenHP == []:
#                 BulgeList = list(range(HP[CurrHP].idx[-1][0] + 1, HP[CurrHP].idx[-1][1]))
#             else:
#                 PairList = sorted(FlattenList(getAllSurroundingBPs(Sequence, HP[CurrHP].idx[-1], 1)))
#                 PairList = getPairsFromList(PairList)
#                 BulgeList = []
#                 for i in PairList:
#                     BulgeList.append(list(range(i[0] + 1, i[1])))
#                 NextHP.insert(0, ChildrenHP)
#             NextHP.remove(CurrHP)
#             NextHP = CleanList(NextHP)
#             BulgeListAll.append(BulgeList)
#             ParentIDList.append([CurrHP])
#             ChildrenIDList.append(ChildrenHP)
#
#         return BulgeListAll, ParentIDList, ChildrenIDList


# def getBulgeListOfBases(Sequence):
#     # Inputs Sequence (for SS) and returns list of all unpaired bases that form a bulge
#     #Reqd Functions:
#     #   getAllSurroundingBPs
#     #   getPairsFromList
#     #   FlattenList
#     # Gets the sequence
#     SS = Sequence[2]
#     #Gets PartnerIndex list
#     PartnerIndexes = getPartnerIdxs(SS)[0]
#     #Identifies all bulged bases
#     BulgedBases = [i for i, x in enumerate(PartnerIndexes) if x == -1]
#     RemainingBases = list(BulgedBases)
#
#     UnpairedBasesInBulgeList = []
#     AllBasesInBulgeList = []
#
#     #iterates through every bulged base
#     for b in BulgedBases:
#         #Makes sure we are not at the start/end of puzzle
#         try:
#             if 0<=b<len(Sequence)-1:
#                 SurroundingBPs = getAllSurroundingBPs(Sequence,[b],1)
#             else:
#                 SurroundingBPs = getAllSurroundingBPs(Sequence,[b],0)
#
#         except:
#             continue
#         #Identifies bases that are between the BPs:
#         RegionsBetweenHPs = getPairsFromList(sorted(FlattenList(SurroundingBPs)))
#         UnpairedBasesInRegions = []
#         AllBasesInRegions = []
#         for i in RegionsBetweenHPs:
#             for j in range(i[0]+1,i[1]):
#                 UnpairedBasesInRegions +=[j]
#             for j in range(i[0], i[1]+1):
#                 AllBasesInRegions += [j]
#             #print(i, UnpairedBasesInRegions)
#         #Adds these to BulgeList containing all bases that make up a unique bulge
#         if UnpairedBasesInRegions!= [] and UnpairedBasesInRegions[0] in RemainingBases:
#             #print("made it",UnpairedBasesInRegions)
#             for k in UnpairedBasesInRegions:
#                 RemainingBases.remove(k)
#             UnpairedBasesInBulgeList += [UnpairedBasesInRegions]
#             AllBasesInBulgeList +=[AllBasesInRegions]
#
#     if RemainingBases!=[]:
#         FinalBulge = list(RemainingBases)
#         try:
#             FinalBulge += getAllSurroundingBPs(Sequence,[RemainingBases[0]],1)
#         except:
#             FinalBulge += getAllSurroundingBPs(Sequence,[RemainingBases[0]],0)
#         FinalBulge = sorted(FlattenList(FinalBulge))
#         UnpairedBasesInBulgeList = [RemainingBases] + UnpairedBasesInBulgeList
#         AllBasesInBulgeList = [FinalBulge] + AllBasesInBulgeList
#
#     #Sorts both of them by their order from the start:
#     AllBasesInBulgeList,UnpairedBasesInBulgeList = [list(t) for t in zip(*sorted(zip(AllBasesInBulgeList,UnpairedBasesInBulgeList)))]
#
#     return UnpairedBasesInBulgeList,AllBasesInBulgeList

# def getPairMap(SS):
#     # Requires Variables:
#     # SS
#     # Requires Functions:
#     # getPartnerIdxs
#     List_Of_Partner_Idx, List_Of_Pairs = getPartnerIdxs(SS)
#     tmp_List_Of_Pairs = List_Of_Pairs
#     List_Of_HPs = []
#     List_Of_PairMaps = []
#     i = 0
#     while i < len(tmp_List_Of_Pairs):  # Identifies pairs that are contiguous
#         tmp_List_Of_HPs = [tmp_List_Of_Pairs[i]]
#         HP_Len_Count = 1  # counting how long the hairpin is
#         New_HP_Flag = 0
#         while New_HP_Flag == 0 and i + 1 < len(tmp_List_Of_Pairs):
#             New_HP_Flag = 1
#             if tmp_List_Of_Pairs[i][0] == tmp_List_Of_Pairs[i + 1][0] - HP_Len_Count and tmp_List_Of_Pairs[i][1] == \
#                     tmp_List_Of_Pairs[i + 1][1] + HP_Len_Count:
#                 tmp_List_Of_HPs.append(tmp_List_Of_Pairs[i + 1])
#                 tmp_List_Of_Pairs.remove(tmp_List_Of_Pairs[i + 1])
#                 HP_Len_Count = HP_Len_Count + 1
#                 New_HP_Flag = 0
#
#         List_Of_HPs.append([tmp_List_Of_HPs[0], tmp_List_Of_HPs[-1]])
#         i = i + 1
#     for i in range(0, len(List_Of_HPs)):  # Re-organizes to be separated by strand
#         List_Of_PairMaps.append(
#             [[List_Of_HPs[i][0][0], List_Of_HPs[i][1][0]], [List_Of_HPs[i][1][1], List_Of_HPs[i][0][1]]])
#     return List_Of_PairMaps  # List of all List_Of_Partner_Idxs, list of all pairs, list of pairs sorted by HP





##### Used by functions in this file:

# def getAllNearbyBPAndBulge(Sequence,StartIdx,Direction_Flag):
#     #print("Start Pair", StartPair)
#     PartnerIdxs = getPartnerIdxs(Sequence[2])[0]
#     BulgedBases,UnpairedBases = getBulgeListOfBases(Sequence)
#
#     if Sequence[2][StartIdx]=='.':
#         BulgeIdx = BulgedBases.index([x for x in BulgedBases if StartIdx in x])
#         BPs = [x for x in PartnerIdxs if x[0] in BulgedBases[BulgeIdx]]
#         return UnpairedBases[BulgeIdx], BPs
#     else:
#         AllPairs = getAllSurroundingBPs(Sequence,StartIdx,Direction_Flag)
#         BulgeIdx = BulgedBases.index([x for x in BulgedBases if CleanList(AllPairs)[0] in x])
#         return UnpairedBases[BulgeIdx], AllPairs

# def getSurroundingFeatureInStructure(Sequence,Index,Direction_Flag,PairMap):
#     # print("Start Pair", StartPair)
#     PartnerIdxs = getPartnerIdxs(Sequence[2])[0]
#
#     try:
#         CurrIndex = Index[0]
#     except:
#         CurrIndex = Index
#     # print("Currinedex:", CurrIndex, type(CurrIndex))
#     AllIndex = []
#     # AllIndex.append(StartPair)
#     if Direction_Flag == 1:
#         if CurrIndex == len(Sequence) - 1:
#             CurrIndex = 0
#
#         CurrIndex += 1
#         print("getIndexOfElement",getIndexOfElement(PairMap,CurrIndex))
#         while CurrIndex < len(Sequence):
#             if getIndexOfElement(PairMap,CurrIndex) == []:
#                 if CurrIndex not in AllIndex and CurrIndex in FlattenList(PairMap):
#                     AllIndex.append(CurrIndex)
#                 CurrIndex += 1
#             else:
#                 NewPair = getPairFromBaseID(Sequence[2], CurrIndex)
#                 CurrIndex = PartnerIdxs[CurrIndex] + 1
#                 if NewPair in AllIndex:
#                     return SortIndexList(AllIndex)
#                 else:
#                     AllIndex.append(NewPair)
#             if CurrIndex >= len(Sequence):
#                 CurrIndex = 0
#     else:
#         if CurrIndex == 0:
#             CurrIndex = len(Sequence) - 1
#         CurrIndex -= 1
#         while CurrIndex >= 0:
#             # print("CID2",CurrIndex,PartnerIdxs[CurrIndex])
#             if PartnerIdxs[CurrIndex] == -1:
#                 if CurrIndex not in AllIndex and CurrIndex in FlattenList(PairMap):
#                     AllIndex.append(CurrIndex)
#                 CurrIndex -= 1
#             else:
#                 NewPair = getPairFromBaseID(Sequence[2], CurrIndex)
#                 CurrIndex = PartnerIdxs[CurrIndex] - 1
#                 # print("found pair",NewPair,AllIndex,CurrIndex)
#                 if NewPair in AllIndex:
#                     return SortIndexList(AllIndex)
#                 else:
#                     AllIndex.append(NewPair)
#             if CurrIndex <= 0:
#                 CurrIndex = len(Sequence) - 1

# def getSurroundingStructure(Sequence,Index,Direction_Flag):
#     # print("Start Pair", StartPair)
#     PartnerIdxs = getPartnerIdxs(Sequence[2])[0]
#
#     try:
#         CurrIndex = Index[0]
#     except:
#         CurrIndex = Index
#     # print("Currinedex:", CurrIndex, type(CurrIndex))
#     AllIndex = []
#     # AllIndex.append(StartPair)
#     if Direction_Flag == 1:
#         if CurrIndex==len(Sequence)-1:
#             CurrIndex = 0
#
#         CurrIndex += 1
#         # print("CID1",CurrIndex,PartnerIdxs[CurrIndex])
#         while CurrIndex < len(Sequence):
#             if PartnerIdxs[CurrIndex] == -1:
#                 if CurrIndex not in AllIndex:
#                     AllIndex.append(CurrIndex)
#                 CurrIndex += 1
#             else:
#                 NewPair = getPairFromBaseID(Sequence[2], CurrIndex)
#                 CurrIndex = PartnerIdxs[CurrIndex] + 1
#                 if NewPair in AllIndex:
#                     return SortIndexList(AllIndex)
#                 else:
#                     AllIndex.append(NewPair)
#             if CurrIndex >= len(Sequence):
#                 CurrIndex = 0
#     else:
#         if CurrIndex==0:
#             CurrIndex = len(Sequence)
#         CurrIndex -= 1
#         while CurrIndex >= 0:
#             # print("CID2",CurrIndex,PartnerIdxs[CurrIndex])
#             if PartnerIdxs[CurrIndex] == -1:
#                 if CurrIndex not in AllIndex:
#                     AllIndex.append(CurrIndex)
#                 CurrIndex -= 1
#             else:
#                 NewPair = getPairFromBaseID(Sequence[2], CurrIndex)
#                 CurrIndex = PartnerIdxs[CurrIndex] - 1
#                 # print("found pair",NewPair,AllIndex,CurrIndex)
#                 if NewPair in AllIndex:
#                     return SortIndexList(AllIndex)
#                 else:
#                     AllIndex.append(NewPair)
#             if CurrIndex <= 0:
#                 CurrIndex = len(Sequence) - 1


# def SortIndexList(IndexList):
#     AllBases = sorted(CleanList(IndexList))
#     FinalIndexList = []
#     for i in AllBases:
#         NextIndex = getIndexOfElement(IndexList,i)[0]
#         if NextIndex not in FinalIndexList:
#             FinalIndexList.append(NextIndex)
#     return([IndexList[x] for x in FinalIndexList])


# def getIndexOfElement(PairMap,Element):
#
#     dep = 0
#     NumLevels = depth(PairMap)
#     CurrPairMap = list(PairMap)
#     NumElements = len(PairMap)
#     IndexList = []
#     Element = FlattenList(Element)[0]
#     # for x in BPMap:
#     #     if Element in FlattenList(x):
#     #         print("found it",x, Element,BPMap.index(x), FlattenList([BPMap.index(x) for x in BPMap if Element in FlattenList(x)]))
#     IndexList.append(FlattenList([PairMap.index(x) for x in PairMap if Element in FlattenList(x)]))
#     IndexList = FlattenList(IndexList)
#     NumLevels -=1
#     CurrPairMap = list(PairMap)
#     if IndexList==[]:
#         return IndexList
#     #print("iondex list here", IndexList,Element,BPMap, IndexList[-1],CurrPairMap[0])
#     CurrPairMap = CurrPairMap[IndexList[-1]]
#     while NumLevels>1:
#         IndexList.append(FlattenList([PairMap.index(x) for x in PairMap if Element in FlattenList(x)]))
#         IndexList = FlattenList(IndexList)
#         NumLevels -= 1
#         CurrPairMap = list(PairMap)
#         CurrPairMap = CurrPairMap[IndexList[-1]]
#
#     return IndexList







# def getBulgeIDFromBase(BaseID,Sequence):
#     BulgedBases, AllBasesInBulge = getBulgeListOfBases(Sequence)
#     if type(BaseID)!=list:
#         return [AllBasesInBulge.index(x) for x in AllBasesInBulge if BaseID in x]
#     else:
#         RetIDs = []
#         for b in BaseID:
#             RetIDs += [AllBasesInBulge.index(x) for x in AllBasesInBulge if b in x]
#         return RetIDs



##### This is where I will be inserting new functions


# def getFGInfo(GroupFeatures):
#
#     ListOfFGTypes, HPsInFG, ListOfNeighborFGs, NumFGs = GroupFeatures
#     ListOfChildren = []
#     ListOfSiblings = []
#     ListOfParents = []
#     ListOfBulgeIndices = []
#     ListisLeaf = []
#     ListIsJxn = []
#     for i in HPsInFG: #Iterating through each FG
#         ParentFG= []
#         SiblingFG = []
#         ChildFG = []
#         if RNAStruct.HPSt[i[-1]].isLeaf=='yes':
#             ListisLeaf+= ['yes']
#         else:
#             ListisLeaf+= ['no']
#
#         FirstHPInFG = i[0]
#         LastHPInFG = i[-1]
#         ParentHP = RNAStruct.HPSt[i[0]].Parent
#         SiblingHPs = RNAStruct.HPSt[i[0]].Siblings
#         ChildHPs = RNAStruct.HPSt[i[-1]].Children
#
#         if len(SiblingHPs) > 1 or len(ChildHPs) > 1:
#             ListIsJxn += ['yes']
#         else:
#             ListIsJxn += ['no']
#
#         if ParentHP != []: #Finding Parent FG
#             ParentFG += [HPsInFG.index(k) for k in HPsInFG if ParentHP[0] in k]
#         if SiblingHPs != []:
#             for j in SiblingHPs:
#                 SiblingFG += [HPsInFG.index(k) for k in HPsInFG if j in k]
#         if ChildHPs != []:
#             for j in ChildHPs:
#                 ChildFG += [HPsInFG.index(k) for k in HPsInFG if j in k]
#         ListOfChildren.append(ChildFG)
#         ListOfSiblings.append(SiblingFG)
#         ListOfParents.append(ParentFG)
#         ListOfBulgeIndices.append([])
#     # print("LOC",ListOfChildren)
#     # print("LOS",ListOfSiblings)
#     # print("LOP",ListOfParents)
#     # print("LOBI",ListOfBulgeIndices)
#
#     return ListOfChildren, ListOfSiblings, ListOfParents, ListOfBulgeIndices,ListisLeaf,ListIsJxn
#
# def getFGIDFromHPList(HPList):
#
#
#     FGIDList = [None] *len(HPList)
#
#     try:
#         HPList[0]
#     except:
#         HPList = [HPList]
#
#     for i in HPList:
#         try:
#             i = i[0]
#         except:
#             print("you messed up?", i)
#
#         FGNo = HPList.index(i)
#         FGIDList[FGNo] = i
#         print("hahah", i, FGNo,FGIDList,HPList)
#
#     return FGIDList
#
# def groupFeaturesByPath(FGSt,HPSt):
#     """ PRETTY MUCH DONE I JUST NEED TO MAKE IT TRANSLATE BETWEEN HAIRPINS AND FGS"""
#
#
#     print("Beginning the grouping of paths")
#     NumFeatures = len(FGSt)
#     AllFGs = []
#     NewGroupFlag = 1
#     GroupedFGs = []
#
#     for CurrFG in range(0,NumFeatures):
#         CurrHPs = FGSt[CurrFG].HPList
#         print(CurrHPs)
#
#         if NewGroupFlag ==1:
#             try:
#                 AllFGs += [newHPs]
#             except:
#                 newFGs = []
#             newHPs = [CurrHPs]
#         print("-hpST ID->", [HPSt[CurrHPs[-1]].ID])
#         HPSt[CurrHPs[-1]].printme()
#         print("-ChildrenID->", [HPSt[CurrHPs[-1]]])
#
#         ChildFGs = getFGIDFromHPList([HPSt[CurrHPs[-1]].ID])
#         print(ChildFGs)
#
#         #print(FG[CurrFG].Children)
#
#     #
#     # for CurrFG in range(0,NumFeatures): #Iterate through each FG
#     #     if NewGroupFlag ==1:
#     #         try:
#     #             AllFGs += [newFGs]
#     #         except:
#     #             newFGs = []
#     #         newFGs = [CurrFG]
#     #
#     #
#     #     ChildFGs = FG[CurrFG].Children[0].HP[0].FGID
#     #     print(ChildFGs)
#     #     CurrFG = CurrFG+1
#     #     print(ChildFGs)
#     #
#     #     print(AllFGs)
#     #     print("Watttt",CurrFG,newFGs,ChildFGs)
#     #     if FG[CurrFG].isLeaf=='no' and len(ChildFGs)==1:
#     #         newFGs += ChildFGs
#     #         CurrFG = ChildFGs[0]
#     #         NewGroupFlag = 0
#     #     else:
#     #         NewGroupFlag = 1
#
#     if NewGroupFlag==1:
#         AllFGs += [newFGs]
#     print(AllFGs)
#
#     return AllFGs



# #Functions for Feature Group Definitions
# def groupHairpinsByFeature(FG):
#     NumHPs = len(FG.HP)
#     IncludedHPs = []
#     NewGroupFlag = 1
#     GroupedHPs = []
#
#     for CurrHP in range(0,NumHPs):
#         if NewGroupFlag ==1:
#             try:
#                 GroupedHPs += [NewGroup]
#             except:
#                 GroupedHPs = []
#             NewGroup = [CurrHP]
#             IncludedHPs += [CurrHP]
#             NewGroupFlag = 0
#         NextHP = FG.HP[CurrHP].Children
#
#         if FG.HP[CurrHP].isLeaf=='no' and len(NextHP)==1:
#             if FG.HP[CurrHP].FeatureType==FG.HP[NextHP[0]].FeatureType:
#                 NewGroup += NextHP
#                 IncludedHPs += NextHP
#                 NewGroupFlag = 0
#             else:
#                 NewGroupFlag = 1
#         else:
#             NewGroupFlag = 1
#
#     if NewGroupFlag==1:
#         GroupedHPs += [NewGroup]
#     print(IncludedHPs)
#
#     return GroupedHPs





# def getBulgeInfo(Sequence,BPMap):
#
#     # Defining length of puzzle, and pair list describing all bulges in the puzzle
#     PuzzleLength = len(Sequence) - 1
#     sortedPairMap = list(sorted(FlattenList([0, BPMap, PuzzleLength])))
#     sortedPairMap = getPairsFromList(list(sorted(FlattenList([0, BPMap, PuzzleLength]))))
#     NumBulges = len(sortedPairMap)
#
#     BulgeList = []
#     BulgeTypeList = []
#     BulgeLengthList = []
#     RemainingBulges = list(sortedPairMap)
#
#     # Iterating through each bulge
#     for i in range(0, NumBulges):
#         CurrBulge = sortedPairMap[i]
#         flag = 1
#         tempBulgeList = [CurrBulge]
#         if CurrBulge in RemainingBulges:  # If we have not seen this one before
#             RemainingBulges.remove(CurrBulge)
#             while flag == 1:
#                 ConnectedBulge = getIndexPartner(BPMap, CurrBulge[1])  # Find the connected bulge
#                 if ConnectedBulge == -1:  # If we are at the end of the puzzle,
#                     flag = 0
#                     NewBulge = []
#                 else:
#                     NewBulge = getPairFromSingleIndex(sortedPairMap, ConnectedBulge, 0)[
#                         0]  # and get connected bulge's indexes
#
#                 if NewBulge == []:  # Make sure we are not at the end of the puzzle
#                     flag = 0
#                 elif NewBulge in RemainingBulges:
#                     tempBulgeList += [NewBulge]  # Add the new bulge to the list
#                     RemainingBulges.remove(NewBulge)
#                     CurrBulge = NewBulge
#                 else:
#                     # print(i,sortedPairMap[i],ConnectedBulge,NewBulge, tempBulgeList)
#                     flag = 0
#             # Add the set of new bulges to the list
#             BulgeList += [tempBulgeList]
#
#             # Get the length of each bulge
#             tempBulgeLengthList = [x[1] - x[0] - 1 for x in tempBulgeList]
#             if tempBulgeList[0][0] == 0 and tempBulgeList[-1][1] == PuzzleLength:  # Correcting outermost bulge
#                 tempBulgeLengthList[0] += 1
#                 tempBulgeLengthList[-1] += 1
#             BulgeLengthList += [tempBulgeLengthList]
#
#             # Gets the type of each bulge
#             if len(tempBulgeList) == 1:
#                 BulgeTypeList += ['leaf']
#             elif len(tempBulgeList) == 2:
#                 BulgeTypeList += ['branch']
#             else:
#                 BulgeTypeList += ['fork']
#     if len(RemainingBulges) > 0:
#         print("Warning: Some bulges were not classified!")
#     return BulgeList, BulgeTypeList, BulgeLengthList


# def getPairFromSingleIndex(BPMap, Index, IndexPosition):
#     return [x for x in BPMap if Index == x[IndexPosition]]

# def getIndexPartner(BPMap, Index):
#     for i in BPMap:
#         FirstPair = i[0]
#         SecondPair = i[1]
#         if FirstPair[0] <= Index <= FirstPair[1]:
#             Diff = Index - FirstPair[0]
#             return SecondPair[1] - Diff
#         if SecondPair[0] <= Index <= SecondPair[1]:
#             Diff = Index - SecondPair[0]
#             return FirstPair[1] - Diff
#     return -1

# def getFamilyList(BPMap,Index):
#
#     Item = FlattenList(getPairFromSingleIndex(BPMap,Index,BPMap[Index]))





##
# def getIndexList(SS):
#     # Requires Variables:
#     # SS
#     # Requires Functions:
#     # getPartnerIdxs
#     List_Of_Partner_Idx, List_Of_Pairs = getPartnerIdxs(SS)
#     tmp_List_Of_Pairs = List_Of_Pairs
#     List_Of_HPs = []
#     List_Of_PairMaps = []
#     i = 0
#     while i < len(tmp_List_Of_Pairs):  # Identifies pairs that are contiguous
#         tmp_List_Of_HPs = [tmp_List_Of_Pairs[i]]
#         HP_Len_Count = 1  # counting how long the hairpin is
#         New_HP_Flag = 0
#         while New_HP_Flag == 0 and i + 1 < len(tmp_List_Of_Pairs):
#             New_HP_Flag = 1
#             if tmp_List_Of_Pairs[i][0] == tmp_List_Of_Pairs[i + 1][0] - HP_Len_Count and tmp_List_Of_Pairs[i][1] == \
#                     tmp_List_Of_Pairs[i + 1][1] + HP_Len_Count:
#                 tmp_List_Of_HPs.append(tmp_List_Of_Pairs[i + 1])
#                 tmp_List_Of_Pairs.remove(tmp_List_Of_Pairs[i + 1])
#                 HP_Len_Count = HP_Len_Count + 1
#                 New_HP_Flag = 0
#         List_Of_HPs.append([tmp_List_Of_HPs[0], tmp_List_Of_HPs[-1]])
#         i = i + 1
#     for i in range(0, len(List_Of_HPs)):  # Re-organizes to be separated by strand
#         List_Of_PairMaps.append(
#             [[List_Of_HPs[i][0][0], List_Of_HPs[i][1][0]], [List_Of_HPs[i][1][1], List_Of_HPs[i][0][1]]])
#     return List_Of_PairMaps  # List of all List_Of_Partner_Idxs, list of all pairs, list of pairs sorted by HP
#
# def getIndexListFromHPList(HPList, ParentFlag, IncludeParentSiblings):
#
#     #   Step 1: Find all HPs surrounding:
#
#     HPList = CleanList(HPList)
#     HPList.append(getAllHPsConnectingList(HPList,ParentFlag)) #Connects all HPs and includes them in the structure
#     HPList = CleanList(HPList)
#
#     MCH = getMinCommonHairpin(HPList)
#     IndexOfChildren = getSiblingsAndChildren(HPList)
#     ParentHP = 'test'
#     #print(HPList,MCH, "MCH", RNAStruct.HPSt[MCH].Parent, "Flags:", ParentFlag, IncludeParentSiblings)
#     #   Find min common hairpin
#     #print("Containing Start", HPList)
#
#     if ParentFlag == 1 and MCH in HPList: #if we want the parent and the MCH is already included, add the parent of the MCH
#         if RNAStruct.HPSt[MCH].Parent != []: #If the parent is a hairpin (i.e. we are not at opening HPs of structure)
#             MCH = RNAStruct.HPSt[MCH].Parent[0] #New MCH is the parent of the old MCH
#             HPList.append(MCH) #adding MCH to the HPList
#             IndexOfChildren = getSiblingsAndChildren(HPList)
#             OpeningHPs = sorted([MCH] + RNAStruct.HPSt[MCH].Siblings)
#         else:
#             OpeningHPs = sorted([x for x in RNAStruct.OpeningHPs if x in HPList + IndexOfChildren])
#     elif ParentFlag == 1 and MCH not in HPList: #If the MCH is not in the HPList
#         HPList.append(MCH)
#         IndexOfChildren = getSiblingsAndChildren(HPList)
#         OpeningHPs = sorted([MCH] + RNAStruct.HPSt[MCH].Siblings)
#     elif ParentFlag ==0 and MCH in HPList:
#         OpeningHPs = sorted([MCH] + RNAStruct.HPSt[MCH].Siblings)
#     elif ParentFlag == 0 and MCH not in HPList:
#         OpeningHPs = sorted(RNAStruct.HPSt[MCH].Children)
#
#     HPList = CleanList(HPList)
#     #print("Containing Midway", HPList, MCH, ParentHP)
#
#     #Remove parent siblings if desired
#     if IncludeParentSiblings ==0:
#         SiblingsToRemove = [x for x in RNAStruct.HPSt[MCH].Siblings if x in HPList] #
#         for i in SiblingsToRemove:
#             HPList.pop(HPList.index(i))
#         SiblingsToRemove = [x for x in RNAStruct.HPSt[MCH].Siblings if x in OpeningHPs]
#         for i in SiblingsToRemove:
#             OpeningHPs.pop(OpeningHPs.index(i))
#         SiblingsToRemove = [x for x in RNAStruct.HPSt[MCH].Siblings if x in IndexOfChildren]
#         for i in SiblingsToRemove:
#             IndexOfChildren.pop(IndexOfChildren.index(i))
#
#     HPList = CleanList(HPList)
#     #print("Containing Final", HPList, IndexOfChildren)
#     FirstHP = OpeningHPs[0]
#     LastHP = OpeningHPs[-1]
#     FirstIndex = getNeighborIndicesOffset(FirstHP, 0)[0]
#     LastIndex = getNeighborIndicesOffset(LastHP, 0)[1]
#     #print("OpeningHPs", OpeningHPs, FirstIndex,LastIndex)
#
#
#     IndexList = [FirstIndex]
#
#     for i in sorted(IndexOfChildren):
#         IdxA = getInsideIndex(i)[0]-1
#         IdxB = getInsideIndex(i)[1]+1
#         IndexList.append([IdxA,IdxB])
#     IndexList.append(LastIndex)
#     #print("IL here", HPList, IndexOfChildren, IndexList)
#     #print(HPList, IndexOfChildren, IndexList, ParentFlag)
#     IndexList = FlattenList(IndexList)
#     #IndexList = sorted(IndexList)
#     FinalList = getPairsFromList(IndexList)
#     #print("Final List", FinalList, HPList, OpeningHPs, FirstIndex,LastIndex)
#
#     return FinalList
#
#
#
# def getCompleteFamilyList(BPMap):  # Inputs SS, Outputs list of isLeaf, Children, Parents, Siblings
#     ##Added
#     # Requires Variables:
#     # BPMap (global)
#     # Requires functions:
#     # getBlankFamilyList
#     # getIsLeaf (getNextStemIdx)
#     # getAllChildren (getNextStemIdx)
#     # getOpeningHPs (getNextStemIdx)
#     # {  Initializing list of empty lists of length # hairpins in SS
#
#     # Dependencies
#     # getNextStemIdx <- getStemIdxFromBaseIdx <- getNextPairedIndex
#     Num_Hairpins = len(BPMap)
#
#     # Step 2a: Get blank FamilyList:
#     List_isLeaf, List_Of_Children, List_Of_Parents, List_Of_Siblings = getBlankFamilyList(Num_Hairpins)
#
#     # Step 2b: Completes Family List
#     for i in range(0, Num_Hairpins):
#
#         # Step 2c: Find all hairpins that are leafs
#         List_isLeaf[i] = getIsLeaf(i)
#         if List_isLeaf[i] == 'yes':
#             continue
#
#         # Step 2d: Find all hairpins that have children
#         List_Of_Children[i] = getAllChildren(i)
#
#         # Step 2e: Assign parents, siblings to all hairpins
#         if List_Of_Children[i] != None:
#             # {  Iterates through each child and assigns its parent and siblings
#             for j in List_Of_Children[i]:
#                 # {  Removes the hairpin of interest from the list of siblings
#                 temp_Siblings = list(List_Of_Children[i])
#                 temp_Siblings.pop(temp_Siblings.index(j))  # removes itself from sibling list
#                 # }
#                 # {  Tests to see if we have assigned a parent before and assigns if we have not
#                 if List_Of_Parents[j] == []:
#                     List_Of_Parents[j] = [i]
#                 # }
#                 # {  Tests to see 1) if have assigned siblings before and 2) whether there are siblings to assign
#                 if List_Of_Siblings[j] == [] and len(List_Of_Children[i]) > 1:
#                     List_Of_Siblings[j] = temp_Siblings
#                 # }
#             # }
#
#     # Step 2f: Assign opening hairpins as siblings
#     OpeningHPs = getOpeningHPs()
#     for i in OpeningHPs:
#         temp_Siblings = list(OpeningHPs)
#         temp_Siblings.pop(temp_Siblings.index(i))
#         List_Of_Siblings[i] = temp_Siblings
#
#     return List_isLeaf, List_Of_Children, List_Of_Parents, List_Of_Siblings
# def main():
#     Sequence = FD.SequenceStruct()
#     StartIdx = 399
#     #getAllSurroundingBPs(FD.Sequence, StartIdx, 0)
#     #getAllSurroundingBPs(FD.Sequence, StartIdx, 1)
#     print("new")
#     print(getAllSurroundingBPs4(FD.Sequence, StartIdx, 0))
#     print(getAllSurroundingBPs4(FD.Sequence, StartIdx, 1))
#     print("old:")
#     print(getAllSurroundingBPs(FD.Sequence, StartIdx, 0))
#     print(getAllSurroundingBPs(FD.Sequence, StartIdx, 1))



#main()



