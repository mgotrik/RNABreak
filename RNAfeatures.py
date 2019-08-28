import input
import basicfunc as bf
global Sequence

class SequenceStruct():

    def __init__(self):

        #Defines the sequences
        self.SeqList = {1: input.PS_Goal, 2: input.SS_Goal, 3: input.LS_Goal, 4: input.IUPAC_Goal}
        #Goes through SS and identifies every bases partner (or -1 for bulge)
        self.PartnerIdxs = self.getPartnerIdxs()[0]
        #Creates a list of all base pairs [[0, 360], [10, 40],  ... ]
        #Also gives the immediate 5p and 3p bulges of every bp
        self.BPMap, self.BPBulge5p, self.BPBulge3p = self._getBPBulgeMap()
        #Creates a list of all HPs, consisting of lists of their base pairs [[[0, 360], [1,359], [2,358]],  ... ]
        #Also gives the immediate 5p and 3p bulges of every HP
        self.HPMap, self.HPBulge5p, self.HPBulge3p = self._getHPBulgeMap()

        #Identifies all neighbors (0 = upstream / towards start; 1 = downstream) of every base pair
        self.BPSurrBP = [[self._getSurroundingBPs(i,0), self._getSurroundingBPs(i,1)] for i in self.BPMap]
        #Identifies all neighbors (0 = upstream / towards start; 1 = downstream) of every hairpin
        self.HPSurrBP = [[self._getSurroundingBPs(i[0],0), self._getSurroundingBPs(i[-1],1)] for i in self.HPMap]

        #Identifies the indices of all bulges, includes [] to signify a 0-length bulge
        #Crawls from base 0 to end, and lists bulges in order of encountering. All bulges are given in order
        self.BulgeMap, self.BulgeParentHPMap = self._getBulgeMap()

        #List of the parent BP of each BP
        self.BPParentMap = self._getBPParentMap()
        #List of the parent HP of each HP
        self.HPParentMap = self._getHPParentMap()

        #List of the child/sibling HP of each HP. Same function could be used for BPParentMap, but isn't needed
        self.HPChildMap = self._getChildMapFromParentMap(self.HPParentMap)
        self.HPSibMap = self._getSiblingMapFromParentMap(self.HPParentMap)

        #Gets information regarding the relative order of bulges. Every HP is defined as having a parent and child bulge, and these maps list it.
        self.BulgeChildHPMap = [bf.FlattenList([self.HPChildMap[x] for x in self.BulgeParentHPMap[y]]) for y in range(0,len(self.BulgeParentHPMap))]
        self.BulgeChildHPMap[0] = self.HPSibMap[0]+[0] #Confirms thats the parent HP and its siblings are assigned the same opening bulge
        self.BulgeParentMap = self._getBulgeParentMap(self.HPParentMap,self.BulgeParentHPMap)
        # for i in range(0,len(self.BulgeParentMap)):
        #     print(i, self.BulgeParentHPMap[i])
        #     print(i, self.BulgeChildHPMap[i])

    #Lets you call Sequence[2] to get the secondary structure
    def __getitem__(self, item):
        try:
            return self.SeqList[item]
        except:
            print('Invalid SequenceStruct Specified', self, item)
            return []
    #len(Sequence) = length of the whole sequence
    def __len__(self):
        return (len(input.SS_Goal))
    # printing sequence gives you the actual sequence
    def __str__(self, item=1):
        return self.SeqList[item]
    def __repr__(self):
        return "SequenceStruct()"

    #These functions define the 'pair maps' i.e. lists of grouped indices. This makes up the core of the puzzle feature definition logic (i.e. defining a HP vs. base pair)
    def getPartnerIdxs(self):  # Input: '...((((...)).))..'; Output: [-1, -1, -1, 14, 13, 11, 10, -1, -1, -1, 6, 5, -1, 4, 3, -1, -1], [[3, 14], [4, 13], [5, 11], [6, 10]]
        # Req'd Variables:
        # SS
        # Req'd Functions:
        # None
        SS = self.SeqList[2]
        List_Of_Pairs = []  # List of all base pair partnerslist
        List_Of_Partner_Idx = [-1] * len(SS)  # List showing which base is paired with which
        SS_Base_Count = 0  # Which base we are looking at
        SS_Open_Idx = []  # Index of last open pair
        SS_Closed_Idx = []  # Index of last closing pair
        for i in list(SS):  # Identifies List_Of_Partner_Idx of each base
            if i == '(':
                SS_Open_Idx.append(SS_Base_Count)
                Open_Pair_Flag = 0
            elif i == ')':
                SS_Closed_Idx.append(SS_Base_Count)
                Open_Pair_Flag = 1
            elif i == '.':
                Open_Pair_Flag = 0

            if Open_Pair_Flag == 1:
                List_Of_Partner_Idx[SS_Open_Idx[-1]] = SS_Base_Count
                List_Of_Partner_Idx[SS_Base_Count] = SS_Open_Idx[-1]
                List_Of_Pairs.append([SS_Open_Idx[-1], SS_Closed_Idx[-1]])
                SS_Open_Idx.remove(SS_Open_Idx[-1])
                SS_Closed_Idx.remove(SS_Closed_Idx[-1])
            SS_Base_Count = SS_Base_Count + 1
        List_Of_Pairs = sorted(List_Of_Pairs)
        return List_Of_Partner_Idx, List_Of_Pairs
    def _getBPBulgeMap(self):

        #Get flattened Pairmap list:
        PartnerIdxs = list(self.PartnerIdxs)
        AllBasePairs = []
        AllBPList = []
        AllBPIdx5pList = []
        AllBPIdx3pList = []
        BPLevelList = []
        BPCount = 0
        BPLevel = 0
        BulgeLen5p = 0
        BulgeLen3p = 0
        BulgeLenList5p = []
        BulgeLenList3p = []
        BulgeIndexList3p = []
        BulgeIndexList5p = []
        BulgeList3p = []
        BulgeList5p = []

        for i in range(0,len(PartnerIdxs)):
            if PartnerIdxs[i]==-1: #if index is a bulge
                BulgeList5p+=[i]
                BulgeList3p+=[i]
                BulgeLen5p+=1
                BulgeLen3p+=1
                continue
            else:
                BasePair = sorted([i,PartnerIdxs[i]])
                if BasePair in AllBasePairs: #If this is the second time we've seen this base pair
                    #AllBPIndex = AllBasePairs.index(BasePair)
                    AllBasePairs.append(BasePair)
                    BPLevelList.append(BPLevel)
                    BPLevel -=1
                else:  #If this is the first time we've seen this base pair
                    BPCount +=1
                    BPLevel+=1
                    AllBasePairs.append(BasePair)
                    BPLevelList.append(BPLevel)
                BulgeLenList5p.append(BulgeLen5p)
                BulgeIndexList5p.append(BulgeList5p)
                if len(BulgeLenList5p)>1:
                    BulgeLenList3p.append(BulgeLen3p)
                    BulgeIndexList3p.append(BulgeList3p)
                BulgeLen5p = 0
                BulgeLen3p = 0
                BulgeList3p = []
                BulgeList5p = []
        BulgeLenList3p.append(BulgeLen3p)
        BulgeIndexList3p.append(BulgeList3p)
        # for i in range(0,len(AllBasePairs)):
        #     print(i, AllBasePairs[i],BPLevelList[i],BulgeLenList5p[i],BulgeLenList3p[i], BulgeIndexList5p[i],BulgeIndexList3p[i])

        #Cleans up the previous list by appending the 5` and 3` bulges

        BPCount = 0
        for i in range(0,len(AllBasePairs)-1):
            CurrBP = AllBasePairs[i]
            if CurrBP!=[]:
                for j in range(i+1,len(AllBasePairs)):
                    if AllBasePairs[j]==CurrBP:
                        AllBPList.append(CurrBP)
                        AllBPIdx5pList.append([BulgeIndexList5p[i]] + [BulgeIndexList3p[j]])
                        AllBPIdx3pList.append([BulgeIndexList3p[i]] + [BulgeIndexList5p[j]])
                        AllBasePairs[j] = []
                        FoundPairFlag = 1
                        BPCount+=1
                        break

        return AllBPList,AllBPIdx5pList, AllBPIdx3pList
    def _getHPBulgeMap(self):

        BPMap = list(self.BPMap)
        HPMap = []
        HPList = []
        BPIDList = list(range(0,len(BPMap)))
        BPBulgeList5p = list(self.BPBulge5p)
        BPBulgeList3p = list(self.BPBulge3p)

        # BPID = 108
        # print("--->",BPBulgeList3p[BPID])
        # print(BPBulgeList5p[BPID])
        breakflag = 0
        #Iterates through each base pair
        for i in BPIDList:
            CurrBP = BPMap[i]
            #
            if CurrBP!=[]: #If we have not already seen this BP before
                HPList = [CurrBP]
                BPMap[i] = []
                if bf.FlattenList(BPBulgeList3p[i]) == []: #If this HP is not length 1 (if it is, will have bulges in 3p)
                    for j in range(i + 1, len(BPMap)): #Iterating through the BPs until we complete the HP
                        NextBP = BPMap[j]
                        HPList.append(NextBP)
                        BPMap[j] = []
                        if bf.FlattenList(BPBulgeList3p[j])!=[]:
                            break
                HPMap.append(HPList)

        BPMap = list(self.BPMap)
        HPBulgeList5p = []
        HPBulgeList3p = []
        for i in range(0,len(HPMap)):
            FirstIdx = BPMap.index(HPMap[i][0])
            LastIdx = BPMap.index(HPMap[i][-1])
            HPBulgeList5p.append(BPBulgeList5p[FirstIdx])
            HPBulgeList3p.append(BPBulgeList3p[LastIdx])
            #print(i, HPMap[i], "Bulges on 5p end", HPBulgeList5p[i], "Bulges on 3p end", HPBulgeList3p[i])

        return HPMap, HPBulgeList5p, HPBulgeList3p
    def _getSurroundingBPs(self,Index,Direction_Flag):
        PartnerIdxs = list(self.PartnerIdxs)
        if Direction_Flag==1:
            Step_Iter = 1
        else:
            Step_Iter = -1

        StartIndex = min(Index)
        if PartnerIdxs[StartIndex]==-1:
            BPList = []
            CurrIndex = StartIndex
        else:
            StartPair = [StartIndex,PartnerIdxs[StartIndex]]
            BPList = [StartPair]
            CurrIndex = StartIndex + Step_Iter


        while CurrIndex!=StartIndex:
            try:
                if PartnerIdxs[CurrIndex]==-1:
                    CurrIndex += Step_Iter
                else:
                    if CurrIndex<=-1:
                        CurrIndex = len(self)-abs(CurrIndex)
                    elif CurrIndex>=len(self):
                        CurrIndex = 0 + abs(CurrIndex-len(self))-1
                    NewPair = sorted([CurrIndex,PartnerIdxs[CurrIndex]])
                    if NewPair in BPList:
                        break
                    else:
                        BPList.append(NewPair)
                    CurrIndex = PartnerIdxs[CurrIndex] + Step_Iter
            except:
                if Direction_Flag==1:
                    CurrIndex = 0
                else:
                    CurrIndex = len(self)-1
        try:
            BPList.remove(StartPair)
        except:
            pass
        return BPList
    def _getBulgeMap(self):

        ##Get Bulges around first HP:
        BulgeMap = []

        SurroundingBPs = self._getSurroundingBPs([0],0)
        SurroundingBPIDs = sorted([self.BPMap.index(x) for x in SurroundingBPs])
        if 0 not in SurroundingBPIDs:
            SurroundingBPIDs+=[0]
        Bulges = [self.BPBulge5p[0][0]]
        for i in SurroundingBPIDs:
            Bulges += [self.BPBulge5p[i][1]]

        BulgeParent = [[]]
        BulgeMap.append(Bulges)

        #Get bulges at the end of every HP
        for HPID in range(0,len(self.HPMap)):
            Index = self.HPMap[HPID][-1]
            SurroundingBPs = self._getSurroundingBPs(Index,1)
            #print("->",SurroundingBPs,Index)
            SurroundingBPIDs = sorted([self.BPMap.index(x) for x in SurroundingBPs])
            ParentBP = self.BPMap.index(Index)
            if SurroundingBPIDs ==[]:
                Bulges = [self.BPBulge3p[ParentBP][0]]
            else:
                Bulges = [self.BPBulge5p[SurroundingBPIDs[0]][0]]
                for i in SurroundingBPIDs:
                    Bulges += [self.BPBulge5p[i][1]]
            BulgeParent.append([HPID])
            BulgeMap.append(Bulges)

        return BulgeMap, BulgeParent

    #These functions define the familial relation of BP/HPs. A parent is the nearest upstream neighbor, children are downstream neighbors, and siblings share the same parent
    def _getBPParentMap(self):

        BPSurrBP = list(self.BPSurrBP)
        ParentSibs = BPSurrBP[0][0]
        ParentHP = [0]
        for Sib in ParentSibs:
            ParentHP += [self.BPMap.index(Sib)]

        #print(ParentHP, self.BPMap[0],BPSurrBP[0])
        ParentList = []

        CurrID = 0
        for i in BPSurrBP:
            #print(i)
            BPIDs = [self.BPMap.index(x) for x in i[0]]
            #print(BPIDs)
            if CurrID not in ParentHP:
                ParentList.append([min(BPIDs)])
            else:
                ParentList.append([])
            CurrID +=1
        return ParentList
    def _getHPParentMap(self):
        HPSurrBP = list(self.HPSurrBP)
        ParentSibs = HPSurrBP[0][0]
        ParentHP = [0]
        for Sib in ParentSibs:
            ParentHP += [self.HPMap.index(x) for x in self.HPMap if Sib in x]

        #print(ParentHP, self.HPMap[0],HPSurrBP[0])
        ParentList = []

        CurrID = 0
        for i in HPSurrBP:
            #print(i)
            #HPIDs = [self.HPMap.index(x) for x in i[0]]
            HPIDs = [self.HPMap.index(x) for x in self.HPMap if any(elem in x for elem in i[0])]
            #print(HPIDs)
            if CurrID not in ParentHP:
                ParentList.append([min(HPIDs)])
            else:
                ParentList.append([])
            CurrID +=1
        return ParentList
    def _getChildMapFromParentMap(self, ParentMap):

        #ChildMap = list([[]]*len(ParentMap))
        ChildMap = [[] for i in range(0,len(ParentMap))]
        for i in range(0,len(ParentMap)):
            ParentID = bf.FlattenList(ParentMap[i])
            ChildID = [i]
            if ParentID!=[]:
                #print(ParentID, ChildID)
                try:
                    ChildMap[ParentID[0]]+=ChildID
                except:
                    ChildMap[ParentID]+=ChildID
        return ChildMap
    def _getSiblingMapFromParentMap(self,ParentMap):

        SiblingMap = [[] for i in range(0,len(ParentMap))]

        #Identify parent siblings:
        count = 0
        for i in ParentMap:
            if i==[] and count!=0:
                SiblingMap[0]+=[count]
                SiblingMap[count]+=[0]
            count+=1

        ChildMap = self._getChildMapFromParentMap(ParentMap)

        for i in range(0,len(ChildMap)):
            ChildIDs = ChildMap[i]
            if len(ChildIDs)>1:
                #print("test",ChildIDs)
                for j in ChildIDs:
                    Siblings = list(ChildIDs)
                    Siblings.remove(j)
                    #print(j, Siblings, SiblingMap[j])
                    SiblingMap[j] += Siblings

        return SiblingMap
    def _getBulgeParentMap(self,HPParentMap,BulgeParentHPMap):

        BPHPMap = list(BulgeParentHPMap)
        PHPMap = list(HPParentMap)
        ParentList = [[] for i in range(0,len(BulgeParentHPMap))]

        for i in range(0,len(BulgeParentHPMap)):
            ParentHP = BPHPMap[i]
            if ParentHP==[] or ParentHP==[[]]:
                continue
            GrandParentHP = PHPMap[ParentHP[0]]
            if GrandParentHP==[] or GrandParentHP==[[]]:
                ParentList[i] += [BulgeParentHPMap.index(x) for x in BulgeParentHPMap if x==GrandParentHP]
            else:
                ParentList[i] += [BulgeParentHPMap.index(x) for x in BulgeParentHPMap if x==GrandParentHP]

        return ParentList

Sequence = SequenceStruct()

#This is the main class that is taken from this file. Allows you to reference all features in a structure as FG.HP[HPID], or FG.Bulge._IDlist, etc.
class FeatureGroup():
    def __init__(self):
        #Sequence = SequenceStruct()
        global BaseMap, BPMap, HPMap, BulgeMap, HP, PathMap, Bulge, HPGroup

        self.FeatureList = {'Base','BP','HP','Bulge','Path','HPGroup'}
        self.Sequence = SequenceStruct()
        self.Base = BaseStruct(Sequence)
        BaseMap = self.Base._pairmap
        self.NumBases = len(self.Base)
        # print("Number of Bases:", len(self.Base), self.Base.length)
        # print([key for key in self.Base.__dict__])
        #print("Number of Bases:", len(self.Base))
        #self.Base._printme()

        # Gets the BasePair (BP) BPMap -> [1, 352] signifying 1-bp-352
        self.BP = BasePairStruct(Sequence)


        self.NumBPs = len(self.BP)
        # print("Number of BP:", len(self.BP), self.BP.length)
        # print([key for key in self.BP.__dict__])

        # Gets the Hairpin (HP) BPMap -> [[1, 14],[2,13]] signifying a 2bp Hairpin
        self.HP = HairpinStruct(Sequence)
        HPMap = self.HP._pairmap
        HP = self.HP
        self.NumHPs = len(self.HP)
        # print("Number of HP:", len(self.HP), self.HP.length)
        # self.HP._printme(0)
        # print([key for key in self.HP.__dict__])

        # Gets the Bulge (HP) Map -> [[1, 14],[2,13]] signifying a 2bp Hairpin
        self.Bulge = BulgeStruct(Sequence)
        BulgeMap = self.Bulge._pairmap
        self.NumBulges = len(self.Bulge)
        Bulge = self.Bulge
        # print("Number of Bulges:", len(self.Bulge), self.Bulge.length)
        # self.Bulge._printme(0)
        # print([key for key in self.Bulge.__dict__])
        #print([key for key in self.Bulge.__dict__])


        #   Searches all HPs for HPs that are of the same type ('kind'), none specifically (None), that are continuous ('near'), and under no ID list constraint (None)
        self.HPGroup = HPGroupStruct(Sequence)
        HPGroupMap = self.HPGroup._pairmap
        self.NumHPGroups = len(HPGroupMap)
        HPGroup = self.HPGroup
        # print("Number of HPGroups:", len(self.HPGroup), self.HPGroup.length, len(self.HPGroup._pairmap))
        # self.HPGroup._printme(0)
        # print([key for key in self.HPGroup.__dict__])
        #print([key for key in self.HPGroup.__dict__])

        # Gets the Path (Path) Map -> [[1, 14],[2,13]] signifying a 2bp Hairpin
        self.Path = PathStruct(Sequence)
        PathMap = self.Path._pairmap
        self.NumPaths = len(self.Path)
        # print("Number of Paths:", len(self.Path), self.Path.length)
        # self.Path._printme(0)
        # print([key for key in self.Path.__dict__])
        #print([key for key in self.Path.__dict__])
        #
        # for i in self.Bulge._IDlist:
        #     print(i, self.Bulge[i].idx,self.Bulge[i].Weight)
        #
        # print("Assigning Group Weights (only Bulge assigned yet)")
        # self.WeighHPGroups()
        # for i in self.HPGroup._IDlist:
        #     print(i, self.HPGroup[i].HPID,self.HPGroup[i].Weight)
        #
        # print("Assigning Path Weights (only Bulge assigned yet)")
        # self.WeighPaths()
        # for i in self.Path._IDlist:
        #     print(i, self.Path[i].HPID,self.Path[i].Weight)
        self.WeighPaths()
        self.HP.PrintFeature()
        self.Bulge.PrintFeature()
        self.HPGroup.PrintFeature()
        self.Path.PrintFeature()


    def WeighPaths(self):
        print("Assigning weights in paths")
        PrevHPID = []

        for HPID in reversed(self.HP._IDlist): #Beginning at the end of the list (we will already have evaluated its children

            ### For every hairpin we need to evaluate its weight based on the following:
            #   1. the identity of its outermost bulge (PrevBulge)
            #       a. Bulge is a loop - 'isleaf'=='yes'
            #       b. Bulge is a bulge - 'isleaf'=='no' && 'siblingID'==[]
            #       c. Bulge is a junction - 'len(childrenID)' > 1
            #   2. The length of its next-inwards HP
            #   2. the length/type of hairpin (ROI - 1-2bp; gate - 3-4bp; stem - 5+ bp)
            #   3. The size of the HP group (large ROI/gate combos get exponentially harder)
            #   4. The type of junction

            #Looking at a hairpin, we want to know what came before it
            #PrevBulgeID = bf.getPropertyOfFeatsWithProperty(self.Bulge, 'parentHPID', 'ID', [HPID])
            #NextBulgeID = bf.getPropertyOfFeatsWithProperty(self.Bulge, 'childHPID', 'ID', [HPID])
            # print("b0",self.Bulge[0].childHPID,self.Bulge[0].parentHPID)
            # print("b1",self.Bulge[1].childHPID,self.Bulge[1].parentHPID)
            # print("b31",self.Bulge[31].childHPID,self.Bulge[31].parentHPID)
            # print("0",self.Bulge[35].childHPID,self.Bulge[35].parentHPID)
            # print("1", HPID,PrevBulgeID,NextBulgeID)
            PrevBulgeID = bf.getAttributeOfFeatFromMemberID(self.Bulge, 'parentHPID',[HPID],'ID','any')
            #print("2", HPID,PrevBulgeID,NextBulgeID)
            NextBulgeID = bf.getAttributeOfFeatFromMemberID(self.Bulge, 'childHPID',[HPID],'ID','any')
            #print("3", HPID,PrevBulgeID,NextBulgeID)
            PrevHPID = list(self.HP[HPID].childrenID)
            CurrLength = self.HP[HPID].length

            #If this hairpin is at the end, we are going to say that its previous length is 0 - i.e. the flexibility
            # allowable by its many children (or lack thereof) essentially makes it easy to stabilize.
            if self.HP[HPID].isleaf=='yes':
                PrevLength=0
            elif self.HP[HPID].isjxn=='yes':
                #print("This is a leaf or junction", HPID, "which was preceded by", PrevHPID)
                #print("Weights:", [self.HP[x].Weight for x in PrevHPID])
                PrevLength = sum([self.HP[x].length for x in PrevHPID])/(len(PrevHPID)+1)
            else: #if the hairpin has no siblings and has 1 child
               # print("This is NOT a leaf or junction", HPID, PrevHPID, self.HP[HPID].isleaf,self.HP[HPID].neighborID)
               try:
                    PrevLength = self.HP[PrevHPID[0]].length
               except:
                   PrevLength=0
            LengthMultiple = (((max(0, 5 - CurrLength)) + (max(0, 5 - PrevLength))))/2 #Makes consecutive short HPs difficult

            #This is the type of HP that we are in...this effects how weights stack
            CurrHPkind = self.HP[HPID].kind

            #We are collecting its neighboring HPGroup (i.e. the HPs that share regional (in)stability)
            HPGroupID = bf.getAttributeOfFeatFromMemberID(self.HPGroup,'HPID',[HPID],'ID','any')
            FGHPList = [self.HPGroup[x].HPID for x in HPGroupID][0]

            ClusterWeight = sum([self.Bulge[x].Weight for x in self.Bulge._IDlist if self.Bulge[x].parentHPID!=[] and self.Bulge[x].parentHPID[0] in FGHPList])


            #Analyzing the HPGroup. Essentially setting a baseline difficulty for the clusterweight, which will then be added to the current HP depending on its location within the group.
            # across all HPs (except weighted) in being applied across every member every time)
            if CurrHPkind == 'ROI':
                #ClusterWeight += 1.5 * (len(FGHPList) + 1) #Weight scales as 1.5x the number of short 1-2bp.
                ClusterWeight += 1.5 * (len(FGHPList) + 1)
            elif CurrHPkind == 'gate':
                #ClusterWeight += sum([5-self.HP[x].length for x in FGHPList]) #Weight scales as a weight +1 or +2 depending on length of ROI.
                ClusterWeight += (sum([5-self.HP[x].length for x in FGHPList]))
            elif CurrHPkind == 'stem':
                ClusterWeight = ClusterWeight **(1/2) + 0.2 * len(FGHPList)

            #This makes it so that split HPGroups (like long connected gates) have their outermost features weighted more heavily, scaling between x and x/len(FG)
            ProgenyOfHPInFG = [x for x in FGHPList if x in [HPID]+self.HP[HPID].progenyID]
            ClusterCorrection = (1+len(ProgenyOfHPInFG))**0.5

            #ClusterWeight = LengthMultiple*self.Bulge[PrevBulgeID[0]].Weight #A high bulge weight on short HPs is hard

            self.HP._weightmap[HPID] += ClusterWeight/ClusterCorrection

            # neighborIDs = list(self.HP[HPID].neighborID)
            # childrenID = list(self.HP[HPID].childrenID)
            # #Applies half of the cluster weight to all of its neighbors, and half to itself.
            # TempWeightMap = 0
            #
            # for i in childrenID:
            #     self.HP._weightmap[i] = self.HP._weightmap[i] + 0.5 * ClusterWeight / ((len(childrenID)))
            #     TempWeightMap += self.HP._weightmap[i]
            # TempWeightMap = (TempWeightMap/(1+len(childrenID)))**(1/(1+len(childrenID)))
            # TempWeightMap += self.HP._weightmap[HPID] + ClusterWeight/2

            childrenID = list(self.HP[HPID].childrenID)
            UpstreamHPWeight = 0
            for i in childrenID:
                UpstreamHPWeight += self.HP[i].Weight**(2) #1 of 3 to uncommend
                #UpstreamHPWeight += self.HP[i].Weight** (2.5 / self.HP[i].length)

            #print(HPID,self.HP[HPID].idx,self.Bulge[0].childHPID,NextBulgeID,PrevBulgeID)
            try:
                AllBulgeWeightNormal = (self.Bulge[PrevBulgeID[0]].Weight**2 + self.Bulge[NextBulgeID[0]].Weight**2)**0.5
            except:
                try:
                    AllBulgeWeightNormal = self.Bulge[PrevBulgeID[0]].Weight
                except:
                    AllBulgeWeightNormal = self.Bulge[NextBulgeID[0]].Weight

            #AllHPWeightNormal = (self.HP._weightmap[HPID] ** 2 + (UpstreamHPWeight) / (len(childrenID) + 1)) ** (1 / 2)
            AllHPWeightNormal = (self.HP._weightmap[HPID] ** (2.5 / CurrLength) + (UpstreamHPWeight) / (len(childrenID) + 1)) ** (1 / 2)
            #CurrHPWeight = (AllBulgeWeightNormal ** (2.5 / CurrLength) + AllHPWeightNormal ** (2.5 / CurrLength)) ** 0.5
            CurrHPWeight = (AllBulgeWeightNormal**2 + AllHPWeightNormal **2) ** 0.5
            self.HP._weightmap[HPID] += CurrHPWeight
        return self
        #
        #     CurrHPkind = self.HP[HPID].kind
        #     CurrHPLength = self.HP[HPID].length
        #     CurrHPWeight = self.HP._weightmap[HPID]
        #
        #     PrevBulgeWeight = self.Bulge[PrevBulgeID].Weight
        #
        #     NextHPLength = self.HP[self.HP[HPID].parentID].length
        #     LengthMultiple = ((max(0, 5 - CurrLength)) + (max(0, 5 - NextLength))) / 2
        #
        #     if self.HP[HPID].children==[]: #If we are looking at a leaf
        #         PrevBulgeWeight = self.Bulge[PrevBulgeID].BulgeWeight*LengthMultiplier
        #         PrevHPWeight = 1
        #         ###ELSE IF WE AREN'T LOOKING AT LEAF
        #         #print(PrevBulgeWeight, FG.Bulge[PrevBulgeID].BulgeWeight)
        #     elif len(self.HP[HPID].childrenID)==1:#len(FG.HP[HPID].neighborID)==2 and len(FG.HP[HPID].siblingID)==0:
        #         PrevHPID = self.HP[HPID].childrenID[0]
        #         PrevHPWeight = self.HP._weightmap[PrevHPID]
        #         PrevBulgeWeight = self.Bulge[PrevBulgeID].BulgeWeight*LengthMultiplier
        #     else:
        #         PrevHPID = self.HP[HPID].childrenID
        #         #print("-->",PrevHPID)
        #         AllHPWeights = 0
        #         for i in PrevHPID:
        #             PrevHPWeight = self.HP[i].Weight
        #             AllHPWeights = AllHPWeights + PrevHPWeight**2
        #         PrevHPWeight = AllHPWeights**(1/2)
        #         PrevBulgeWeight = self.Bulge[PrevBulgeID].BulgeWeight*LengthMultiplier
        #
        #     if CurrHPkind=='stem':
        #         #CurrHPWeight = (CurrHPWeight**2+PrevHPWeight**2)**(1/2)*PrevBulgeWeight
        #         CurrHPWeight = (PrevBulgeWeight*CurrHPWeight)**(1/2)
        #         #FG.Bulge[PrevBulgeID].Weight = (FG.Bulge[PrevBulgeID].BulgeWeight*PrevHPWeight)**1/2
        #         #self.Bulge._weightmap[PrevBulgeID] = (PrevBulgeWeight*PrevHPWeight)**(1/2)
        #     else:
        #         #FG.Bulge[PrevBulgeID].Weight = PrevHPWeight*PrevBulgeWeight
        #         #self.Bulge._weightmap[PrevBulgeID] = PrevHPWeight*PrevBulgeWeight
        #         #print(CurrHPWeight,PrevHPWeight,PrevBulgeWeight,LengthMultiplier)
        #         CurrHPWeight = ((CurrHPWeight**2+PrevHPWeight**2)+(PrevBulgeWeight)**2)**(1/2)
        #     self.HP[HPID].Weight = CurrHPWeight
        #     print("weighing HP:", HPID, CurrHPWeight, "vs",StartHPWeights[HPID], "vars:", PrevBulgeWeight, PrevHPWeight, LengthMultiplier)
        # self.Path._weightmap[PathID] = CurrHPWeight
        #         # print("PrevBulgeWeight", PrevBulgeWeight)
        #         # print("PrevHPWeight", PrevHPWeight)
        #         # print("CurrHPWeight", CurrHPWeight)
        #         # print("FG.HP[HPID+1].Weight", FG.HP[HPID+1].Weight)
        # return self

        #print(FG.HP._weightmap)


        # for HPID in FG.HP._IDlist:
        #     print(HPID, FG.HP[HPID].Weight)



        #PathsByHPIDs = bf.getPropertyOfFeatsWithProperty(FG.Path, 'isjxn', 'HPID', 'no')

        #for PathIDs in OutsidePathIDs
        #PathIDs = OutsidePathIDs[0]
        #HPIDs = FG.HP


        #LeafList = bf.getPropertyOfFeatsWithProperty(FG.HP, 'isleaf', 'ID', 'yes', FGList)

        #print(FGList,LeafList)


        #for


#These definitions are all "features" of a structure. a feature is just something that can be described by a list of indices given as a 'pairmap'

class Feature():

    def __init__(self,PairMap,name):

        ### Feature definition: Any combination of bases forming a coherent linkage or network. i.e.
        #   Entry   Name        Base layout     how to reference
        #   1.  Base (.)        - PairMap = [   0, 1, 2,                                    . = "1,"    ;   ]-[ = "],["
        #   2.  BasePair (BP)   - PairMap = [   [0,100],[1, 99],                            [-.-.-]
        #   3.  Hairpin (HP)    - PairMap = [   [[0, 100],[1, 99]],                         [-[..]-[..]-]
        #   4.  Bulge           - PairMap = [   [[0, 1, 2], [45], [98, 99]],                [-[[...][.][..]]-]
        #   5.  HPGroup         - PairMap = [   [[[0,100], [1, 99]], [[4, 96]]],            [-[Bulge-[BP-BP-BP]-Bulge-[BP-BP-BP]-Bulge-[BP-BP-BP]-Bulge]-]
        #   6.  Path            - PairMap = [   [[[[0,100], [1, 99]], [[4, 96]][[..]]]],    [-[HPG-HPG-HPG-]-]

        self.name = name #Name of the structure
        self.length = len(PairMap)
        self.parentseq = Sequence
        self._pairmap = PairMap #PairMap
        print(name," PairMap initialized ->",self._pairmap)
        self._IDlist = list(range(0, len(PairMap)))
        self._allbases = [sorted(bf.FlattenList(self._pairmap[x])) for x in list(range(0,len(self._pairmap)))]
        if self.name=='Base':
            return

        #Identifies the families based on the functions _getparentIDmap and _getchildrenIDmap which has a small amount of feature-specific code
        self._parentIDmap = [Feature._getparentIDmap(self, x) for x in self._IDlist]
        self._childIDmap = Sequence._getChildMapFromParentMap(self._parentIDmap)
        self._siblingIDmap = Sequence._getSiblingMapFromParentMap(self._parentIDmap)
        self._parentmap = [[self._pairmap[x] if x != [] else [] for x in y] for y in self._parentIDmap]
        self._childmap = [[self._pairmap[x] if x != [] else [] for x in y] for y in self._childIDmap]
        self._siblingmap = [[self._pairmap[x] if x != [] else [] for x in y] for y in self._siblingIDmap]
        self._neighborIDmap = self._getneighborIDmap(self._parentIDmap, self._childIDmap,self._siblingIDmap)
        self._neighbormap = [[self._pairmap[x] if x != [] else [] for x in y] for y in self._neighborIDmap]
        self._progenyIDmap = [self._getAllProgenyID(x) for x in self._IDlist]
        self._leafmap = [self._isLeaf(x) for x in self._IDlist]
        self._jxnmap =  [self._isJunction(x) for x in self._IDlist]
        self._openmap = [self._isOpening(x) for x in self._IDlist]
        self._weightmap = [0]*len(PairMap)
    def __getitem__(self, item,seqno=None):
        if item<0 or item>len(self._pairmap)-1:
            item = len(self._pairmap)-1
        if seqno is None:
            seqno = 1
        self.ID = item
        self.idx = self._pairmap[item]
        if self.name =='Base':
            return self
        self.allbases = self._allbases[item]
        self.parent = self._parentmap[item]
        self.children = self._childmap[item]
        self.sibling = self._siblingmap[item]
        self.parentID = self._parentIDmap[item]
        self.childrenID = self._childIDmap[item]
        self.siblingID = self._siblingIDmap[item]
        self.neighborID = self._neighborIDmap[item]
        self.progenyID = self._progenyIDmap[item]
        self.isleaf = self._leafmap[item]
        self.isjxn = self._jxnmap[item]
        self.isopening = self._openmap[item]
        self.Weight = self._weightmap[item]
        self.length = len(self._pairmap[item])
        return self
    def __len__(self):
        return len(self._IDlist)
    def __iter__(self):
        self.index = 0
        return self
    def __next__(self):
        if self.index==len(self):
            #print("Iteration Done")
            raise StopIteration
        self.index = self.index+1
        return self[self.index-1]
    def __str__(self):
        try:
            return str(self.ID)
        except:
            pass

    def _getAllProgenyID(self,Index):

        IncludedID = []
        RemainingIDs = [Index]
        while RemainingIDs != []:
            CurrID = RemainingIDs.pop(0)
            IncludedChildren = list(self._childIDmap[CurrID])
            IncludedID += [x for x in IncludedChildren if x not in IncludedID]
            RemainingIDs += IncludedChildren
        IncludedID = list(sorted(IncludedID))
        return IncludedID
    def _getparentIDmap(self, Index):
        #Inputs an ID in the pairmap. Must search through pairmap to identify parent.
        if self.name=='Base':
            return []
        elif self.name=='BP':
            ParentBP = Sequence.BPParentMap[Index]
            #print(self.name,"Parent of", Index,ParentBP)
            if ParentBP ==[] or ParentBP==[[]]:
                return []
            else:
                return ParentBP

        elif self.name=='HP':
            ParentHP = (Sequence.HPParentMap[Index])
            if ParentHP ==[] or ParentHP ==[[]]:
                return []
            else:
                return ParentHP

        elif self.name=='Bulge':

            ParentBulge = Sequence.BulgeParentMap[Index]
            if ParentBulge ==[] or ParentBulge==[[]]:
                return []
            else:
                return ParentBulge

            return ParentBulge

        elif self.name=='HPGroup' or self.name == 'Path':
            HPIDsInGroup = list(self._HPIDmap[Index])
            HPParentMap = list(Sequence.HPParentMap)
            ParentHP = HPParentMap[HPIDsInGroup[0]]
            return [self._HPIDmap.index(x) for x in self._HPIDmap if any(elem in ParentHP for elem in x)]
    def _getneighborIDmap(self,ParentIDMap,ChildrenIDMap,SiblingIDMap):
        NeighborIDMap = []
        #print(len(ParentIDMap),len(ChildrenIDMap),len(SiblingIDMap))
        for i in self._IDlist:
            Neighbors = sorted(ParentIDMap[i]+ChildrenIDMap[i]+SiblingIDMap[i])
            NeighborIDMap.append(Neighbors)
        return NeighborIDMap
    def _isLeaf(self, item):
        if self._childmap[item]==[]:
            return 'yes'
        else:
            return 'no'
    def _isJunction(self,item):
        if len(self._childmap[item])>1:#or len(self._siblingmap[item])>0:
            return 'yes'
        else:
            return 'no'
    def _isOpening(self,item):
        if self._parentmap[item]==[] or self._siblingIDmap[item]!=[]:
            return 'yes'
        else:
            return 'no'
    def PrintFeature(self):
        print(self.name, 'Summary')
        allvars = [key for key in self.__dict__]
        itemvars = [var for var in self[0].__dict__]
        print(itemvars)
        print('{:<4}'.format('ID'), '{:<15}'.format('length'), '{:<10}'.format('kind'), '{:<10}'.format('isleaf'),'{:<10}'.format('isjxn'),
              '{:<10}'.format('isopening'),
              '{:<12}'.format('parentID'),
              '{:<20}'.format('childrenID'),
              '{:<20}'.format('siblingID'),
              '{:<20}'.format('Weight'),
              '{:<50}'.format('Bases'))
        count = 0
        for i in self:
            print('{:<4}'.format(i.ID), '{:<15}'.format(str(i.length)), '{:<10}'.format(str(i.kind)),
                  '{:<10}'.format(str(i.isleaf)), '{:<10}'.format(str(i.isjxn)), '{:<10}'.format(str(i.isopening)),
                  '{:<12}'.format(str(i.parentID)),
                  '{:<20}'.format(str(i.childrenID)),
                  '{:<20}'.format(str(i.siblingID)),
                  '{:<20}'.format(str(i.Weight)),
                  '{:<50}'.format(str(i.idx)))

            count = count + 1
    #Dated - can delete _printme
    def _printme(self,item=None):

        if item==None:
            IDList = list(self._IDlist)
        else:
            IDList = [item]
        allvars = [key for key in self.__dict__]
        itemvars = [var for var in self[0].__dict__]
        itemvars = [x for x in itemvars if x not in allvars]
        try:
            itemvars.remove('parent')
            itemvars.remove('sibling')
            itemvars.remove('children')
        except:
            pass
        print(itemvars)
        for i in IDList:
            print([[x,getattr(self[i],x),"\n"] for x in itemvars])
class BaseStruct(Feature):

    def __init__(self,Sequence):
        name = 'Base'
        BaseMap = list(range(0, len(Sequence)))
        self.BaseMap = BaseMap
        Feature.__init__(self,BaseMap,name)
    def __getitem__(self, item,seqno=None):
        Feature.__getitem__(self,item)
        if item<0 or item>len(self._pairmap)-1:
            item = len(self._pairmap)-1
        if  item +1 > len(self._pairmap) - 1:
            self.next = 0
        else:
            self.next = item+1
        if item-1<0:
            self.prev = len(self._pairmap)-1
        else:
            self.prev = item-1
        return self
class BasePairStruct(Feature):

    def __init__(self,Sequence):
        #PartnerIdxs, self.BPMap = bf.getPartnerIdxs(Sequence[2])
        self.BPMap = Sequence.BPMap
        Feature.__init__(self,self.BPMap,'BP')
    def __getitem__(self, item,seqno=None):
        Feature.__getitem__(self,item)
        return self
class HairpinStruct(Feature):

    def __init__(self,Sequence):
        #HPMap = bf.getHPMap(Sequence[2])
        HPMap = Sequence.HPMap
        Feature.__init__(self,HPMap,'HP')
        self._leafmap = [self._isLeaf(x) for x in range(0,len(self._pairmap))]
        self._jxnmap =  [self._isJunction(x) for x in range(0,len(self._pairmap))]
        self._kindmap = [self._getHPkind(x) for x in range(0,len(self._pairmap))] #Defines how hairpins are later grouped into HPGroups
        self._bulge5pmap = list(Sequence.HPBulge5p)
        self._bulge3pmap = list(Sequence.HPBulge3p)
    def __getitem__(self, item,seqno=None):
        if item<0 or item>len(self._pairmap)-1:
            item = len(self._pairmap)-1
        Feature.__getitem__(self,item)
        self.length = len(self._pairmap[item])
        self.kind = self._kindmap[item]
        self.HPID = [self._pairmap.index(self.idx)]
        self.Bulge5p = self._bulge5pmap[item]
        self.Bulge3p = self._bulge3pmap[item]
        return self

    def _getHPkind(self,item):
        HPLength = len(self._pairmap[item])
        if HPLength<=2:
            return 'ROI'
        elif HPLength<=4:
            return 'gate'
        else:
            return 'stem'
class BulgeStruct(Feature):

    def __init__(self,Sequence):
        self._parentHPIDmap = Sequence.BulgeParentHPMap
        self._childHPIDmap = Sequence.BulgeChildHPMap
        Feature.__init__(self,list(Sequence.BulgeMap),"Bulge")
        #self._numopeningmap = [len(x) for x in self.BulgeInfo[0]]
        self._numopeningmap = [len(x) for x in self._pairmap]
        self._kindmap = [self._getBulgekind(x) for x in self._IDlist]
        self._lengthordermap = [[len(x) for x in self._pairmap[y]] for y in self._IDlist]
        self._numbasesmap = [len(bf.FlattenList(self._pairmap[x])) for x in self._IDlist]
        self._bulgeweightmap = [self._getBulgeWeight(x) for x in self._IDlist]
        self._weightmap = list(self._bulgeweightmap)
    def __getitem__(self, item,seqno=None):
        if item<0 or item>len(self._pairmap)-1:
            item = len(self._pairmap)-1
        Feature.__getitem__(self,item)
        self.NumOpenings = self._numopeningmap[item]
        self.LengthOrder = self._lengthordermap[item]
        self.NumBases = self._numbasesmap[item]
        self.kind = self._kindmap[item]
        #Bulge weight is the primary factor guiding pathfinding in deconstruction. See BulgeLengthvsWeight.png to visualize how they're weighted
        #This logic is the bulk of how difficult parts of the puzzle are determined. It may seem arbitrary but it isn't
        self.BulgeWeight = self._bulgeweightmap[item]
        self.parentHPID = self._parentHPIDmap[item]
        self.childHPID = self._childHPIDmap[item]
        self.Weight = self._weightmap[item]
        return self

    def getBulgesConnectingHPs(self, HP, HPList):
        # print("Hmmm", HPList)
        ChildBulges = bf.getAttributeOfFeatFromMemberID(Bulge, 'parentHPID', HPList, 'ID', 'any')
        ParentBulges = bf.getAttributeOfFeatFromMemberID(Bulge, 'childHPID', HPList, 'ID', 'any')
        return bf.CleanList([ChildBulges + ParentBulges for x in HPList])

    def _getBulgekind(self,item):
        NumOpenings = self._numopeningmap[item]
        if NumOpenings==1:
            return 'loop'
        elif NumOpenings==2:
            return 'bulge'
        else:
            return 'jxn'
    def _getBulgeWeight(self,item):
        AdjFactor = 3 #arbitrary scaling factor
        BulgeType = self._kindmap[item]
        LO = self._lengthordermap[item] #length order
        NumBases = self._numbasesmap[item]

        NumBulges = len(LO)
        if BulgeType=='loop':
            if NumBases==4:
                Weight = 0.2
            else:
                Weight = 0.15*NumBases*(1 + 2/NumBases + 4/NumBases**2)
            Weight = Weight
        elif BulgeType=='bulge':
            ##Rationale:
            #   Strain: Approaches DOA of asymmetric loops, but is zero when symmetric
            #   Relief: For an asymmetric loop, value approaches the degree of asymmetricity when min = 0 (DOA). Is <1 for most symmetric loops, or asymmetric ones that are longer
            #   Roll : Value is 1 for all symmetric loops, approach 0 linearly for high DOA, proportional to the difference
            #   Normal: Even more pronounced trend than Relief. Basically makes highly asymmetric loops very costly. As a divisor it will approximately normalize strain*relief

            Strain = max(LO) - min(LO) #if its asymmetrical                     0 (symm)    < - > Inf           (asymm)
            Relief = Strain / (min(LO) + 1) #Grants relief unless 0-n bulge     0 (min>>0)  < - > 1             (min==0)
            #Roll starts to weight long bulges as difficult the longer they get. This is especially important if the closing HPs are short
            Roll = (min(LO)+1)/(Strain+1) #essentially 1/(1+Relief).                0 (asymm, min==0) < - > Inf     (long, ~symm)
            Normal = 2 * Relief **2 + 1 #Becomes very large as symmetry
            Ramp = 1/((NumBases+1))**(1/2)
            BaseWeight = ((Relief*Strain*Roll)**2/(Normal)**(1/2))
            Weight = BaseWeight*Ramp + 0.1/Ramp #This landscape looks like a bird about to strike prey, with its body representing x=y and wings at x=0, y=0, which are maximum and peak at ~2-3

            Weight = Weight*AdjFactor
            #Weight maxes out at about ~0.4 for 2-0 bulge
        elif BulgeType=='jxn': #Bulge is a junction
            CombWeight = 0
            for Bulge1 in range(0,NumBulges-1):
                for Bulge2 in range(Bulge1+1,NumBulges):
                    LO1 = LO[Bulge1]
                    LO2 = LO[Bulge2]
                    if LO1==LO2==0: #This is the penalty for having zero-length bulges connecting two hairpins...they are very difficult!
                        BaseWeight=3
                    else:
                        LODuo = [LO1,LO2]
                        Strain = max(LODuo) - min(LODuo)  # if its asymmetrical                     0 (symm)    < - > Inf           (asymm)
                        Relief = Strain / (min(LODuo) + 1)  # Grants relief unless 0-n bulge     0 (min>>0)  < - > 1             (min==0)
                        # Roll starts to weight long bulges as difficult the longer they get. This is especially important if the closing HPs are short
                        Roll = (min(LODuo) + 1) / (Strain + 1)  # essentially 1/(1+Relief).                0 (asymm, min==0) < - > Inf     (long, ~symm)
                        Normal = 2 * Relief ** 2 + 1  # Becomes very large as symmetry
                        Ramp = 1 / (NumBases + 1)
                        BaseWeight = ((Relief * Strain * Roll) ** 2 / Normal) ** (1 / 2)
                    CombWeight +=BaseWeight
            Weight = max(0,(CombWeight)/NumBulges - 0.2*NumBulges)
        else:
            Weight = 0.5
        return Weight
class HPGroupStruct(Feature):

    def __init__(self,Sequence):
        [self._HPIDmap, self._kindmap] = list(bf.groupAllFeaturesBasedOnAttribute(HP, 'kind',None, 'near', None))
        self._bulgeIDmap = [Bulge.getBulgesConnectingHPs(HP,x) for x in self._HPIDmap]
        HPGroupPairMap = [[HPMap[x] for x in y] for y in self._HPIDmap]
        Feature.__init__(self,HPGroupPairMap,'HPGroup')
    def __getitem__(self, item,seqno=None):
        Feature.__getitem__(self,item)
        if item<0 or item>len(self._pairmap)-1:
            item = len(self._pairmap)-1
        self.NumHPs = len(self._pairmap[item])
        self.HPID = self._HPIDmap[item]
        self.BulgeID = self._bulgeIDmap[item]
        self.kind = self._kindmap[item]
        return self



class PathStruct(Feature):

    def __init__(self,Sequence):

        [self._HPIDmap, self._ismultijxnmap] = list(bf.groupAllFeaturesBasedOnAttribute(HP, 'isjxn', None, 'near', None))
        self._HPGroupMap = [bf.getAttributeOfFeatFromMemberID(HPGroup,'HPID',x,'ID','any') for x in self._HPIDmap]

        self._bulgeIDmap = [Bulge.getBulgesConnectingHPs(HP, x) for x in self._HPIDmap]
        PathMap = [[HPMap[z] for z in x]+[BulgeMap[q] for q in y] for [x,y] in zip(self._HPIDmap,self._bulgeIDmap)]
        Feature.__init__(self,PathMap,'Path')
        self._kindmap = ["-".join([HPGroup[x].kind for x in y]) for y in self._HPGroupMap]
    def __getitem__(self,item,seqno=None):
        Feature.__getitem__(self,item)
        if item<0 or item>len(self._pairmap)-1:
            item = len(self._pairmap)-1
        self.NumHPs = len(self._HPIDmap[item])
        self.NumBulges = len(self._bulgeIDmap[item])
        self.HPID = self._HPIDmap[item]
        self.BulgeID = self._bulgeIDmap[item]
        self.isMultiJxn = self._ismultijxnmap[item]
        self.kind = self._kindmap[item]
        return self

    def _getWeight(self,FG,item):
        Bulges = self._bulgeIDmap[item]
        HPs = self._HPIDmap[item]
        kind = self._ismultijxnmap[item] #This is 'yes' if its a 3WJ++
        print(HPs,Bulges)

        if kind=='yes':
            BulgeWeights = sum([FG.Bulge[y].Weight**2 for y in Bulges])
            HPWeights = sum([FG.HP[y].Weight**2 for y in Bulges])
            Weight = (BulgeWeights+HPWeights)**(1/2)
        else:
            BulgeWeights = sum([FG.Bulge[y].Weight for y in Bulges])
            HPWeights = sum([FG.HP[y].Weight for y in Bulges])
            Weight = (BulgeWeights**2 +HPWeights**2)**(1/2)
        if FG.HP[HPs[0]].kind=='stem':
            Weight = Weight**(1/2)



        return Weight







        #
        # HPsToReplaceByStableHP = [x for x in HPsToReplaceByStableHP if x not in StableStems]
        # print("Will add the following stable stems back in",StableStems)
        # print("and replace the following HPs", HPsToReplaceByStableHP)
        # HPLevel1 = bf.CleanList([x for x in HPLevel1 + StableStems if x not in HPsConnectingList])
        #
        # HPListByLevel.append(HPLevel1)
        # HPListByLevel.append(HPsToReplaceByStableHP)
        #
        # print("HPLBL222",HPListByLevel)
        #
        #
        #
        #
        #
        # print("testing GAOFFMID", getAttributeOfFeatFromMemberID(FG.Path, 'HPID', [3,4,5,6], 'ID', 'any', None))
        #
        # ## Identify HPs in the list that are
        # print("completed", getAttributeOfFeatFromMemberID(FG.HPGroup, 'HPID', HPLevel1, 'ID', 'all', None))
        # print("completed", bf.FlattenList([FG.HPGroup[x].HPID for x in [3,5]]))
        #
        # HPGroupsToFinish = self._IdentifyIncompleteFeaturesInHPList(HPsConnectingList, 'HPGroup')
        #
        #
        #
        #
        # PathsToFinish = self._IdentifyIncompleteFeaturesInHPList(HPsConnectingList, 'Path')
        #
        # print([FG.HPGroup[x].HPID for x in HPGroupsToFinish])
        # print([FG.Path[x].HPID for x in PathsToFinish])
        #
        # print("completing HPgroups:", HPGroupsToFinish)
        #
        # HPsConnectingList += self._CompleteFeaturesInHPList(HPsConnectingList, 'HPGroup')
        #
        #
        # #Find the hairpins that connect this list, minimally:
        # HPsConnectingList = bf.getAllHPsConnectingList(FG,HPList)
        #
        #
        # # #Identify HPs that are too unstable to function by themselves
        # # if CompleteUnfinishedJxnPathsFlag==1:
        # #     # Identify Incomplete Paths:
        # #     CompletedPaths = getAttributeOfFeatFromMemberID(FG.Path, 'HPID', HPList, 'ID', 'all', None)
        # #     AllPaths = getAttributeOfFeatFromMemberID(FG.Path, 'HPID', HPList, 'ID', 'any', None)
        # #     IncompletePaths = [x for x in AllPaths if x not in CompletedPaths and FG.Path[x].isjxn=='yes']
        # #
        # #     #Get nearby HPs that are not already in the list
        # #     OtherHPsToConsider = bf.FlattenList([FG.Path[x].HPID for x in IncompletePaths])
        # #     OtherHPsToConsider = [x for x in OtherHPsToConsider if x not in HPList]
        # #
        # #     #Ignore HPs in this list if they're long (i.e. stems)
        # #     HPsToRemove = [x for x in OtherHPsToConsider if FG.HP[x].kind=='stem'] #Identifies any stem stopping points in the path
        # #     HPsToRemove += [FG.HP._getAllProgenyID(x) for x in HPsToRemove if x not in HPList] #Identifies children of 'HPs that we need to remove
        # #
        # #     #   Removes from the HP List any HPs that belong in the above incomplete paths
        # #     OtherHPsToConsider = [x for x in OtherHPsToConsider if x not in HPsToRemove]
        # #     HPsConnectingList = [x for x in HPsConnectingList if x not in HPsToRemove]
        # #     HPsConnectingList = bf.CleanList(HPsConnectingList+OtherHPsToConsider)
        # #     HPsConnectingList = bf.getAllHPsConnectingList(FG, HPsConnectingList)
        # HPsConnectingList += self._CompleteFeaturesInHPList(HPsConnectingList, 'HPGroup')
        # HPsConnectingList = bf.getAllHPsConnectingList(FG,HPsConnectingList)
        #
        # LeafHPs = bf.getLeafHPsinHPList(FG,HPsConnectingList)
        # #LeafHPs = [x for x in LeafHPs if FG.HP[x].childrenID!=[] and x not in HPList]
        # LeafHPs = [x for x in LeafHPs if FG.HP[x].childrenID!=[] and x not in HPList]
        # #print("-->",HPList, LeafHPs)
        #
        # #Identify those that are likely too unstable to include as leafs:
        # HPsTooUnstable = [x for x in LeafHPs if FG.HP[x].length < MinStableHPLength and x not in HPList]
        # HPsToRemoveFromLevel1 = bf.CleanList([FG.HP[x].childrenID + [x] for x in HPsTooUnstable])
        # print(HPsConnectingList, "removing HPs", HPsToRemoveFromLevel1)
        #
        # #1. Identify all HPs that are within the same bulges etc. of HPList (i.e. minimum necessary HPs for puzzle similarity):
        # HPsConnectingList = [x for x in HPsConnectingList+HPList if x not in HPsToRemoveFromLevel1]
        # HPListByLevel.append(sorted([x for x in HPsConnectingList if x not in HPList]))
        #
        # #2. Identify HPs belonging to all sequence that are in the original list (i.e. those we will replace by stable HP)
        # HPsToReplaceByStableHP = bf.FlattenList([FG.HP[x].childrenID for x in HPList]+HPsTooUnstable)
        # HPListByLevel.append(sorted([x for x in HPsToReplaceByStableHP if x not in bf.FlattenList([HPListByLevel[1]+HPList])]))
        # HPsToReplaceByBulge = bf.FlattenList([FG.HP[x].childrenID for x in bf.FlattenList(HPListByLevel[1]+HPListByLevel[2])])
        # HPsToReplaceByBulge = [x for x in HPsToReplaceByBulge if x not in bf.FlattenList(HPListByLevel[0]+HPListByLevel[1]+HPListByLevel[2])]
        # HPListByLevel.append(sorted(HPsToReplaceByBulge))
        # return HPListByLevel


#main2()




 #
    # def _getSurroundingFeature(self,Direction_Flag): #MOSTLY WORKS, I JUST NEED TO MAKE IT RECOGNIZE ALL DOWNSTREAM FEATURES INSTEAD OF JUST ONE
    #
    #     Index=bf.FlattenList(self.idx)[0]
    #
    #     PairMap = self.PairMap
    #     PartnerIdxs = bf.getPartnerIdxs(Sequence[2])[0]
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
    #         #print("getIndexOfElement", bf.getIndexOfElement(BPMap, CurrIndex))
    #         while CurrIndex < len(Sequence):
    #             if bf.getIndexOfElement(PairMap, CurrIndex) == []:
    #                 if CurrIndex not in AllIndex and CurrIndex in bf.FlattenList(PairMap):
    #                     AllIndex.append(CurrIndex)
    #                 CurrIndex += 1
    #             else:
    #                 NewPair = bf.getPairFromBaseID(Sequence[2], CurrIndex)
    #                 CurrIndex = PartnerIdxs[CurrIndex] + 1
    #                 if NewPair in AllIndex:
    #                     try:
    #                         Indexes = [bf.getIndexOfElement(PairMap,x[0])[0] for x in AllIndex]
    #                     except:
    #                         Indexes = [bf.getIndexOfElement(PairMap,x)[0] for x in AllIndex]
    #                     print(Indexes)
    #                     return [PairMap[x] for x in Indexes]
    #                     #return bf.SortIndexList(AllIndex)
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
    #                 if CurrIndex not in AllIndex and CurrIndex in bf.FlattenList(PairMap):
    #                     AllIndex.append(CurrIndex)
    #                 CurrIndex -= 1
    #             else:
    #                 NewPair = bf.getPairFromBaseID(Sequence[2], CurrIndex)
    #                 CurrIndex = PartnerIdxs[CurrIndex] - 1
    #                 # print("found pair",NewPair,AllIndex,CurrIndex)
    #                 if NewPair in AllIndex:
    #                     try:
    #                         Indexes = [bf.getIndexOfElement(PairMap,x)[0] for x in AllIndex]
    #                         return [PairMap[x] for x in Indexes]
    #                     except:
    #                         Indexes = [bf.getIndexOfElement(PairMap,x) for x in AllIndex]
    #                         # Indexes = [bf.getIndexOfElement(BPMap,x)[0] for x in AllIndex]
    #
    #                     return bf.SortIndexList(AllIndex)
    #                 else:
    #                     AllIndex.append(NewPair)
    #             if CurrIndex <= 0:
    #                 CurrIndex = len(Sequence) - 1
    #     # AllBases = sorted(bf.CleanList(self.idx))
    #     # #AllBases = list(bf.find_ranges(AllBases))
    #     # AllOtherFeatures = list(self.BPMap)
    #     # AllOtherFeatures.remove(self.idx)
    #     # SurroundingFeatures = []
    #     #
    #     # if DirectionFlag==0:
    #     #     SurroundingPairs = sorted(bf.getAllSurroundingBPs(Sequence, AllBases[0], DirectionFlag))
    #     #     bp = SurroundingPairs.pop(0)
    #     #     bp = [AllBases[0], AllBases[-1]]
    #     #     foundflag = 0
    #     #     try:
    #     #         count = 1
    #     #         while foundflag == 0:
    #     #             idxF = bp[0]-count
    #     #             idxR = bp[1]+count
    #     #             count += 1
    #     #             for x in AllOtherFeatures:
    #     #                 if idxF in bf.FlattenList(x) or idxR in bf.FlattenList(x):
    #     #                     if foundflag == 0 and x not in SurroundingFeatures:
    #     #                         SurroundingFeatures.append(x)
    #     #                     foundflag = 1
    #     #     except:
    #     #         pass
    #     #
    #     #     # print("allbases",AllBases[0], self.idx)
    #     #     # for idx in reversed(range(0,AllBases[0])):
    #     #     #     if idx in bf.FlattenList(AllOtherFeatures):
    #     #     #         SurroundingFeatures = [x for x in AllOtherFeatures if idx in bf.FlattenList(x)]
    #     #     #         return SurroundingFeatures
    #     # else:
    #     #     SurroundingPairs = sorted(bf.getAllSurroundingBPs(Sequence,AllBases[0],DirectionFlag))
    #     #     if AllBases[0] in SurroundingPairs[0]:
    #     #         SurroundingPairs.pop(0)
    #     #     for bp in SurroundingPairs:
    #     #         foundflag = 0
    #     #         try:
    #     #             count = 0
    #     #             while foundflag ==0:
    #     #                 idxF = bp[0]+count
    #     #                 idxR = bp[1]-count
    #     #                 count+=1
    #     #                 for x in AllOtherFeatures:
    #     #                     if idxF+1 in bf.FlattenList(x) or idxR-1 in bf.FlattenList(x):
    #     #                         if foundflag==0 and x not in SurroundingFeatures:
    #     #                             SurroundingFeatures.append(x)
    #     #                         foundflag = 1
    #     #         except:
    #     #             pass
    #     #     print("Right here:-Dir=1->",SurroundingFeatures)
    #     #
    #     #     return SurroundingFeatures
    #     #
    #     #
    #     #         #
    #     #         #
    #     #         #
    #     #         #     for idx in range(bp[0],bp[1]):
    #     #         #         if idx in bf.FlattenList(AllOtherFeatures) and foundflag==0:
    #     #         #             foundflag=1
    #     #         # except:
    #     #         #     for idx in range(bp[0],bp[1]):
    #     #         #         if idx in bf.FlattenList(AllOtherFeatures) and foundflag==0:
    #     #         #             SurroundingFeatures.append([x for x in AllOtherFeatures if idx in bf.FlattenList(x) and x not in SurroundingFeatures])
    #     #         #             foundflag=1
    #     #
    #     #
    #     #
    #     #     #
    #     #     # print(AllBases,"SurroundingPairs,", SurroundingPairs)
    #     #     # BaseRanges = list(bf.find_ranges(AllBases))
    #     #     # print(BaseRanges)
    #     #     # SurroundingFeatures = []
    #     #     # for bprng in range(0,len(BaseRanges)-1):
    #     #     #     foundflag = 0
    #     #     #     print("trouble->",bprng,BaseRanges)
    #     #     #     try:
    #     #     #         for idx in range(BaseRanges[bprng][1],BaseRanges[bprng+1][0]):
    #     #     #             if idx in bf.FlattenList(AllOtherFeatures) and foundflag==0:
    #     #     #                 SurroundingFeatures.append([x for x in AllOtherFeatures if idx in bf.FlattenList(x)])
    #     #     #                 foundflag=1
    #     #     #     except:
    #     #     #         for idx in range(BaseRanges[bprng],BaseRanges[bprng+1]):
    #     #     #             if idx in bf.FlattenList(AllOtherFeatures) and foundflag==0:
    #     #     #                 SurroundingFeatures.append([x for x in AllOtherFeatures if idx in bf.FlattenList(x)])
    #     #     #                 foundflag=1
    #     #     # print("Right here:-->",SurroundingFeatures)
    #     #     # return SurroundingFeatures
    #     #
    #     #     #print(self,self.idx,SurroundingPairs,BaseRanges[0])
    #
    # def _getSurroundingFeatures(self,DirectionFlag,IndexToReturn=None):
    #     print("Looking at", self.idx, "And going in dir", DirectionFlag)
    #
    #     FeatureName = self.name
    #     print("Looking at", FeatureName, "#", self.ID, self.idx, "And going in dir", DirectionFlag)
    #     BasesInFeature = self.idx
    #     SurroundingFeatures = []
    #     FeatureDepth = bf.depth(self.idx)
    #     AllOtherFeatures = list(self.PairMap)
    #     AllOtherFeatures.remove(self.idx)
    #     if FeatureName=='Base':
    #         return []
    #     elif FeatureName=='BP':
    #         SurroundingBases = bf.getAllSurroundingBPs(Sequence,bf.CleanList(BasesInFeature)[0],DirectionFlag)
    #         return SurroundingBases
    #     elif FeatureName=='HP':
    #         if DirectionFlag==0:
    #             StartingPair = self.idx[0]
    #         else:
    #             StartingPair = self.idx[-1]
    #         SurroundingBases = bf.getAllSurroundingBPs(Sequence,StartingPair,DirectionFlag)
    #         if SurroundingBases ==None:
    #             return []
    #         # for BasePair in SurroundingBases:
    #         #     for Feature in AllOtherFeatures:
    #         #         if BasePair in Feature:
    #         #             SurroundingFeatures.append(Feature)
    #         # return(SurroundingFeatures)
    #         for Base in bf.FlattenList(SurroundingBases):
    #             SurroundingFeatures.append([x for x in AllOtherFeatures if Base in bf.CleanList(x) and x not in SurroundingFeatures])
    #         return(SurroundingFeatures)
    #        #else: return([x for x in AllOtherFeatures if
    #     elif FeatureDepth==1:
    #
    #         AllOtherFeatureRoutes = [bf.getAllSurroundingBPs(Sequence, x, DirectionFlag) for x in AllOtherFeatures]
    #         PossibleRoutes = sorted(bf.getAllSurroundingBPs(Sequence,BasesInFeature,DirectionFlag))
    #         print(DirectionFlag,"Entering:",AllOtherFeatureRoutes)
    #         print("PossibleRoutes before:", PossibleRoutes)
    #         if DirectionFlag==0:
    #             PossibleRoutes = [PossibleRoutes[0]]
    #         else:
    #             PossibleRoutes.pop(0)
    #         print("PossibleRoutes after:", PossibleRoutes)
    #
    #         for CurrRoute in PossibleRoutes:
    #             try:
    #                 SurroundingFeatures.append(AllOtherFeatures[AllOtherFeatureRoutes.index([x for x in AllOtherFeatureRoutes if CurrRoute in x][0])])
    #                 print("Found feature", SurroundingFeatures)
    #             except:
    #                 print("getAllNearbyBPAndBulge", bf.getAllNearbyBPAndBulge(Sequence,CurrRoute[0],DirectionFlag))
    #                 while CurrRoute is not None:
    #                     try:
    #                         SurroundingFeatures.append(AllOtherFeatures[AllOtherFeatureRoutes.index([x for x in AllOtherFeatureRoutes if CurrRoute in x][0])])
    #                         CurrRoute = None
    #                     except:
    #                         print("Current Route after not finding", CurrRoute)
    #                         CurrRoute = bf.getAllSurroundingBPs(Sequence,CurrRoute,DirectionFlag)
    #                         print("Current Route after not finding", CurrRoute)
    #
    #                         print("Current Route after not finding", [x for x in AllOtherFeatureRoutes if CurrRoute in x])
    #                         #SurroundingFeatures.append([])
    #         print("Returning:", SurroundingFeatures)
    #         return SurroundingFeatures
    #
    #
    #     #AllOtherBases = list(self.BPMap).remove(BasesInFeature)
    #
    #     #print(SurroundingBases)

    # def _getBulgeMap(self,HP):
    #
    #     # We are starting at the first HP, so the first bulge has children 0
    #     ParentIDList = [[]]
    #     NextHP = [0]
    #     ChildrenIDList = [[0]]
    #     BulgeList = [list(range(0, HP[0].idx[0][0])), list(range(HP[0].idx[0][1] + 1, len(Sequence)))]
    #     IDlist = list(HP._IDlist)
    #     BulgeListAll = [BulgeList]  # List containing all bulges separately
    #
    #     while NextHP != []:  # While we don't find any more to add:
    #         CurrHP = NextHP[0]
    #         ChildrenHP = HP[CurrHP].childrenID
    #         if ChildrenHP == []:
    #             BulgeList = list(range(HP[CurrHP].idx[-1][0] + 1, HP[CurrHP].idx[-1][1]))
    #         else:
    #             PairList = sorted(bf.FlattenList(bf.getAllSurroundingBPs(Sequence, HP[CurrHP].idx[-1], 1)))
    #             PairList = bf.getPairsFromList(PairList)
    #             BulgeList = []
    #             for i in PairList:
    #                 BulgeList.append(list(range(i[0] + 1, i[1])))
    #             NextHP.insert(0, ChildrenHP)
    #         NextHP.remove(CurrHP)
    #         NextHP = bf.CleanList(NextHP)
    #         BulgeListAll.append(BulgeList)
    #         ParentIDList.append([CurrHP])
    #         ChildrenIDList.append(ChildrenHP)
    #
    #     return BulgeListAll, ParentIDList, ChildrenIDList
    # #
    # # def _getSplitBulge(self, item):
    # #     BaseIDs = list(self._pairmap[item])
    # #     AllGroup = []
    # #     BaseGroup = []
    # #     base = BaseIDs[0]
    # #     while base in BaseIDs:
    # #         tempbase = base
    # #         while tempbase+1 in BaseIDs: #if we want to continue adding to baseIDs
    # #             BaseIDs.remove(tempbase+1)
    # #             BaseGroup.append(tempbase+1)
    # #             tempbase=tempbase+1
    # #         BaseIDs.remove(base)
    # #         BaseGroup.append(base)
    # #         BaseGroup = sorted(BaseGroup)
    # #         AllGroup.append(BaseGroup)
    # #         try:
    # #             base = BaseIDs[0]
    # #         except:
    # #             break
    # #         BaseGroup = []
    # #
    # #     return AllGroup
    # #
    # # def _findMissingBulge(self,item):
    # #
    # #     SplitBulge = self._getSplitBulge(item)
    # #     SurroundingStructures = bf.getSurroundingStructure(Sequence,SplitBulge[0],1)
    # #     MissingIndexList = []
    # #     for i in range(0,len(SurroundingStructures)-1):
    # #         if type(SurroundingStructures[i])==list and type(SurroundingStructures[i-1])==list:
    # #             MissingIndexList.append(i-1)
    # #     CountCorrect = 0
    # #     #print("Before",item, SplitBulge,MissingIndexList,SurroundingStructures,CountCorrect)
    # #     for newidx in MissingIndexList:
    # #         if newidx-CountCorrect<0:
    # #             SplitBulge.append([])
    # #         else:
    # #             SplitBulge.insert(newidx-CountCorrect,[])
    # #         CountCorrect = CountCorrect+1
    # #
    # #     #print("After",MissingIndexList, SurroundingStructures, CountCorrect, SplitBulge)
    # #     return SplitBulge
    # #     #
    # #     # for i in reversed(range(0,len(idx)-1)):
    # #     #     try:
    # #     #         print("1",idx,idx[i], [SurroundingStructures.index(x) for x in idx])
    # #     #     except:
    # #     #         print("2",i, idx, SurroundingStructures, [SurroundingStructures.index(x) for x in idx])
    # #     #     #SplitBulge.insert(SurroundingStructures.index(idx[i]),[])
    # #     #     #Surr
    # #     #
    # #     # print(idx, ":final:", SplitBulge)
    # #     # return SplitBulge
    # #     # print(SurroundingStructures, idx, MissingIndexList)
    # #     #for i in range(0,len(SplitBulge)-1):
#print(
# FirstBaseInFeature = self._pairmap[Index][0]
# SurroundingHPs = bf.getAllSurroundingHPs(Sequence, FirstBaseInFeature, 1)
# ParentHP = sorted(SurroundingHPs)[0]
# if ParentHP==[]:
#     return []
# FirstBaseInParentHP = HPMap[ParentHP][0]
# SurroundingStruct = bf.getSurroundingStructure(Sequence,FirstBaseInParentHP,0)
# ParentBases = [x for x in SurroundingStruct if type(x)!=list]
# #print(IndexList, "Surrounding Pairs", SurroundingPairs)
# if SurroundingHPs==[]:
#     return []
# else:
#     return []
# return []

    #print(FG.Bulge._findMissingBulge(12))
    #print(FG.Bulge._findMissingBulge(13))
    #print(FG.Bulge._findMissingBulge(15))
    #print(FG.Bulge[15].LengthOrder, FG.Bulge[15].NumBases,FG.Bulge[15].seq,FG.Bulge[15].BaseList,FG.Bulge[15].idx)

    #def _getPropertyFromIndexList(self,IndexList,property):
        #print("Receiving index list", IndexList, property)
        #return[getattr(self,property) for x in Feature._getIDfromindex(self, IndexList)]

    #def _getIDfromindex(self, IndexList):
        # if IndexList == []:
        #     return []
        #
        # try:
        #     return bf.CleanList([bf.getIndexOfElement(self._pairmap, x) for x in IndexList])
        #     # return [self.HP.PairMap.index(x) for x in IndexList]
        # except:
        #     try:
        #         return bf.getIndexOfElement(self._pairmap, IndexList)
        #     except:
        #         return []

# def getPartnerIdxs(SS):  # Input: '...((((...)).))..'; Output: [-1, -1, -1, 14, 13, 11, 10, -1, -1, -1, 6, 5, -1, 4, 3, -1, -1], [[3, 14], [4, 13], [5, 11], [6, 10]]
#     # Req'd Variables:
#     # SS
#     # Req'd Functions:
#     # None
#     List_Of_Pairs = []  # List of all base pair partnerslist
#     List_Of_Partner_Idx = [-1] * len(SS)  # List showing which base is paired with which
#     SS_Base_Count = 0  # Which base we are looking at
#     SS_Open_Idx = []  # Index of last open pair
#     SS_Closed_Idx = []  # Index of last closing pair
#     for i in list(SS):  # Identifies List_Of_Partner_Idx of each base
#         if i == '(':
#             SS_Open_Idx.append(SS_Base_Count)
#             Open_Pair_Flag = 0
#         elif i == ')':
#             SS_Closed_Idx.append(SS_Base_Count)
#             Open_Pair_Flag = 1
#         elif i == '.':
#             Open_Pair_Flag = 0
#
#         if Open_Pair_Flag == 1:
#             List_Of_Partner_Idx[SS_Open_Idx[-1]] = SS_Base_Count
#             List_Of_Partner_Idx[SS_Base_Count] = SS_Open_Idx[-1]
#             List_Of_Pairs.append([SS_Open_Idx[-1], SS_Closed_Idx[-1]])
#             SS_Open_Idx.remove(SS_Open_Idx[-1])
#             SS_Closed_Idx.remove(SS_Closed_Idx[-1])
#         SS_Base_Count = SS_Base_Count + 1
#     List_Of_Pairs = sorted(List_Of_Pairs)
#     return List_Of_Partner_Idx, List_Of_Pairs
#
#
# def main():
#
#     global Sequence
#     Sequence = SequenceStruct()
#
#     # Defines the first feature - the base (length of sequence):
#     BaseMap = list(range(0,len(Sequence)))
#     Base = BaseStruct(BaseMap, 'Base')
#     print("Number of Bases:", len(Base))
#     print("Number of Bases:", len(Base))
#
#     # Gets the BasePair (BP) BPMap -> [1, 352] signifying 1-bp-352
#     PartnerIdxs, BPMap = bf.getPartnerIdxs(Sequence[2])
#     BP = BasePairStruct(BPMap, 'BP')
#     print("Number of BPs:", len(BPMap))
#     print("Number of BPs:", len(BP))
#
#     # Gets the Hairpin (HP) BPMap -> [[1, 14],[2,13]] signifying a 2bp Hairpin
#     HPMap = bf.getHPIdxList(Sequence[2])
#     HP = HairpinStruct(HPMap, 'HP')
#     print("Number of hairpins:", len(HPMap))
#     print("Number of hairpins:", len(HP))
#
#     BulgedBases = bf.getBulgeListOfBases(Sequence)[0]
#     Bulge = Feature(BulgedBases, 'Bulge')
#     print("Number of bulges:", len(BulgedBases))
#     print("Number of bulges:", len(Bulge))
#
#     Stem = Feature([x for x in HPMap if len(x)>=5], 'Stem')
#     Gate = Feature([x for x in HPMap if len(x) in [3,4]], 'Gate')
#     ROI = Feature([x for x in HPMap if len(x)<3], 'ROI')
#
#     #for x in Bulge:
#         #print(len(bf.getAllSurroundingHPs(Sequence, x.idx, 1))>2)
#
#     ForkList = [x for x in BulgedBases if len(bf.getAllSurroundingHPs(Sequence, x, 1))>2]
#     print(ForkList[0])
#     print("FORK LIST:",ForkList)
#     Fork = Feature(ForkList, 'Fork')
#     #Fork = Feature([x for x in BulgedBases if ])
#
#     NumBases = len(Sequence)
#     NumBPs = len(BP)
#     NumHPs = len(HP)
#     NumBulges = len(Bulge)
#     BaseID = random.randint(0,NumBases-1)
#     BPID = random.randint(0,NumBPs-1)
#     HPID = random.randint(0,NumHPs-1)
#     BulgeID = random.randint(0,NumBulges-1)
#     #BulgeID =33
#     #HPID = 39
#
#     print("ID:", Base[BaseID].ID,BP[BPID].ID,HP[HPID].ID,Bulge[BulgeID].ID)
#     print("idx:", Base[BaseID].idx,BP[BPID].idx,HP[HPID].idx,Bulge[BulgeID].idx)
#     print("name:", Base[BaseID].name,BP[BPID].name,HP[HPID].name,Bulge[BulgeID].name)
#     print("seq:", Base[BaseID].seq,BP[BPID].seq,HP[HPID].seq,Bulge[BulgeID].seq)
#
#     print("Testing the getSurroundingFeatures Function:\n")
#     print("ID:", Base[BaseID].ID, BP[BPID].ID,HP[HPID].ID,Bulge[BulgeID].ID)
#     print("Family:", "parent:", HPID, HP[HPID].Parent, "siblings:", HP[HPID].Siblings,"children:",  HP[HPID].Children)
#     print("     :", "parent:","siblings:","children:")
#     print("HP", HPID, HP.idx)
#     print("Family:", [HP[x].seq for x in Feature._getIDfromindex(HP,HP[HPID].Parent)], [HP[x].seq for x in Feature._getIDfromindex(HP,HP[HPID].Siblings)], [HP[x].seq for x in Feature._getIDfromindex(HP,HP[HPID].Children)])
#     print("      :",HP._getPropertyFromIndexList(HP[HPID].Parent, "seq"),HP._getPropertyFromIndexList(HP[HPID].Siblings, "seq"),HP._getPropertyFromIndexList(HP[HPID].Children, "seq"))
#     print("      :",HP.Parent, HP.Siblings,HP.Children)
#     print("      :",HP.ParentID, HP.SiblingsID,HP.ChildrenID) #Notice these give the same, because when HPID is defined, its global?
#     print("      :",HP[HPID].ParentID, HP[HPID].SiblingsID,HP[HPID].ChildrenID) #Notice these give the same, because when HPID is defined, its global?
#     print("      :",HP[HP[HPID].ParentID[0]].seq) #Gets parent seq
#
#
#     print("BP", BPID, BP[BPID].idx)
#     print("Family:", [BP[x].seq for x in Feature._getIDfromindex(BP,BP[BPID].Parent)], [BP[x].seq for x in Feature._getIDfromindex(BP,BP[BPID].Siblings)], [BP[x].seq for x in Feature._getIDfromindex(BP,BP[HPID].Children)])
#     print("      :",BP._getPropertyFromIndexList(BP[BPID].Parent, "seq"),BP._getPropertyFromIndexList(BP[BPID].Siblings, "seq"),BP._getPropertyFromIndexList(BP[BPID].Children, "seq"))
#     print("      :",BP.Parent, BP.Siblings,BP.Children)
#     print("      :",BP.ParentID, BP.SiblingsID,BP.ChildrenID) #Notice these give the same, because when HPID is defined, its global?
#     print("      :",BP[BPID].ParentID, BP[BPID].SiblingsID,BP[BPID].ChildrenID) #Notice these give the same, because when HPID is defined, its global?
#     print("      :",BP[BP.ParentID[0]].seq)
#
#     print(BP.idx, BPID, BP.ID, BP.Parent,BP.seq,BP.ParentID,BP.Children,BP.ChildrenID)
#
#     BPID = 10
#     BP[BPID]
#     print(BP.idx, BPID, BP.ID, BP.Parent,BP.seq,BP.ParentID,BP.Children,BP.ChildrenID)
#
#
#     HP[0]
#     while HP.isLeaf=='no':
#         HPID = HP.ChildrenID[0]
#         HP[HPID]
#
#     print("ID       length      kind        ParentID        SibID       ChildID")
#     for i in range(0,len(HPMap)-1):
#         HP[i]
#         print(HP.ID, HP.length, HP.seq, HP.kind,HP.ParentID, HP.SiblingsID,HP.ChildrenID, HP.idx)
#
#
#     # #   THIS IS A BASIC PATHFINDER. BEGIN WITH HP I AND, IF YOU MEET THE CONDITION, PROCEED.
#     # i = 1
#     # while i in range(0,len(HPMap)-1):
#     #     print("Looking here:")
#     #     print(HP[HPID], HPID)
#     #     Testing = HP
#     #    #Testing = FeatureMinion(HP)
#     #     #minionHP = HairpinStruct(HPMap)
#     #     minionHP = copy.deepcopy(HP)
#     #     HPID = HPID+1
#     #     Testing[HPID]
#     #     print(HP, HPID, HP.HPidx, Testing[HPID].HPidx)
#     #     print("Entering into HP:",i, HPID,"NumChildren: ", len(HP.ChildrenID))
#     #     print(HP.Children,HP.ChildrenID)
#     #
#     #
#     #     if len(HP.ChildrenID)>1:
#     #         ChildLengths = []
#     #         for ChildHP in HP.ChildrenID:
#     #             ChildLengths.append(HP[ChildHP].length)
#     #             HP[i] #Need to somehow get child lengths without having to call HP...
#     #         Shortest = ChildLengths.index(min(ChildLengths))
#     #         nextidx = HP.ChildrenID[Shortest]
#     #         print("Children IDs:", HP.ChildrenID,ChildLengths,"shortest:",Shortest, "->",nextidx)
#     #         i = nextidx
#     #     elif len(HP.ChildrenID)==0:
#     #         break
#     #     else:
#     #         i = i+1
#
#     print("Testing Pathfindin", HPID,HP._getPathFromIndex(2,-1))
#
#     i = 0
#     CurrFeature = []
#     FeaturePairMap = []
#     CurrFeature = [i]
#     Minion = copy.deepcopy(HP)
#     Minion[i]
#     while i in range(0,Minion.NumHPs-1):
#         idxstart = i
#         idxend = HP._getPathFromIndex(i,1)
#         i = idxend
#         print("s:e",idxstart,idxend)
#         if idxstart==idxend:
#             FeaturePairMap.append(CurrFeature)
#             CurrFeature = []
#         else:
#             for feat in range(idxstart,idxend):
#                 CurrFeature.append(feat)
#         i = i+1
#
#     print(FeaturePairMap)
#
#     # i = 0
#     # CurrFeature = []
#     # FeaturePairMap = []
#     # CurrFeature = [i]
#     # Minion = copy.deepcopy(HP)
#     # FeatureKind = Minion[0].kind
#     # Minion[i]
#     # FGPairMap = []
#     # count = 0
#     #
#     # while i in range(0, Minion.NumHPs - 1):
#     #     count=count+1
#     #     ChildID = Minion[i].ChildrenID
#     #     #print("CHILDI,",ChildID,[x for x in ChildID], Minion[ChildID[0]].kind==FeatureKind)
#     #     NextHPs = [x for x in ChildID if Minion[x].kind==FeatureKind]
#     #     SibIDs = Minion[i].SiblingsID
#     #     NextHPs2 = [x for x in SibIDs if Minion[x].kind==FeatureKind]
#     #     NewFeatures = [x for x in bf.FlattenList(NextHPs+NextHPs2) if x not in CurrFeature and x not in bf.FlattenList(FGPairMap)]
#     #     CurrFeature.append(NewFeatures)
#     #     while len(NewFeatures)>0 and Minion._getPathFromIndex(i,1)!=[i]:
#     #         count=count+1
#     #         #print("ree")
#     #         ChildID = bf.FlattenList([Minion[x].ChildrenID for x in NewFeatures])
#     #         NextHPs = [x for x in ChildID if Minion[x].kind==FeatureKind]
#     #         NewFeatures =  [x for x in NextHPs if x not in CurrFeature and x not in bf.FlattenList(FGPairMap)]
#     #         CurrFeature.append(NewFeatures)
#     #     FGPairMap.append(sorted(bf.FlattenList(CurrFeature)))
#     #     flag = 1
#     #     i = i + 1
#     #     CurrFeature = [i]
#     #     FeatureKind = Minion[i].kind
#     #     while flag:
#     #         count=count+1
#     #         if i not in bf.FlattenList(FGPairMap) or i>=Minion[i].NumHPs - 1:
#     #             i = i - 1
#     #             Minion[i]
#     #             flag = 0
#     #         i = i + 1
#     #     CurrFeature = [i]
#     #     FeatureKind = Minion[i].kind
#     #
#     # print(count)
#     # print("FGPairMap",FGPairMap)
#     # print(Minion[29].kind, Minion[30].kind,Minion[1].kind)
#
#     FG = HPGroupStruct(Feature, 'FG')
#     print(FG.PairMap)
#
#
#
#
#
#
#     print("first hp:", HP.ID, HP.isLeaf,HP.seq,HP.idx)
#     #HP = HairpinStruct(HPMap,'HP')
#
#
#
#
#     HPID = 1
#     BulgeID = 3
#
#     # HP[HPID]._getSurroundingFeature(0)
#     # Bulge[BulgeID]._getSurroundingFeature(0)
#     # ROI[0]._getSurroundingFeature(1)
#     # ROI[-1]._getSurroundingFeature(0)
#
#     # print(bf.getSurroundingFeatureInStructure(Sequence,bf.FlattenList((ROI[0].idx[0])),0,ROI.PairMap))
#     # SurroundingFeat = bf.getSurroundingFeatureInStructure(Sequence,bf.FlattenList((ROI[0].idx[0])),0,ROI.PairMap)
#     # print(bf.getIndexOfElement(SurroundingFeat,ROI.PairMap))
#     #
#     # print(Fork[2].idx, "Looking at this fork <-")
#     # Fork[2]._getSurroundingFeature(1) #TROUBLE IS THAT I DONT GET ALL THE SURROUNDING ,ONLY ONE
#     #
#     # print(bf.getSurroundingFeatureInStructure(Sequence,363,1,ROI.PairMap))
#     #
#     #
#     # getFGPairMap(Sequence)







