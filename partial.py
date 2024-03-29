import RNAfeatures as FD
import basicfunc as bf
import input
import argparse

class HairpinList():

    def __init__(self, FG, HPList, minsubHPlength, subHPlength, maxlength):

        assert(minsubHPlength is not None)
        #Initialization will give it a list of hairpins and the FG that it belongs to.
        #This will then be "finished" by adding in missing HPs and tracking available routes, current sequence, length, etc.
        try:
            HPList[0]
        except:
            HPList = [HPList]
        
        self.maxlength = maxlength
        self.minsubHPlength = minsubHPlength
        self.subHPlength = subHPlength
        assert(self.minsubHPlength is not None)
        self.FG = FG

        #print(FG.HP[8].progenyID)
        if HPList==[[]]:
            self.HPList = []
            self.StartHPList = []
            self.PrevHPList = []
            self.AddedLast = 'na'
            return
        self.StartHPList = HPList
        self.PrevHPList = list([self.StartHPList])
        self.HPList = self._getAllConnectingHPs(HPList,0)
        #
        # print("all",self.StartHPList)
        self.categorizeHPs()
        if len(self) <= self.maxlength:
            self.AddedLast='yes'
        else:
            self.AddedLast='no'
            self._resetHPList()

#Used to categorize all of the hairpins in self.HPList
    def _addHPList(self,HPListToAdd):
        """
        Takes either an int (converted to a list in place) or a list of ints
        Filters that list to ensure it doesn't overlap with HPList
        """

        if type(HPListToAdd)==int:
            HPListToAdd = [HPListToAdd]
        HPListToAdd = bf.CleanList([x for x in HPListToAdd if x not in self.HPList])
        if HPListToAdd==[]:
            return
        
        # PrevHPList is a list of all previous HPLists held by this HairpinList
        self.PrevHPList.append(list(self.HPList))
        
        # Clean the sum of InnerHPs and the new, non-intersecting-with-HPList stuff to be added
        self.HPList = bf.CleanList(self.InnerHPs+HPListToAdd)

        self.HPList = self._getAllConnectingHPs(self.HPList,0)
        self.categorizeHPs()

        if len(self)>self.maxlength: #Checks to see if the length is too long
            self._resetHPList()
            self.AddedLast = 'no'
        else:
            self._checkHPList() #Checks things like whether we added new HPs or not
        # if self.AddedLast=='yes':
        #     print(self._getSeq(1))
        #     print(self._getSeq(2))
        # AMW: don't return self; we never use it and that's confusing. This mutates the object.
        #return self

    def _resetHPList(self):
        self.HPList = list(self.PrevHPList[-1])
        self.categorizeHPs()
        try:
            self.PrevHPList.pop(-1)
        except:
            self.PrevHPList = list(self.StartHPList)

        return
    def _checkHPList(self):
        #print(all(elem in self.PrevHPList[-1] for elem in self.HPList), self.HPList,self.PrevHPList[-1])
        if all(elem in self.PrevHPList[-1] for elem in self.HPList): #Checking to see if we added any new HPs
            self._resetHPList()
            self.AddedLast = 'no'
        else:
            self.lastAddedHP = [x for x in self.HPList if x not in self.PrevHPList[-1]]
            self.AddedLast = 'yes'

#Used to identify common features of HP groups
    def update_parent_HP(self):
        try:
            self.ParentHP = self.findParentHP()
        except:
            self.ParentHP = self.findParentHP(start=True)
    
    def findParentHP(self, start=False):
        temp_HPList = None
        if start: temp_HPList = self.StartHPList
        else: temp_HPList = self.HPList
        if type(temp_HPList) == int:
            temp_HPList = [int]
        temp_HPList = bf.CleanList(temp_HPList)
        for i in temp_HPList:
            progenyIDs = list(self.FG.HP[i].progenyID + [i])
            if all(elem in progenyIDs for elem in temp_HPList):
                return i
        ParentIDs = [self.FG.HP[x].parentID for x in temp_HPList]
        ParentIDs = [x for x in ParentIDs if x != []]
        minParent = min(ParentIDs)[0]
        for i in list(reversed(range(0, minParent + 1))):
            progenyIDs = list(self.FG.HP[i].progenyID + [i])
            if all(elem in progenyIDs for elem in HPList):
                return [i]
        return []

    def _getAllConnectingHPs(self, HPList,IncludeParentFlag=0):
        FG = self.FG
        ParentHP = self.findParentHP(HPList)
        if type(ParentHP)==int:
            ParentHP = [ParentHP]

        if ParentHP not in HPList and IncludeParentFlag ==1 and ParentHP!=[]:
            HPList.insert(0, ParentHP)
        # Gets all connecting HPs, searches for their children, and then repeats to be sure

        HPsConnectingList = bf.FlattenList([HPList])
        # for j in range(0,1): #Repeats twice
        for i in range(0, len(HPList) - 1):
            for j in range(i, len(HPList)):
                HPsConnectingList.append(self.getPathBetweenHPs([HPList[i], HPList[j]], 0))
                HPsConnectingList = bf.CleanList(HPsConnectingList)

        # for x in HPsConnectingList:
        #     print(":", any(elem in FG.HP[x].childrenID for elem in HPsConnectingList))
        # print("LeafHPs", HPList, "are as follows:", LeafHPs)
        NeighboringHPs = [FG.HP[x].neighborID for x in HPsConnectingList]
        HPsConnectingList = bf.CleanList(HPsConnectingList + NeighboringHPs)
        return list(sorted(HPsConnectingList))
    # {  Gives a list of HPs to traverse between any two hairpins
    def getPathBetweenHPs(self, PairOfHPs, IncludeParentFlag):
        FG = self.FG

        # print("TSIng2", ListOfHPs, getPathToHome(ListOfHPs[0]), getPathToHome(2))
        PathToHomeOne = self.getPathToHome(PairOfHPs[0])
        PathToHomeTwo = self.getPathToHome(PairOfHPs[1])
        UnionOfPaths = set(PathToHomeOne) & set(PathToHomeTwo)
        # Intersect - Union gives only those that are unique

        # print(PairOfHPs,"Paths:", PathToHomeOne,PathToHomeTwo)
        HPsinPath = list((set(PathToHomeOne) | set(PathToHomeTwo)) - (set(PathToHomeOne) & set(PathToHomeTwo)))

        # print("TSIng", ListOfHPs, HPsinPath, PathToHomeOne, PathToHomeTwo)
        if IncludeParentFlag == 1 and UnionOfPaths != set():
            HPsinPath.append(max(UnionOfPaths))
        if HPsinPath == []:
            return (PairOfHPs)
        return HPsinPath
    # {  Gives a List of HPs to traverse in order to get home
    def getPathToHome(self, HPID):
        ParentID = self.FG.HP[HPID].parentID
        PathToHome = list()
        PathToHome.append(HPID)
        while ParentID != []:
            HPID = ParentID[0]
            PathToHome.append(HPID)
            ParentID = self.FG.HP[HPID].parentID
        # if len(ParentHP)==1:
        #     return HPNo
        return PathToHome

#Identifies the special types of HPs needed to define a puzzle
    def categorizeHPs(self):
        if self.HPList==[]:
            self.ParentHP = []
            self.InnerHPs = []
            self.OuterHPs = []
            self.StapleHPs = []
            self.SubHPs = []
            self.ShellHPs = []
        else:
            self.update_parent_HP()
            self.InnerHPs = self.getInnerHPs()
            self.OuterHPs = self.getOuterHPs()
            self.StapleHPs = self.getStapleHPs()
            self.SubHPs = self.getSubHPs()
            self.ShellHPs = self.getShellHPs()
        return self
    def getInnerHPs(self):
        #Want to return all HPs that have their children, siblings, parents in the list
        ConnectedHPs = self.HPList
        #for HPID in ConnectedHPs:
        #print(ConnectedHPs)
        InnerHPs = [x for x in ConnectedHPs if all(elem in ConnectedHPs for elem in self.FG.HP[x].neighborID) or self.FG.HP[x].isleaf=='yes']
        #print("StarT",self.StartHPList,"Connected",ConnectedHPs,"inside",InnerHPs)
        return InnerHPs
    def getOuterHPs(self):
        #Returns all HPs that are dependent on substitutions (i.e. are substituted, or neighbors of subbed HPs)
        #InnerHPs = self.InnerHPs
        #print("StarT",self.StartHPList,"Connected",self.HPList,"inside",InnerHPs, "outer",[x for x in self.HPList if x not in self.InnerHPs])
        return [x for x in self.HPList if x not in self.InnerHPs]
    def getStapleHPs(self):
        #Returns all HPs that are included in the minimal list but that are being replaced due to being too short (or too long and we don't want all their children)
        return [x for x in self.OuterHPs if self.FG.HP[x].length <= self.minsubHPlength or self.FG.HP[x].length > 5]
    def getSubHPs(self):
        #Returns all HPs that are included in the minimal list but that are being replaced due to being too short
        return [x for x in self.OuterHPs if x not in self.StapleHPs]
    def getShellHPs(self):
        #Returns HPs that are being replaced by bulges
        BulgeHPs = bf.CleanList([self.FG.HP[x].neighborID for x in self.SubHPs])
        BulgeHPs = [x for x in BulgeHPs if x not in self.InnerHPs+self.OuterHPs]
        return BulgeHPs


#Functions used to assemble the sequence from the HPList
    def _getBulgeIdx(self):
        FG = self.FG
        InnerHPs = list(self.InnerHPs) #HPs that are included and are sufficiently isolated
        OuterHPs = list(self.OuterHPs)
        SubHPs = list(self.SubHPs) #Outer HPs that are sufficiently long to directly replace safely
        ShellHPs = list(self.StapleHPs) #Outer HPs that are not long enough to directly replace
        BulgeHPs = list(self.ShellHPs) #OuterHPs that
        Bulges = [FG.HP[x].Bulge5p+FG.HP[x].Bulge3p for x in InnerHPs+SubHPs]
        Bulges += [FG.HP[x].Bulge5p for x in BulgeHPs+ShellHPs]
        return list(sorted(bf.CleanList(Bulges)))
        #
        #
        #
        # if ParentHP==[]:
        #     DownwardBulges += [bf.getPropertyOfFeatsWithProperty(FG.Bulge, 'parentHPID', 'ID', [x]) for x in InnerHPs+SubHPs]
        #     UpwardBulges += [bf.getPropertyOfFeatsWithProperty(FG.Bulge, 'childHPID', 'ID', [x]) for x in StapleHPs]
        # else:
        #     if ParentHP in StapleHPs:
        #         StapleHPs.remove(ParentHP)
        #         UpwardBulges += Sequence.HPBulge5p[ParentHP[0]]
        #     DownwardBulges += [bf.getPropertyOfFeatsWithProperty(FG.Bulge, 'parentHPID', 'ID', [x]) for x in
        #                       InnerHPs + SubHPs]  # going from start outwards
        #     UpwardBulges += [bf.getPropertyOfFeatsWithProperty(FG.Bulge, 'childHPID', 'ID', [x]) for x in StapleHPs]  # going towards start
        #
        # return DownwardBulges+UpwardBulges
        # ChildBulges =  [bf.getPropertyOfFeatsWithProperty(FG.Bulge, 'parentHPID', 'ID', [x]) for x in HPList]
        # ParentBulges =  [bf.getPropertyOfFeatsWithProperty(FG.Bulge, 'childHPID', 'ID', [x]) for x in HPList]
        # return bf.FlattenList(ChildBulges+ParentBulges)
    def _getHPIdx(self):
        return list(sorted(bf.CleanList([self.FG.HP[x].idx for x in self.InnerHPs+self.SubHPs])))
    def _getSubstituteHP(self,seqno,HPID,insideorout = 'out'):
        SubsHP = input.SubsHP[seqno]
        HPLen = self.FG.HP[HPID].length

        if insideorout=='in':
            idxlist = self.FG.HP[HPID].idx
            idxlist5p = [x[0] for x in self.FG.HP[HPID].idx]
            idxlist3p = [x[1] for x in self.FG.HP[HPID].idx]
            Sub5p = None
            if seqno == "master_numbering":
                Sub5p = [it for it in SubsHP[0:max(0,self.subHPlength-HPLen)]]
                Sub5p.extend(bf.FlattenList([self.FG.Sequence.SeqList[seqno][x] for x in sorted(bf.FlattenList(idxlist5p))]))
            else:
                Sub5p = SubsHP[0:max(0,self.subHPlength-HPLen)] + "".join([self.FG.Sequence.SeqList[seqno][x] for x in sorted(bf.FlattenList(idxlist5p))])
            Sub3p = None
            if seqno == "master_numbering":
                Sub3p = [it for it in [self.FG.Sequence.SeqList[seqno][x] for x in sorted(bf.FlattenList(idxlist3p))]]
                Sub3p.extend(bf.FlattenList(SubsHP[min(0,-self.subHPlength+HPLen):len(SubsHP)]))
            else:
                Sub3p = "".join([self.FG.Sequence.SeqList[seqno][x] for x in sorted(bf.FlattenList(idxlist3p))]) + SubsHP[min(0,-self.subHPlength+HPLen):len(SubsHP)]
            #Sub3p = "".join([Sequence.SeqList[seqno][x] for x in sorted(bf.FlattenList(idxlist3p))]) + SubsHP[0:max(0,5-HPLen)]
            #Sub3p = "".join([Sequence.SeqList[seqno][x] for x in sorted(bf.FlattenList(idxlist3p))] + list(reversed(SubsHP))[0:max(0, 5 - HPLen)])
            #print("Sub5p", Sub5p)
            #print("Sub3p", Sub3p)
            return Sub5p, Sub3p
        else:
            i = self.FG.HP[HPID].idx[0]
            SubsHPtemp =  None
            if seqno == "master_numbering":
                #print(self.FG.Sequence[seqno][i[0]:i[0] + HPLen])
                #print(SubsHP[min(HPLen, self.subHPlength):-min(HPLen, self.subHPlength)])
                #print(self.FG.Sequence[seqno][i[1] - HPLen + 1:i[1] + 1])
                SubsHPtemp = [it for it in self.FG.Sequence[seqno][i[0]:i[0] + HPLen]]
                SubsHPtemp.extend(SubsHP[min(HPLen, self.subHPlength):-min(HPLen, self.subHPlength)])
                SubsHPtemp.extend(self.FG.Sequence[seqno][i[1] - HPLen + 1:i[1] + 1])
                #print(SubsHPtemp)
            else:
                SubsHPtemp = self.FG.Sequence[seqno][i[0]:i[0] + HPLen] + SubsHP[min(HPLen, self.subHPlength):-min(HPLen, self.subHPlength)] + self.FG.Sequence[seqno][i[1] - HPLen + 1:i[1] + 1]
            return SubsHPtemp
    def _getSubstituteParentHP(self,seqno,HPID):
        SubsHP = input.SubsHP[seqno]
        HPLen = self.FG.HP[HPID].length
        idxlist = self.FG.HP[HPID].idx
        idxlist5p = [x[0] for x in self.FG.HP[HPID].idx]
        idxlist3p = [x[1] for x in self.FG.HP[HPID].idx]
        Sub5p = None
        if seqno == "master_numbering":
            #Sub5p = SubsHP[0:max(0, self.subHPlength - HPLen)].extend(bf.FlattenList([self.FG.Sequence.SeqList[seqno][x] for x in sorted(bf.FlattenList(idxlist5p))]))
            Sub5p = [it for it in SubsHP[0:max(0, self.subHPlength - HPLen)]]
            Sub5p.extend(bf.FlattenList([self.FG.Sequence.SeqList[seqno][x] for x in sorted(bf.FlattenList(idxlist5p))]))
        else:
            Sub5p = SubsHP[0:max(0, self.subHPlength - HPLen)] + "".join([self.FG.Sequence.SeqList[seqno][x] for x in sorted(bf.FlattenList(idxlist5p))])
        Sub3p = None
        if seqno == "master_numbering":
            Sub3p = [it for it in [self.FG.Sequence.SeqList[seqno][x] for x in sorted(bf.FlattenList(idxlist3p))]]
            Sub3p.extend(bf.FlattenList(list(reversed(SubsHP))[0:max(0,self.subHPlength-HPLen)]))
        else:
            Sub3p = "".join([self.FG.Sequence.SeqList[seqno][x] for x in sorted(bf.FlattenList(idxlist3p))+list(reversed(SubsHP))[0:max(0,self.subHPlength-HPLen)]])
            
        if seqno == "master_numbering":
            Sub3p = [it for it in [self.FG.Sequence.SeqList[seqno][x] for x in sorted(bf.FlattenList(idxlist3p))]]
            Sub3p.extend(bf.FlattenList((SubsHP[min(0,-self.subHPlength+HPLen):len(SubsHP)])))
        else:
            Sub3p = "".join([self.FG.Sequence.SeqList[seqno][x] for x in sorted(bf.FlattenList(idxlist3p))]) + SubsHP[min(0,-self.subHPlength+HPLen):len(SubsHP)]

        return Sub5p,Sub3p
    def _getSeq(self,seqno): #Logic for what to substitute for which kinds of bases (unchanged, hp->bulge (ShellHPs), bulge/other_hp ->hp (SubbedHPs)
        ShellHPs = list(self.StapleHPs)
        BulgeHPs = list(self.ShellHPs)
        ParentHP = self.ParentHP
        self.NativeIdx = list(sorted(bf.FlattenList(self._getHPIdx() + self._getBulgeIdx())))
        SubbedHPIdx = bf.FlattenList([self.FG.HP[x].idx[0] for x in BulgeHPs])

        AllIdx = list(range(0,len(self.FG.Sequence)))
        AllSeq = ['' for x in AllIdx]

        #Inserts the native bases
        for idx in self.NativeIdx:
            AllSeq[idx] = self.FG.Sequence[seqno][idx]

        #Inserts the bases for all HPs that are being subbed in
        for idx in SubbedHPIdx:
            AllSeq[idx] = input.GenSub[seqno] #This defines what to replace hairpins with when going from a HP->bulge base transition

        #Inserts the bases for all HPs that are being replaced
        for HPID in ShellHPs:
            HPIdx = self.FG.HP[HPID].idx[0]
            FirstIdx = HPIdx[0]
            if HPID==ParentHP and any(elem in self.HPList for elem in self.FG.HP[ParentHP].progenyID):
                LastIdx = HPIdx[1]
                Sub5p,Sub3p = self._getSubstituteHP(seqno,HPID,'in')
                AllSeq[FirstIdx] = Sub5p 
                AllSeq[LastIdx] = Sub3p
            else:
                AllSeq[FirstIdx] = self._getSubstituteHP(seqno,HPID,'out')

        if seqno == "master_numbering":
            return [it for it in bf.FlattenList(AllSeq) if it != '']
        else:
            return "".join(AllSeq)

    def printSeq(self):

        print(self._getSeq('sequence'))
        print(self._getSeq('secstruct'))
        print(self._getSeq('lock'))
        print(self._getSeq('iupac'))
        print(self._getSeq('master_numbering'))

        return

    def __len__(self):
        return len(self._getSeq('sequence'))
    def __iter__(self):
        self.index = len(self.HPList)
        return self
    def __next__(self):
        if self.index==0:
            #print("Iteration Done")
            raise StopIteration
        self.index = self.index-1
        return self.HPList[self.index]

class PartialPuzzle():

    def __init__(self,FG,maxlength,minsubHPlength, subHPlength, RequiredHPs):
        self.FG = FG
        self.maxlength = maxlength
        assert(minsubHPlength is not None)
        self.minsubHPlength = minsubHPlength
        self.subHPlength = subHPlength
        assert(self.minsubHPlength is not None)
        #self.RequiredHPs = self.FG.HP._IDlist
        self.RequiredHPs = list(RequiredHPs)
        self.OverWriteCandidateFlag = 0
        #HairpinList.__init__(self,FG,[],self.maxlength)
        self.HPList = [] #Initializes as an empty list, but in initHPlist is
        self._updateIHPList() #IHP - inverse hairpinlist. Essentially identifies what is remaining in the overall puzzle
        self._initHPList(self.Candidate)
        self.FailedToAddCandidate = 0
        # print("IHPList",self.IHPList)
        # print("HPList",self.HPList.HPList)
        # print("LeafHP",self.LeafHP)
        # print("Candidate",self.Candidate)
        while self.FailedToAddCandidate==0:
            if self.OverWriteCandidateFlag==0:
                self._updateIHPList('hard',1)
            self._CrawlPath()
        if self.HPList.InnerHPs !=[]:
            print("Included HPs:", self.HPList.HPList)
            print("Included Idxs:", list(bf.find_ranges(list(self.HPList.NativeIdx))))
            self.HPList.printSeq()
        else:
            print("HP does not fit within length requirements")

    def update_parent_HP(self):
        self.ParentHP = []
        for HPID in self.IHPList:
            if self.FG.HP[HPID].parentID==[] or not any(elem in self.FG.HP[HPID].parentID for elem in self.IHPList):
                self.ParentHP.append(HPID)

    def _updateIHPList(self,hardoreasy='hard',Normalize=0):
        self.IHPList = self.getIHPList()
        
        self.update_parent_HP()
        #self.ParentHP = self.getParentHP()
        
        self.LeafHP = self.getLeafHP()
        self.Candidate = self._sortFeatListByWeight('HP', list(self.LeafHP), hardoreasy, Normalize)

    def getIHPList(self):
        return [x for x in self.RequiredHPs if x not in self.HPList]
   
    def getParentHP(self):
        ParentHP = []
        for HPID in self.IHPList:
            if self.FG.HP[HPID].parentID==[] or not any(elem in self.FG.HP[HPID].parentID for elem in self.IHPList):
                ParentHP.append(HPID)
        return ParentHP
    def getLeafHP(self):
        return [x for x in self.IHPList if self.FG.HP[x].isleaf=='yes' or not any(elem in self.FG.HP[x].progenyID for elem in self.IHPList)]
    def _updateCandidates(self):
        self.IHPList = self.getIHPList()
        update_parent_HP(self)
        self.LeafHP = self.getLeafHP()
        try:
            self.LastAddedHP = self.HPList.lastAddedHP
            NeighborIDs = [self.FG.HP[x].neighborID for x in self.HPList]
            self.Candidate = self._sortFeatListByWeight('HP', list(self.LeafHP), 'easy', 1)
            self.Candidate = [x for x in NeighborIDs if x in self.Candidate]
        except:
            self.Candidate = self._sortFeatListByWeight('HP', list(self.LeafHP), 'hard', 0)
        return self

    def _sortFeatListByWeight(self,FeatName,FeatList,EasyOrHard='hard',NormalizeFlag=0):
        SortedIDs = {}

        Feat = getattr(self.FG, FeatName)
        Weights = [Feat[x].Weight for x in FeatList]
        Lengths = [Feat[x].length for x in FeatList]
        if NormalizeFlag==1:
            Weights = [x/(y**(0/5)) for x,y in zip(Weights,Lengths)]
        if EasyOrHard=='hard':
            SortedIDs = list(reversed([y for x,y in sorted(zip(Weights,FeatList))]))
        else:
            SortedIDs = [y for x,y in sorted(zip(Weights,FeatList))]
        return SortedIDs

    def _initHPList(self,CandidateList):
        for HPID in CandidateList:
            self.HPList = HairpinList(FG = self.FG, HPList=[HPID], minsubHPlength=self.minsubHPlength, subHPlength=self.subHPlength, maxlength=self.maxlength)
            if self.HPList.AddedLast=='yes':
                self._updateIHPList()
                return

    def _CrawlPath(self):
        for i in self.Candidate:
            self.HPList._addHPList([i])
            if self.HPList.AddedLast=='yes':
                self.OverWriteCandidateFlag = 0
                return
        if len(self.HPList)<0.6*self.maxlength and self.OverWriteCandidateFlag==0:
            self.Candidate += [self.FG.HP[x].neighborID for x in self.HPList]
            self.OverWriteCandidateFlag = 1
            self.FailedToAddCandidate = 0
        else:
            self.HPList.maxlength = 1.1*self.HPList.maxlength
            for i in self.HPList.StapleHPs: #Attempt to add any remainders (i.e. replace staple HPs with
                self.HPList._addHPList(list(self.FG.HP[i].progenyID))
                self.HPList._addHPList(list(self.FG.HP[i].parentID))
            self.FailedToAddCandidate = 1

#maxlength = 300
#minsubHPlength = 3 #HPs <=3 will be substituted
#subHPlength = 6  #will be replaced with HPs of len <


def main():
    parser = argparse.ArgumentParser(description='Make subpuzzles')
    parser.add_argument('--maxlength', type=int, help='max puzzle length', default=300)
    parser.add_argument('--minsubHPlength', type=int, help='min length of a hairpin below which we substitute', default=3)
    parser.add_argument('--subHPlength', type=int, help='length of sub hairpin', default=6)
    parser.add_argument('--spec', type=str, help='file with puzzle specification', default='16s.puzzle_spec')

    args = parser.parse_args()

    PS_Goal, SS_Goal, LS_Goal, IUPAC_Goal = None, None, None, None
    with open(args.spec) as f:
        for line in f:
            if line[0] == '#': continue
            if line.split()[0] == 'PS_Goal:':
                PS_Goal = line.split()[1].strip().replace("'", "")
            if line.split()[0] == 'SS_Goal:':
                SS_Goal = line.split()[1].strip().replace("'", "")
            if line.split()[0] == 'LS_Goal:':
                LS_Goal = line.split()[1].strip().replace("'", "")
            if line.split()[0] == 'IUPAC_Goal:':
                IUPAC_Goal = line.split()[1].strip().replace("'", "")

    FG = FD.FeatureGroup(PS_Goal, SS_Goal, LS_Goal, IUPAC_Goal, [str(ii+1) for ii in range(len(PS_Goal))])
    count = 1
    print("Puzzle #", count)
    IHP = PartialPuzzle(FG, minsubHPlength=args.minsubHPlength, subHPlength=args.subHPlength, maxlength=args.maxlength,RequiredHPs=FG.HP._IDlist)
    RequiredHPs = [x for x in IHP.RequiredHPs if x not in IHP.HPList.HPList]
    RequiredHPs = [x for x in RequiredHPs if len(FG.HP[x].childrenID) < 5]

    while RequiredHPs!=[]:
        count+=1
        print("Puzzle #", count)
        IHP = PartialPuzzle(FG, minsubHPlength=args.minsubHPlength, subHPlength=args.subHPlength, maxlength=args.maxlength, RequiredHPs=RequiredHPs)
        RequiredHPs = [x for x in IHP.RequiredHPs if x not in IHP.HPList.HPList]
        if RequiredHPs==IHP.RequiredHPs:
            break

    print(PS_Goal)
    print(SS_Goal)


    # print(bf.getParentHP(FG,[1,2,33]))
    # print(bf.getAllHPsConnectingList(FG,[1,2,33]))
    # print(HPList._getAllConnectingHPs([1,2,33]))


main()
