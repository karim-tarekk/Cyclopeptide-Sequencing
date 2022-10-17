# Developer:
# Kareem Mohamed Wardany
# ---------------------------------------------------------

weightDict = {}
def Readfile():
    # This function to fill in the weight dictionary
    file = open("weight.txt")
    # read weights file and save data in dictionary
    for kj in file:
        line = kj.rstrip()
        words = line.split(" ")
        wit = int(words[1])
        amino = words[0]
        weightDict[wit] = amino

def GetWeight(letter):
    # This Function will retrun the Weight of any entered aa
    for i in weightDict :
        if weightDict[i] == letter:
            return i

def calcFrequency(l):
    # This Funtion will find that the aa is repeated how many times in the spectrum
    temdir = {}
    for i in l:
        x = 0
        for j in l:
            if i == j:
                x += 1
        temdir[i] = x
    return temdir

def checkFrequency(seq , InitDic):
    # This Funtion will check if that the aa will be found how many times in sequence
    sv = True
    seqfreq = calcFrequency(seq)
    for key, value in seqfreq.items():
        if InitDic[key] < value:
            sv = False
    return sv

def linear_spectrum(pept):
    # This function will generate a spectrum for a peptide
    seqList = [] 
    seqWeights = []
    # 2 loops to make all combinations of peptide
    for i in range(len(pept)):
        for j in range(i , len(pept)):
            seqList.append(pept[i:j+1])
    # Getting weights for spectrum
    for i in seqList:
        seqweight =0
        for j in i:
            seqweight = seqweight + GetWeight(j)
        seqWeights.append(seqweight)
    seqWeights.sort()
    return seqWeights

def isConsistent(subpeptide, spec, specDic): # Spectrum is the input spectrum
    generatedspec = linear_spectrum(subpeptide) # here will generate a spectrum for this subpeptide
    # check the generated spectrum with the original spectrum
    rtn = False
    for i in generatedspec:
        if i not in spec:
            return rtn
    # calcFrequency for the generated spectrum
    generatedSpecFreq = calcFrequency(generatedspec)
    # Check for the applicability of the frequencies with the original spectrum
    for key, value in generatedSpecFreq.items():
        if specDic[key] < value: # check if the weight of aa and subpeptide has only the same value of frequency with the original spectrum
            return rtn
    return True

def Initial_List(spectrum):
    # this function will get the all aa of the antibaiotic 
    initialList = []
    for i in spectrum:
        if i <= 186 and i != 0:
            initialList.append(weightDict[i])
    return initialList

def extend(InitList , TempList , InitDic):
    # This function will extend the items in 'TempList' with all aa in 'Initial_L' and found the compination
    # that will continue to the next iteration and add it to temp list and at the end return temp list
    rz = []
    for s in TempList:
        for g in InitList:
            _wd_ = s + g
            accepted = checkFrequency(_wd_ , InitDic) # check if any duplicate found it will be allowed or not just by checking the frequencey of each aa in the word
            # Accepts the extends that will help to final solution only others will not be added to solution
            if(accepted):
                rz.append(_wd_)
    return rz

def iterate (Initial_L, theoriticalSpectrum, spectrumDic, InitialDic):
    # The function will iterate and extands the subpeptide and generate spectrum for each time it extends and compare it with the original 
    # spectrum and check frequencey of generate spectrum with original spectrum if it has same value and then add it to consistantList that 
    # repeate the previous steps 
    temp = Initial_L.copy() # start first iteration with the aa that form antibiotic
    k = 1
    while k < len(Initial_L): 
        temp = extend(Initial_L, temp, InitialDic) # extends the peptides in temp list
        consistantList = []
        for i in temp:
            if isConsistent(i, theoriticalSpectrum , spectrumDic): # check if each peptide's spectrum equal in frequencey and value so it is accepted
                consistantList.append(i)
    # Stop when all are inconsistant
        if len(consistantList) == 0: 
            break
        temp = set(consistantList.copy())
        k += 1
    return temp

if __name__ == "__main__":

    Readfile()
    # Test case : 0 97 97 99 101 103 196 198 198 200 202 295 297 299 299 301 394 396 398 400 400 497
    iptSpec = list(input("Enter the spectrum list: ").split(" "))
    
    # Convert string to int
    iptSpec = [int(i) for i in iptSpec]

    # calcFrequency for theoritical spectrum
    iptSpecFreq = calcFrequency(iptSpec)

    # Generate Initial List
    InitialList = Initial_List(iptSpec)

    # calcFrequency for Initial List
    InitialListFreq = calcFrequency(InitialList)

    # Generate the output
    result = iterate(InitialList, iptSpec, iptSpecFreq, InitialListFreq)

    # Print the output
    print(result)