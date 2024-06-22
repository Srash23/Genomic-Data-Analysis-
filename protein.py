
def readFile(filename):
    
    file = open(filename, 'r')
    
    sequence = file.read()
    
    file.close()
    
    organized_sequence = sequence.replace('\n', '')
    
    return organized_sequence



def dnaToRna(dna, startIndex):
    
    RNA = dna[startIndex:].replace('T', 'U')

    stop_codons = ["UAA", "UAG", "UGA"]

    Codons = []
    
    for r in range(0, len(RNA), 3):
        
        if r + 3 <= len(RNA):
            
            Codon = RNA[r:r + 3]
            
            Codons.append(Codon)
            
            if Codon in stop_codons:
                
                break
    
    return Codons


def makeCodonDictionary(filename):
    import json
    
    with open(filename, 'r') as file:
        
        aa_to_Codon = json.load(file)
        
    Codon_to_aa = {}
    
    for aa, Codons in aa_to_Codon.items():
        
        for Codon in Codons:
            
            RNA_Codon = Codon.replace('T', 'U')
            
            Codon_to_aa[RNA_Codon] = aa
            
    return Codon_to_aa


def generateProtein(codons, codonD):
    
    Protein = []
 
    for c, codon in enumerate(codons):
       
        if c == 0 and codon == "AUG":
            
            Protein.append("Start")
        
        else:
            aa = codonD[codon]
            
            if aa == "Stop":
                
                Protein.append("Stop")
                
                break
            
            else:
                
                Protein.append(aa)

    return Protein


def synthesizeProteins(dnaFilename, codonFilename):
    
    DNA = readFile(dnaFilename)
    
    codonD = makeCodonDictionary(codonFilename)
    
    Proteins = []
    
    Bases_unused = 0
    
    p = 0
    
    while p < len(DNA):
        
        if DNA[p:p+3] == 'ATG':
            
            RNA = dnaToRna(DNA, p)
            
            Protein = generateProtein(RNA, codonD)
            
            Proteins.append(Protein)
            
            p += 3 * len(RNA)
            
        else:
            
          Bases_unused = Bases_unused + 1
          
          p = p + 1
          
    print('Total number of bases : ' + str(len(DNA)))
    print('Total number of unused bases : ' + str(Bases_unused))
    print('Total number of proteins synthesized : ' + str(len(Proteins)))
    
    return Proteins


def commonProteins(proteinList1, proteinList2):
   
    Unique_Proteins_List_1 = set(tuple(p) for p in proteinList1)
    
    Unique_Proteins_List_2 = set(tuple(p) for p in proteinList2)

    Set_of_Common_Proteins = Unique_Proteins_List_1.intersection(Unique_Proteins_List_2)

    List_of_Common_Proteins = [list(p) for p in Set_of_Common_Proteins]

    return List_of_Common_Proteins


def combineProteins(proteinList):
   
    Combined_protein_list = []

    for p in proteinList:
        
        Combined_protein_list.extend(p)

    return Combined_protein_list


def aminoAcidDictionary(aaList):
    
    aa_Dictionary = {}

    for a in aaList:
       
        if a in aa_Dictionary:
            
            aa_Dictionary[a] = aa_Dictionary[a] + 1
       
        else:
            
            aa_Dictionary[a] = 1

    return aa_Dictionary


def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    
    aa_List_1 = combineProteins(proteinList1)
    
    aa_List_2 = combineProteins(proteinList2)
    

    aa_count_in_List_1 = aminoAcidDictionary(aa_List_1)
    
    aa_count_in_List_2 = aminoAcidDictionary(aa_List_2)
    

    Total_aa_List_1 = len(aa_List_1)
    
    Total_aa_List_2 = len(aa_List_2)
    

    Difference_in_aa_frequency = []
    
    
    Unique_aa = set(aa_List_1 + aa_List_2) - {'Start', 'Stop'}
    
    
    for u in Unique_aa:
        
        Frequuency_of_aa_List_1 = aa_count_in_List_1.get(u, 0) / Total_aa_List_1
        
        Frequuency_of_aa_List_2 = aa_count_in_List_2.get(u, 0) / Total_aa_List_2
        
        
        if abs(Frequuency_of_aa_List_1 - Frequuency_of_aa_List_2) > cutoff:
            
            Difference_in_aa_frequency.append([u, Frequuency_of_aa_List_1, Frequuency_of_aa_List_2])
            
    
    return Difference_in_aa_frequency



def displayTextResults(commonalities, differences):
    print("--------------------------")
    
    print("List of Common Proteins:")
    
    for c in commonalities:
        
        Proteins = []
        
        for a in c:
            
            if a not in ["Start", "Stop"]:
                
                Proteins.append(a)
                
        if Proteins:
            
            Proteins = "- ".join(Proteins)        

        print("- " + Proteins)
        
    print()
    
    print("----------------------------")
    
    print()
    
    print("List of Most Different Amino Acids:")
    
    for d in differences:
        
        aa = d[0]
        
        frequency_in_aa_List_1 = round(d[1] * 100, 2)
        
        frequency_in_aa_List_2 = round(d[2] * 100, 2)
        
        print(aa + " : " + str(frequency_in_aa_List_1) + "% in Seq1, " + str(frequency_in_aa_List_2) + "% in Seq2")
        
    print("----------------------------------")



def makeAminoAcidLabels(proteinList1, proteinList2):
    
    Unique_aa = []

    for p in proteinList1:
        
        for aa in p:
            
            if aa not in Unique_aa:
                
                Unique_aa.append(aa)

    for p in proteinList2:
        
        for aa in p:
            
            if aa not in Unique_aa:
                
                Unique_aa.append(aa)

    Unique_aa.sort()

    return Unique_aa


def setupChartData(labels, proteinList):
    
    aa_in_protein_list = []
    
    for p in proteinList:
        
        for aa in p:
            
            aa_in_protein_list.append(aa)

    
    Total_aa = len(aa_in_protein_list)

    Freq_aa = []

    for l in labels:
        
        count = 0
        
        for aa in aa_in_protein_list:
            
            if aa == l :
                
                count = count + 1  

        Freq = count / Total_aa
        
        Freq_aa.append(Freq) 

    return Freq_aa


def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt

    No_of_bars_aa = len(xLabels)

    Position_bar_1 = range(No_of_bars_aa)
    
    Position_bar_2 = [b + 0.3 for b in Position_bar_1]
    
    plt.figure(figsize=(8, 4)) 

    plt.bar(Position_bar_1, freqList1, width = 0.3, label = label1, color = '#0F52BA', edgecolor = edgeList)
    
    plt.bar(Position_bar_2, freqList2, width = 0.3, label = label2, color = '#FFC000', edgecolor = edgeList)

    plt.xlabel('Amino Acids')
    
    plt.ylabel('Frequency')
    
    plt.title('Amino Acid Frequency Distribution between Ex1 and Ex2')

    plt.xticks([aa + 0.3 for aa in range(No_of_bars_aa)], xLabels, rotation = 45)

    plt.legend()

    plt.tight_layout() 
    
    plt.show()


def makeEdgeList(labels, biggestDiffs):
   
    aa_with_biggest_diff = [i[0] for i in biggestDiffs]

    Edge_list = []
    
    for l in labels:
        
        if l in aa_with_biggest_diff:
            
            Edge_list.append('black') 
            
        else:
            
            Edge_list.append('white') 

    return Edge_list

def runFullProgram():
    test.week1Tests()
    runWeek1()
    humanDNA = 'data/human_p53.txt'
    elephantDNA = 'data/elephant_p53.txt'
    codonFilename = 'data/codon_table.jason'
    
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)
    
    labels = makeAminoAcidLabels(humanProteins, elephantProteins)
    humanFreqList = setupChartData(labels, humanProteins)
    elephantFreqList = setupChartData(labels, elephantProteins)
    edgeList = makeEdgeList(labels, differences)
    createChart(labels, humanFreqList, 'Human', elephantFreqList, 'Elephant', edgeList)
       
    return
