''' This script is written by Angad Garg. This script allows strand specific identification of annotations
    with increased, proximal or distal polyA sites from RNA-seq data.
    This is compatible with Illumina TruSeq stranded RNA-Seq library preparation
    and determination of poly(A) peaks from HOMER as described in Schwer et al 2021. (can be modified based on end-user needs) 
    Input Files: Peak (putative polyA site); co-ordinate file (tab-delimited text; format given below); and
    GFF file (tab delimited text; format given below) '''


## use in current directory with all accessory files in the same folder

## Add Peak File (ammend the following text below - your_filename)
## [format: .txt (tab-delimited) with no header, column1: peakID; column2: Chromosome number; column3: peak Start coordinate; column4: peak End coordinate; column5: Strand (+/-); can add your favorite columns after this) 

polyA_data = open ('your_filename.txt', 'r')
polyADataMemory = polyA_data.readlines()
polyA_data.close()


## Add feature file (ammend the following text below - your_gff_file)
## [format: .txt (tab-delimited) file with no header, column1: Chromosome number (format same as peak file); column2: Feature start coordinate; column3: Feature end coordinate; column4: Strand (+/-); can add your favorite columns after this]

chr_GFF = open ('your_gff_file.txt', 'r')
chr_GFF_data = chr_GFF.readlines()
chr_GFF.close()

polyAPeak = []
polyAData = []
chrGFF = []


for line in polyADataMemory:
    line = line.rstrip('\n')
    polyAData.append(line.split('\t'))


for line in chr_GFF_data:
    line = line.rstrip('\n').rstrip(' ')
    chrGFF.append(line.split('\t'))




## File name of resulting file(s); ammend output filenames below

match_file = open('your_annotations_upstream_of_peaks.txt', 'w')
match_file2 = open('peaks_within_annotations.txt', 'w')

for line in polyAData:
    polyAPeak =((int(line[2])+int(line[3]))/2)
    
    for featLine in chrGFF:
        if line[1] == featLine[0]:
            ## library prep stranded information is antisense to the features, therefore '+' strand of peak determined from HOMER is '-' strand of feature 
            if line[4] == '+':
                if featLine[3] == '-':
                    if polyAPeak >= int(featLine[1]):
                        if polyAPeak <= int(featLine[2]):
                            match_file2.write('\t'.join(line) + '\t' + '\t'.join(featLine) + '\n')
                    if polyAPeak < int(featLine[1]):
                        #  int(featLine[1]) - polyAPeak <= number can be varied to determine how far upstream an annotation is (current is 500 nucleotide) 
                        if int(featLine[1]) - polyAPeak <= 500:
                            match_file.write('\t'.join(line) + '\t' + '\t'.join(featLine) + '\n')
                        
            else:
                ## library prep stranded information is antisense to the features,  therefore '-' strand of peak determined from HOMER is '+' strand of feature
                if line[4] == '-':
                    if featLine[3] == '+':
                        if polyAPeak >= int(featLine[1]):
                            if polyAPeak <= int(featLine[2]):
                                match_file2.write('\t'.join(line) + '\t' + '\t'.join(featLine) + '\n')
                        if polyAPeak > int(featLine[2]):
                            #  polyAPeak - int(featLine[2]) <= number can be varied to determine how far upstream an annotation is (current is 500 nucleotide)
                            if polyAPeak - int(featLine[2]) <= 500:
                                match_file.write('\t'.join(line) + '\t' + '\t'.join(featLine) + '\n')

match_file.close()
match_file2.close()
