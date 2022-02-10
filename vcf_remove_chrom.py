## This script takes a vcf file created by ipyrad (with unique #CHROM names for each RAD locus), removes the unique names and changes the position of each SNP so that no position
## appears twice. Such a vcf file can be converted to .bed format using plink without a problem.


## open vcf file
with open("clean_90_69.vcf", 'r') as vcf_f:
    ## create to sub-categories: the header and all lines containing SNP information
    lines = []
    header = []
    for line in vcf_f:
        ## assign all lines that belong to the header to list "header"
        if line.startswith('#'):
            
            header.append(line)
            lines.append("X") ## placeholder for header lines
        ## assign all lines that contain SNP information to list "lines"    
        else:
            lines.append(line)

## create position counters 
prev_line = 0
curr_line = 0
new_line = 0
pos = []
rest = []

## There are 11 lines beginning with '#' in a vcf file created by ipyrad. The 11th line is the header of the SNP table.
for i in range(0,len(lines)):
  ## the first 11 lines do not need to be manipulated 
    if i <11:
        pos.append(0)
    ## the 12th line (python starts counting at 0) is the first line containing information about SNPs (Chromosome name, position etc.)    
    elif i == 11:
        ## this is the position argument of the first relevant line
        ## Note: vcf file entries are separated by tabulators ("\t")
        curr_line = lines[i].split("\t")[1]
        
        ## this stores all other information of this line in the list "rest"
        rest.append(lines[i].split("\t")[2:])
        
    else:
        ## current line is assigned to new line variable
        curr_line = new_line
        
        ## previous line with its position argument
        prev_line = lines[i-1].split("\t")[1]
        
        ## the position of the SNP on this line is the sum of its current position plus the position of the SNP preceding it
        new_line = int(curr_line) + int(prev_line)
        
        ## add the position of this line to the "pos" list
        pos.append(new_line)
        
        ## this stores all other information of this line in the list "rest"
        rest.append(lines[i].split("\t")[2:])
        
## Now, we have a list "rest" with all information that is not chromosome name and position for each SNP, and a list "pos" 
## that stores the updated position of all SNPs after removing unique RAD locus names.

## The vcf file in question has 572,000 lines. It is possible to assign an arbitrary number of artifical chromosomes to all the SNPs, but this is not necessary. 
#import itertools
#lst = range(1,1000)
#ls = list(itertools.chain.from_iterable(itertools.repeat(x, 573) for x in lst))

## Open your vcf file
## Here, we reassemble what we have broken up to create a new vcf file that can be read by plink
new_vcf_file = open("../../../../clean_90_69_python.vcf", "a")

## The list "pos" is as long as the original vcf file had lines
for i in range(0, len(pos)):
  
    ## For the first 11 lines, just write the lines that are stored in the list "header"
    if i < 11:
      
        new_vcf_file.write(header[i])
        
    else:
      
        ## For the additional lines, add the new artifical contigs (0s would also work), the updated positions and the rest of the vcf file back together.
        # new_vcf_file.write("contig" + str(ls[i]) + "\t" + str(pos[i]) + "\t" + '\t'.join(rest[i][0:])) ## for an arbirary number of artifical chromosomes
        new_vcf_file.write(0 + "\t" + str(pos[i]) + "\t" + '\t'.join(rest[i][0:])) ## for one artifical chromosome
        
new_vcf_file.close()
