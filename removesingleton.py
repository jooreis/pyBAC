import pysam
bam=pysam.Samfile("filter4-3.sorted.bam",'rb') #Load input file 
bam_out = pysam.Samfile("filter4-4.bam", 'wb', template=bam) #Create output file 
names=[] #Keep list of names to match to mate 
output=[]
for read in bam.fetch(): #Loop over all reads in file 
    if read.is_duplicate == False: #Filter out PCR duplicates 
        name=read.query_name 
        names.append(name) #Add read name to list
for nam in names: #Loop over each read name, reads with no mate will only have one read for each name, will pairs will have 2 (or more)
    counter=0
    for mate in bam.fetch(): #For each name check BAM-file  
        if mate.query_name==nam: #Match name from list to name from BAM-file 
            counter+=1
    if counter!=1: #If more than one read is given for name, put it list output 
        output.append(nam)
for out in output: #Convert the raw names in output to the full read and app to new BAM-file 
    for item in bam.fetch():
        if item.query_name==out:
            bam_out.write(item)
print "Total number of reads: " + str(len(names))
print "Number of reads after filtering out single reads: " + str(len(output))
bam.close()
bam_out.close()
print "done"