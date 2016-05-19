import pysam
import os
bam=pysam.Samfile("filter4-2.sorted.bam",'rb') #Load input file
bam_out = pysam.Samfile("filter4-3-1.bam", 'wb', template=bam) #Create output file 
seen = []    #Create empty lists to keep track of unique reads
uniq = []
names={}
position=[] 
uniqpos=[]
counter=0
for read in bam.fetch(): #Loop through each read 
    if read.query_name not in seen or read.reference_start not in position: #Check list for unique name and position
        uniq.append(read.query_name) #Count unique reads 
        uniqpos.append(read.reference_start)
        bam_out.write(read)
        if not read.is_duplicate: #Filter out PCR duplicates 
            name=read.query_name #Make dictionary with all names
            if not name in names:
                names[name] = 1
            else:
                names[name] += 1
    seen.append(read.query_name) #Add name and position to list of seen items, so it won't be added again
    position.append(read.reference_start)
for item in bam.fetch(): #Check each read
    if item.query_name in names and names[item.query_name] > 1: #Only take reads that have a name with more than 1 read
        counter+=1
        bam_out.write(item)
os.system("./sambamba_v0.6.0 sort filter4-4-1.bam") #Sort and index file
print "Number of reads before filtering: " + str(len(seen)) #Print stats 
print "Number of unique reads after filtering: " + str(len(uniq))
print "Number of reads that have mate: " + str(counter)
bam.close()
bam_out.close()
print "done"
