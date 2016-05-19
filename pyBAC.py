import pysam
import os
samplenumber=4
bam=pysam.Samfile("/data/MALBAC/160311_NS500813_0100_AHTCNHBGXX/IPS-chromothripsis-0" + str(samplenumber) + "/mapping/IPS-chromothripsis-0" + str(samplenumber) + "_dedup.realigned.bam",'rb') #open BAM-file
bam_out = pysam.Samfile("filter" + str(samplenumber) + "-1.bam", 'wb', template=bam) #create output file 
total=0
discordant=0
filt=0
for read in bam.fetch(): #loop through bamfile
    total+=1 #count each read to give total number of reads 
    if read.reference_name != read.next_reference_name or abs(read.template_length)>10000: #Only take reads that have discordant mate or indel bigger than 10.000
        discordant +=1 #Count each discordant read
        if read.mapping_quality >=50 and read.query_length>=145: #Only take reads that meet certain mapping quality paramaters 
            filt+=1
            bam_out.write(read) #Write reads that meets these selection criteria to new BAM-file 
            
fraction=((discordant*100)/total) #Calculate ratio of passed reads 
update=((filt*100)/total)            
print "Total number of reads: \t" + str(total)            
print "Total number of discordant reads: \t" + str(discordant)
print "Percentage total discordant reads: \t" + str(fraction)+ "%"
print "Total number of discordant reads after quality filter: \t" + str(filt)
print "Percentage filtered discordant reads: \t" + str(update)+ "%"
os.system("./sambamba_v0.6.0 sort filter" + str(samplenumber) + "-1.bam") #Sort and index file
bam.close()
bam_out.close() #Close files

#------------------------------------------------------------------------------

bam=pysam.Samfile("filter" + str(samplenumber) + "-1.sorted.bam",'rb') #Open BAM-file that only has discordant reads of high quality
breaks=open("bin.txt", "r") #Open list of bin locations
bam_out = pysam.Samfile("filter" + str(samplenumber) + "-2.bam", 'wb', template=bam) #Create output file 
new=0
filtered=0
lines=breaks.read()
samples=lines.split("\n") #is list of bin coordinates
for sample in samples[1:-1]: #Loop through each bin 
    location=sample.split("\t") 
    binned=[]                       #for each bin coordinate make an empty bin
    for read in bam.fetch(location[0], (int(location[1])), (int(location[2]))):  #give all reads in bin coordinate location
        new+=1
        binned.append(read) #add that read to bin 
    if len(binned)!=0: 
        mate=[]  #if bin is not empty make new list to give all mate chromosomes
        for bin in binned:
            mate.append(bin.next_reference_name) #add mate chromosome number
        diction={}
        for reference in mate:
            if reference in diction:
                diction[reference]+=1
            else:
                diction[reference]=1        #count each chromosome number per bin
        frequent = dict((k, v) for k, v in diction.items() if v >= 2) #Create new dictionary with all mate chromosome number that occur more than 2 times
        if len(frequent) ==1 or len(frequent) ==2: #and not longer than 2
            frequencies = frequent.keys() #Only take the chromosome number, not the value i.e. count
            for frequency in frequencies: 
                for bin in binned: #only take the reads that have a mate chromosome number that is frequent (more than 2). 
                    if bin.next_reference_name == frequency:
                        bam_out.write(bin) #write to new bam file
                        filtered+=1
print "Total number of reads: " + str(new)
print "Number of reads after filtering: " + str(filtered)
os.system("./sambamba_v0.6.0 sort filter" + str(samplenumber) + "-2.bam") #Sort and index file
bam.close()
bam_out.close()

#---------------------------------------------------------------------------

bam=pysam.Samfile("filter" + str(samplenumber) + "-2.sorted.bam",'rb') #Load input file
bam_out = pysam.Samfile("filter" + str(samplenumber) + "-3.bam", 'wb', template=bam) #Create output file 
seen = []    #Create empty lists to keep track of unique reads
uniq = []
position=[] 
uniqpos=[]
for read in bam.fetch(): #Loop through each read 
    if read.query_name not in seen or read.reference_start not in position: #Check list for unique name and position
        uniq.append(read.query_name) #Count unique reads 
        uniqpos.append(read.reference_start)
        bam_out.write(read) #Add unique read to new BAM-file 
    seen.append(read.query_name) #Add name and position to list of seen items, so it won't be added again
    position.append(read.reference_start)
print "Number of reads before filtering: " + str(len(seen)) #Print stats 
print "Number of unique reads after filtering: " + str(len(uniq))
os.system("./sambamba_v0.6.0 sort filter" + str(samplenumber) + "-3.bam") #Sort and index file
bam.close()
bam_out.close()

#-----------------------------------------------------------------------------

bam=pysam.Samfile("filter" + str(samplenumber) + "-3.sorted.bam",'rb') #Load input file 
bam_out = pysam.Samfile("filter" + str(samplenumber) + "-4.bam", 'wb', template=bam) #Create output file
names={} #Keep list of names to match to mate 
counter=0
total=0
for read in bam.fetch(): #Loop over all reads in file
    total+=1
    if not read.is_duplicate: #Filter out PCR duplicates 
        name=read.query_name #Make dictionary with all names
        if not name in names:
            names[name] = 1
        else:
            names[name] += 1

for item in bam.fetch(): #Check each read
    if item.query_name in names and names[item.query_name] > 1: #Only take reads that have a name with more than 1 read
        counter+=1
        bam_out.write(item)
os.system("./sambamba_v0.6.0 sort filter" + str(samplenumber) + "-4.bam") #Sort and index file
print "Total number of reads: " + str(total)
print "Number of reads after filtering out single reads: " + str(counter)
bam.close()
bam_out.close()

#-------------------------------------------------------------------------------

bam=pysam.Samfile("filter" + str(samplenumber) + "-4.sorted.bam",'rb') #Open BAM-file that only has discordant reads of high quality
breaks=open("bin.txt", "r") #Open list of bin locations
bam_out = pysam.Samfile("filter" + str(samplenumber) + "-5.bam", 'wb', template=bam) #Create output file 
new=0
filtered=0
lines=breaks.read()
samples=lines.split("\n") #is list of bin coordinates
for sample in samples[1:-1]: #Loop through each bin 
    location=sample.split("\t") 
    binned=[]                       #for each bin coordinate make an empty bin
    for read in bam.fetch(location[0], (int(location[1])), (int(location[2]))):  #give all reads in bin coordinate location
        new+=1
        binned.append(read) #add that read to bin 
    if len(binned)!=0: 
        mate=[]  #if bin is not empty make new list to give all mate chromosomes
        for bin in binned:
            mate.append(bin.next_reference_name) #add mate chromosome number
        diction={}
        for reference in mate:
            if reference in diction:
                diction[reference]+=1
            else:
                diction[reference]=1        #count each chromosome number per bin
        frequent = dict((k, v) for k, v in diction.items() if v >= 2) #Create new dictionary with all mate chromosome number that occur more than 2 times
        if len(frequent) ==1 or len(frequent) ==2: #and not longer than 2
            frequencies = frequent.keys() #Only take the chromosome number, not the value i.e. count
            for frequency in frequencies: 
                for bin in binned: #only take the reads that have a mate chromosome number that is frequent (more than 2). 
                    if bin.next_reference_name == frequency:
                        bam_out.write(bin) #write to new bam file
                        filtered+=1
print "Total number of reads: " + str(new)
print "Number of reads after filtering: " + str(filtered)
os.system("./sambamba_v0.6.0 sort filter" + str(samplenumber) + "-5.bam") #Sort and index file
bam.close()
bam_out.close()

#--------------------------------------------------------------------------------

bam=pysam.Samfile("filter" + str(samplenumber) + "-5.sorted.bam",'rb') #Load input file
bam_out = pysam.Samfile("filter" + str(samplenumber) + "-6.bam", 'wb', template=bam) #Create output file 
seen = []    #Create empty lists to keep track of unique reads
uniq = []
position=[] 
uniqpos=[]
for read in bam.fetch(): #Loop through each read 
    if read.query_name not in seen or read.reference_start not in position: #Check list for unique name and position
        uniq.append(read.query_name) #Count unique reads 
        uniqpos.append(read.reference_start)
        bam_out.write(read) #Add unique read to new BAM-file 
    seen.append(read.query_name) #Add name and position to list of seen items, so it won't be added again
    position.append(read.reference_start)
print "Number of reads before filtering: " + str(len(seen)) #Print stats 
print "Number of unique reads after filtering: " + str(len(uniq))
os.system("./sambamba_v0.6.0 sort filter" + str(samplenumber) + "-6.bam") #Sort and index file
bam.close()
bam_out.close()
print "Done"