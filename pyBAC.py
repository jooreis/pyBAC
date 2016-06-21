#pyBAC is a filtering pipeline that reports breakpoints in MALBAC generated sequence data.
#This pipeline requires the pysam module and Sambamba. In addition, you have to provide a BED-file with bin coordinates
#At the end of this script you can tune the settings. pyBAC will give a BAM-file for every filtering step and a text-file as output with called breakpoints.
#Run time is dependent on BAM-file size, On our files it took between 15 and 30 minutes each.

import os
import pysam

def getdiscordant(samplenumber):
    bam=pysam.Samfile("/data/MALBAC/160311_NS500813_0100_AHTCNHBGXX/IPS-chromothripsis-0" + str(samplenumber) + "/mapping/IPS-chromothripsis-0" + str(samplenumber) + "_dedup.realigned.bam",'rb') #open BAM-file
    bam_out = pysam.Samfile("filter" + str(samplenumber) + "-1.bam", 'wb', template=bam) #create output file 
    total=0
    discordant=0
    filt=0
    for read in bam.fetch(): #loop through bamfile
        total+=1 #count each read to give total number of reads 
        if read.reference_name != read.next_reference_name or abs(read.template_length)>20000: #Only take reads that have discordant mate or indel bigger than 20.000
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
    bam.close()
    bam_out.close() #Close files
    os.system("./sambamba_v0.6.0 sort filter" + str(samplenumber) + "-1.bam") #Sort and index file
    
#--------------------------------------------------------------------------------------------------

def gethomogenousbins1(samplenumber):
    bam=pysam.Samfile("filter" + str(samplenumber) + "-1.sorted.bam",'rb') #Open BAM-file that only has discordant reads of high quality
    breaks=open("bin.txt", "r") #Open list of bin locations coordinates, file generated with sambamba; ./sambamba_v0.6.0 depth window -w 500 -t 4 bam > bin.txt
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
            frequent = dict((k, v) for k, v in diction.items() if v >= 5) #Create new dictionary with all mate chromosome number that occur more than 5 times
            if 1 <= len(frequent) <= 2: #Select bins that have between 1 and 2 different discordant locations, a measure for homozygosity 
                frequencies = frequent.keys() #Only take the chromosome number, not the value i.e. count
                for frequency in frequencies: 
                    for bin in binned: #only take the reads that have a mate chromosome number that passed the filter in line 28 (more than 2). 
                        if bin.next_reference_name == frequency:
                            bam_out.write(bin) #write to new bam file
                            filtered+=1
    print "Number of reads after filtering: " + str(filtered)
    bam.close()
    bam_out.close()
    os.system("./sambamba_v0.6.0 sort filter" + str(samplenumber) + "-2.bam") #Sort and index file
    
#---------------------------------------------------------------------------------------------------

def removedup1(samplenumber):
    bam=pysam.Samfile("filter" + str(samplenumber) + "-2.sorted.bam",'rb') #Load input file
    bam_out = pysam.Samfile("filter" + str(samplenumber) + "-3.bam", 'wb', template=bam) #Create output file 
    seen = []
    uniq = [] #Create empty lists that will track all unique names in the BAM-file
    position=[]
    uniqpos=[]
    for read in bam.fetch(): #Loop over each read in BAM-file
        if read.query_name not in seen or read.reference_start not in position: #Only keep unique reads, as indentified by their name or position (as reads with the same name but different position are not duplicates)
            uniq.append(read.query_name) #Append read info to list 
            uniqpos.append(read.reference_start)
            bam_out.write(read) #Add read to new BAM-file
        seen.append(read.query_name) #Add read to seen list, to know which are double
        position.append(read.reference_start)        
    print "Number of reads after removing dupplicates: " + str(len(uniq)) #Print number of unique reads 
    bam.close()
    bam_out.close()
    os.system("./sambamba_v0.6.0 sort filter" + str(samplenumber) + "-3.bam") #Sort and index file
    
#---------------------------------------------------------------------------------------------------

def removesingelton(samplenumber): 
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
    print "Number of reads after filtering out single reads: " + str(counter)
    bam.close()
    bam_out.close()
    os.system("./sambamba_v0.6.0 sort filter" + str(samplenumber) + "-4.bam") #Sort and index file
    
#--------------------------------------------------------------------------------------------------------

def gethomogenousbins2(samplenumber):
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
    print "Number of reads after filtering: " + str(filtered)
    bam.close()
    bam_out.close()
    os.system("./sambamba_v0.6.0 sort filter" + str(samplenumber) + "-5.bam") #Sort and index file
    
#------------------------------------------------------------------------------------------------------

def removedup2(samplenumber):
    bam=pysam.Samfile("filter" + str(samplenumber) + "-5.sorted.bam",'rb') #Load input file
    bam_out = pysam.Samfile("filter" + str(samplenumber) + "-6.bam", 'wb', template=bam) #Create output file 
    seen = []
    uniq = []#Create empty lists that will track all unique names in the BAM-file
    position=[]
    uniqpos=[]
    for read in bam.fetch(): #Loop over each read in BAM-file
        if read.query_name not in seen or read.reference_start not in position: #Only keep unique reads, as indentified by their name or position (as reads with the same name but different position are not duplicates)
            uniq.append(read.query_name) #Append read info to list 
            uniqpos.append(read.reference_start)
            bam_out.write(read) #Add read to new BAM-file
        seen.append(read.query_name) #Add read to seen list, to know which are double
        position.append(read.reference_start)          
    print "Number of reads after removing dupplicates: " + str(len(uniq)) #Print number of unique reads 
    bam.close()
    bam_out.close()
    os.system("./sambamba_v0.6.0 sort filter" + str(samplenumber) + "-6.bam") #Sort and index file
    
#------------------------------------------------------------------------------------------------------

def getoutput():
    breaks=open("bin.txt", "r") #Open list of bin locations
    temp=pysam.Samfile("filter1-6.sorted.bam",'rb') #Load template file
    bam_out = pysam.Samfile("fullout.bam", 'wb', template=temp) #create output file to combine breakpoints
    lines=breaks.read()
    samples=lines.split("\n") #is list of bin coordinates
    for i in range(1,9): #Loop over each BAM-File
        bam=pysam.Samfile("filter" + str(i) + "-6.sorted.bam",'rb') #Load input BAM-file for each sample     
        for sample in samples[1:-1]: #Loop through each bin 
            location=sample.split("\t") 
            binned=[]    #for each bin coordinate make an empty bin
            for read in bam.fetch(location[0], (int(location[1])), (int(location[2]))):  #give all reads in bin coordinate location
                binned.append(read) #add that read to bin 
            if len(binned)!=0:
                bam_out.write(read) #add bin to BAM-file
        bam.close() #Close current BAM-file to allow next sample to be processed 
    breaks.close()
    bam_out.close() #Close all used files
    os.system("./sambamba_v0.6.0 sort fullout.bam") #Sort and index BAM-file
    
#--------------------------------------------------------------------------------------------------------

def getoverlap():
    breaks=open("bin.txt", "r") #Open list of bin locations to define overlap locations
    bam = pysam.Samfile("fullout.sorted.bam", 'rb') #Open BAM-file with all final reads
    out = open("fullout.txt", 'w') #Create output text file
    out.write("chrom\tstart\tnext chrom\tnext start\tsample\n") #Create header in text file
    lines=breaks.read()
    samples=lines.split("\n") #Make list of bin coordinates
    for sample in samples[1:-1]: #Loop through each bin 
        location=sample.split("\t") 
        binned=[]                       #for each bin coordinate make an empty bin
        for read in bam.fetch(location[0], (int(location[1])), (int(location[2]))):  #give all reads in bin coordinate location
            binned.append(read) #add that read to bin 
        if len(binned)!=0: #Only select bins that contain reads 
            diction={} #Make dictionary to count number of samples per bin
            for i in binned: #Loop over each sample in bin
               for j in i.tags: #Take tag from each read (tag contains sample name)
                    for n in j: #Order of tag atributes differs from read to read so this--with line 18--looks for sample name
                        if str(n).find("IPS") >=0: 
                            if n[20] in diction: #Count the number of different samples in each bin
                                diction[n[20]]+=1
                            else:
                                diction[n[20]]=1
                            
            if len(diction) >=4: #Only select bins that have more than 4 different samples i.e. Breakpoint is seen in more than 4 cells
                out.write(str(binned[0].reference_name) + "\t") #Write different atributes of breakpoint to text file 
                out.write(str(binned[0].reference_start) + "\t")
                out.write(str(binned[0].next_reference_name) + "\t")
                out.write(str(binned[0].next_reference_start) + "\t")
                for k in diction:
                    out.write(str(k) + "\t") #Report by which samples the breakpoints are covered 
                out.write("\n")
    bam.close()
    out.close()
    print "Done"
    
#-------------------------------------------------------------------------------------------------------

def main(sample):
    getdiscordant(sample)
    gethomogenousbins1(sample)
    removedup1(sample)
    removesingelton(sample)
    gethomogenousbins2(sample)
    removedup2(sample)
    
def final():
    getoutput()
    getoverlap()
    
for i in range(1,9):
    main(i)
final()