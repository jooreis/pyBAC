def gethomogenousbins(samplenumber):
    import pysam
    import os
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