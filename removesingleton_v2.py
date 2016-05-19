def removesingelton(samplenumber):
    import pysam
    import os 
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