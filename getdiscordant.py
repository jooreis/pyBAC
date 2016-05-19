def getdiscordant(samplenumber):
    import pysam
    import os
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
    bam.close()
    bam_out.close() #Close files
    os.system("./sambamba_v0.6.0 sort filter" + str(samplenumber) + "-1.bam") #Sort and index file
