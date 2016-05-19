def removedup(samplenumber):
    import pysam
    import os
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