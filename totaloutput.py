import pysam
breaks=open("bin.txt", "r") #Open list of bin locations
out = open("finaloutpyBAC.txt", 'w') #Create output file
out.write("chrom\tstart\tsample\n")
new=0
filtered=0
lines=breaks.read()
samples=lines.split("\n") #is list of bin coordinates
for i in range(1,9):
    print i
    bam=pysam.Samfile("filter" + str(i) + "-6.sorted.bam",'rb') #Load input file 
    for sample in samples[1:-1]: #Loop through each bin 
        location=sample.split("\t") 
        binned=[]                       #for each bin coordinate make an empty bin
        for read in bam.fetch(location[0], (int(location[1])), (int(location[2]))):  #give all reads in bin coordinate location
            new+=1
            binned.append(read) #add that read to bin 
        if len(binned)!=0:
            out.write(str(binned[0].reference_name) + "\t")
            out.write(str(binned[0].reference_start) + "\t")
            out.write(str(i)+ "\n")
    bam.close()
                