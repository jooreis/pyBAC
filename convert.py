import pysam
bam=pysam.Samfile("filter4-2-old.sorted.bam",'rb')
bam_out = pysam.Samfile("filter4-3-old.bam", 'wb', template=bam) #Create output file 
seen = []
uniq = []
position=[]
uniqpos=[]
for read in bam.fetch():
    if read.query_name not in seen or read.reference_start not in position:
        uniq.append(read.query_name)
        uniqpos.append(read.reference_start)
        bam_out.write(read)
    seen.append(read.query_name)
    position.append(read.reference_start)        
print len(uniq)
print len(seen)

    