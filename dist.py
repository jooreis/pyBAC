import pysam
bam=pysam.Samfile("fullfilt2-7-1.sorted.bam",'rb') #Open BAM-file that only has discordant reads of high quality
dictionary={}
count=0
for read in bam.fetch():
    count+=1
    if read.reference_name in dictionary:
        dictionary[read.reference_name]+=1
    else:
        dictionary[read.reference_name] =1
print count
print dictionary