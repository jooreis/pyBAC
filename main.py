import getdiscordant
import gethomogeniousbins_1
import removedup
import removesingleton_v2
import gethomogenousbins_2
import removedup_2
import getoutput

def main(sample):
    getdiscordant.getdiscordant(sample)
    gethomogeniousbins_1.gethomogenousbins(sample)
    removedup.removedup(sample)
    removesingleton_v2.removesingelton(sample)
    gethomogenousbins_2.gethomogenousbins(sample)
    removedup_2.removedup(sample)
    getoutput.getoutput(sample)
    
for i in range(1,9):
    main(i)