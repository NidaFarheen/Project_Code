#This file contains all the user defined modules
import random
import math
import collections
import sys
import numpy
from udf import *
from string import maketrans

color_space={'AA':'0', 'CC':'0', 'GG':'0', 'TT':'0', 'AC':'1', 'CA':'1', 'GT':'1', 'TG':'1', 'AG':'2', 'CT':'2', 'GA':'2', 'TC':'2', 'AT':'3', 'CG':'3', 'GC':3, 'TA': 3} #For Solid Platform Color Space Coding

revcomp=maketrans('AGTC', 'TCAG') #Used for finding reverse complement

def substitute(base, l):
    
    '''l is a list -for eg.  ['A', 'G', 'T', 'C']
    Base : One of list elements
    Function returns any random element from the list excluding the base'''

    l.remove(base) #Remove that base
    i=random.randrange(0, 3) #Choose any other base
    return(l[i])


def readg(seq, readlen, genome_type,type1,insertsize=500):

    '''Seq : Chromosome sequence
    Readlen : Readlength
    Genome_Type :Linear/Circular
    Type1 determines SIngle/Paired End Sequencing
    Function chooses a random position in the sequence and generates a read of the required readlength starting from the chosen position
    Returns read and the location of read'''
    
    if (type1==1):
        if(genome_type=='l'): #Generate Read1
            loc=random.randrange(0, len(seq)-readlen+1) #Random Location selected (loc)
            read=seq[loc:loc+readlen]  #Read Generation
        else: #Circular Genome
             loc=random.randrange(0, len(seq))
             if(loc>len(seq)-readlen):
                r1=seq[loc:]
                r2=seq[:readlen-len(r1)]
                read=r1+r2
             else:
                read=seq[loc:loc+readlen]

        return read, loc

    if(type1==2):
        if (genome_type=='l'):
            loc=random.randrange(0,(len(seq)-insertsize)+1)
            read3_5=seq[loc:loc+readlen]
            read5_3=seq[(loc+insertsize)-1:(loc+insertsize)-readlen+1:-1]
            read5_3=read5_3.translate(revcomp)  #Reverse complement

        if(genome_type=='c'):
            loc=random.randrange(0,len(seq))
            if (loc>len(seq)-insertsize) : #Case1 

                if(loc>len(seq)-readlen): #Read1 goes out of loop 1 range #Subcase
                    r1=seq[loc:] #Region in loop 1
                    r2=seq[:readlen-len(r1)] #Region in loop 2
                    p1=insertsize-len(r1) #Determine the end position of fragment
                    read3_5=r1+r2 #generated first read 
                    read5_3=seq[p1-1:p1-readlen+1:-1]
                    read5_3=read5_3.translate(revcomp)  #Reverse complement

                if(insertsize<loc <= len(seq)-readlen): #Read 1 is well within loop 1 and read2 goes out #SubCase
                    read3_5=seq[loc:loc+readlen] #No complexity with Read3_5 generation
                    p2=loc-(len(seq)-insertsize) #This offset will be added in loop 2

                    if(p2>=readlen): #The read2 will completely lie in Loop 2
                        read5_3=seq[p2-1:p2-readlen+1:-1]
                        read5_3=read5_3.translate(revcomp) 

                    if(p2<readlen): #Half part of read2 will lie in loop 1 and half part in loop 2
                        p3=readlen-p2 #This much would be taken from  loop 1 again
                        read5_3_partA=seq[p2-1:p1-readlen+1:-1]
                        read5_3_partA=read5_3_partA.translate(revcomp) 
                        read5_3_partB=seq[len(seq)-1:len(seq)-p3+1:-1]
                        read5_3_partB=read5_3partB.translate(revcomp)
                        read5_3=read5_3_partA + read5_3_partB


            else: #Case 2 - Both reads are well in loop 1 - Normal Case
                 read3_5=seq[loc:loc+readlen]
                 read5_3=seq[(loc+insertsize)-1:(loc+insertize)-readlen+1:-1]
                 read5_3=read5_3.translate(revcomp)  #Reverse complement

        
        return read3_5, read5_3, loc


def illumina_errorprofile(read, platform):
    '''Inputs are read, platform (This function is used for both SOLiD and Illumina profile) 
    Function returns modified read, quality, pseudo_cigar_format. 
    Quality of reads follow illumina error profile'''

    qual='' # Quality
    sub=0 
    match=0 
    for j in range(len(read)):
        y1=35.0/(1+math.exp(0.6*(j-(len(read))))) # A variant of Logistic function used
        stdev=0.5+(j/10.0) #Error Bars for each position  
        y2=random.normalvariate(0, stdev)
        q=y1+y2 
        while (q<0 and q>40):
            y2=random.normalvariate(0, stdev)
            q=y1+y2   #Final Quality

        prob=pc(q) #Calling probability calculator

        r=random.random()   #Generate a random number between 0 and 1
        if(r<prob):    #Poor Quality!
            if(platform==1): #Illumina
                newbase=substitute(read[j], ['A', 'G', 'T', 'C'])
            else:  #If Platform==SOLiD
                newbase=substitute(read[j], ['0', '1', '2', '3'])

            read=read[0:j]+newbase+read[j+1:] #Read Modified
            sub=sub+1
        else:
            match=match+1
        qual=qual+chr(int(round(q))+33) #Offset of 33+ ascii

    pseudo_cigar_format=str(sub)+'S'+'0IOD'+str(match)+'M'
    return [read, qual, pseudo_cigar_format]


def error(option):
    """If meaningless command is passed, then error() function is called. The input is the command line option which was not acceptable."""
    
    print "ERROR Found\nExecution stopped\nValue assigned to option -", option, "is wrong"
    if(option=='a'):
        print('Type -a 1 for Single End Sequencing\nType -a <2,insertsize>  :  Paired End Sequencing. Eg. -a 2,500 implies that user has chosen Paired End Sequencing with a insert size of 500')
    if(option=='p'):
        print('Type -p 1 for Illumina\nType -p 2 for SOLiD\nType -p 3 for PacBio RS')
    if(option=='s'): 
        print('Enter either Coverage or Number of Reads\nFor entering coverage type -s c<insert coverage>\nFor entering Number of Reads type -s n<insert number of reads>')
    if(option=='r'):
        print('For entering read length, type -r <insert readlength>\nIf platform is PacBio RS, type -r mean,standard_dev')
    if(option=='x' or option=='y' or option=='z'):
        print('For single end sequencing, to enter the output file path,  type -x <filepath>\nFor paired end sequencing, to enter output file paths, type -y <filepath1> & -z <filepath2>')
    if(option=='i'):
        print('To enter input file path, type -i <insert input filepath>')
    if(option=='g'):
        print('Type -g l for linear genome sequencing\nType -g c for circular genome sequencing')
    if(option=='v'):
        print('Type -v True for Verbose mode on\nType -v False for verbose mode off')
    
    print ('Type help for more details.')
    sys.exit()


def writer(fh, idn, read, qual):
    '''fh is file handle. The data is written on file in FASTQ Format
    idn: Chrmosome Name/Id, qual is quality of reads'''
    
    try:
        fh.write(idn)
        fh.write(read)
        fh.write('\n+\n')
        fh.write(qual)
        fh.write('\n')
    except IOError:  #If the data couldnt be written due to some issue
        print('The output file cannot be written to the directory mentioned')
        sys.exit() #Stop execution

def qualc(p):
    '''Quality Calculator:Takes the probability and generates the quality'''
    q=-10*(math.log10(p))
    return q

def pc(q):
    '''Probability Calculator: Takes the quality and generates probability'''
    p=math.pow(10, -0.1*q)
    return p

def bool_parse(s):
    '''Bool_parse takes a string s and checks if s== 'True'. If it is, 'True' is converted to boolean True.'''
    if(s=='True'):
        s=True
    elif(s=='False'):
        s=False
    return s




def fasta(file_handle):
    '''COnverts a fasta file into {chrm_name:seq} format. Returns the dictionary'''
    d={}
    d=collections.OrderedDict()
    for line in file_handle:
        line=line.strip()
        if line.startswith('>'):
            seq=''
            name=line[1:]
        else:
            seq=seq+line
            d[name]=seq
    return(d)

def colorspace(dna):
    '''Used for encoding a DNA Sequence into color space - Used for SOLiD platform'''
    cs=''
    for j in range(len(dna)-1):
        dibase=dna[j:j+2] #Dibase refers to AG and GT in 'AGT'
        cs=cs+str(color_space[dibase]) #cs : colorspace coding
    return cs


def pacb(read):
    '''Takes a read and simulates insertion, deletion and substitution errors.
    Returns the modified read, quality and pseudo_cigar_format'''

    qual=''
    changes=[] #Makes a list of changes
    offset=0
    match=0
    sub=0
    ins=0
    delete=0

    for position, letter in enumerate(read):
        q=numpy.random.normal(30, 10)  #Assign a quality to each base. Quality is taken randomly from a gaussian distribution
        while(q<0 or q>40):    #Quality range between  0-50
           q=random.gauss(30, 10)

        #pdel=30* (rdel/(rsub+rins+rdel))
        #pins=pc(q) * (rsub/(rsub+rins+rdel))
        #psub=pc(q) * (rins/(rsub+rins+rdel))
        psub=0.1
        pins=0.2
        pdel=0.3
        r=random.random()   #Generate a random number between 0 and 1
        if (0<r<=psub):  #Do substitution
            changes.append((position, "Sub"))
            qual=qual+chr(int(round(q))+33)
            sub=sub+1
        if(psub<r<=pins): #Do Insertion
            changes.append((position, "Insertion"))
            qual=qual+chr(33)+chr(int(round(q))+33)  #Keep the insert base Quality fixed
            ins=ins+1
        if(pins<r<=pdel): #Do Deletion
            changes.append((position, "Delete"))
            delete=delete+1
        if(r>pdel):
            qual=qual+chr(int(round(q))+33)
            match=match+1
 
    pseudo_cigar_format=str(sub)+'S'+str(ins)+'I'+str(delete)+'D'+str(match)+'M'

    for each_change in changes:
        pos=each_change[0]
        operation=each_change[1]
        if(operation=='sub'):
            read[pos+offset]=substitute(read[pos+offset], ['A', 'G', 'T', 'C'])
        if(operation=='delete'):
            read.remove(read[pos+offset])
            offset=offset-1
        if(operation=='insert'):
            read.insert(pos+offset, read[pos+offset])
            offset=offset+1

    read="".join(read)
    return(read, qual, pseudo_cigar_format)

def pacb_readlen(mean, stddev):
    '''Generates a random readlength from a log normal distribution.'''

    y=numpy.random.lognormal(mean, stddev)
    readlen=math.exp(y)
    readlen=int(round(readlen))
    return readlen

