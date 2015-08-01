import random
import math
import sys
import os
import getopt
import collections
from string import maketrans
from udf import * #udf.py contains all the user defined functions


workdir=os.getcwd()  #Get the path of working directory

optlist, arg=getopt.getopt(sys.argv[1:], 'a:p:s:r:x:y:z:i:g:v:') #CommandLine Parsing

if('help' in sys.argv): #Type help to open help window
    obj=open(workdir+str('/help'), 'r')
    text=obj.read()
    print(text)
    sys.exit() #Stop Execution

config=dict(optlist)  #Convert into Dictionary {option:value} format

#Setting Default values

config.setdefault("-a", '1')  #Single End Sequencing
config.setdefault("-p", 1)  #Illumina platform
if(config['-p']=='3'): #PacBio
    config.setdefault("-r", '1.95, 0.3') #readlength follows log normal distribution. Default mean : 1.95  stddev 0.1
else:
    config.setdefault('-r', 40)  #Default Readlength for Illumina, SOLiD

config.setdefault('-s', 'c10') #Default coverage 10X
config.setdefault("-i", workdir+str('/Ecoli.fna')) #Default genome phi-x

if(config['-a']==1 or config['-a']=='1'): #Single End Sequencing
    config.setdefault('-x', workdir+'/0.fq')  #output file goes to current Directory
if(config['-a'][0]==2 or config['-a'][0]=='2'): #Paired End Sequencing Value is (2,insertsize)
     config.setdefault('-y', workdir+'/1.fq')
     config.setdefault('-z', workdir+'/2.fq')

config.setdefault('-g', 'l') #Genome Type
config.setdefault('-v', 'True') #Verbose mode: yes 


#Finding Errors
for option, value in config.items(): 

    if(option=='-a' and value[0] not in ('1', '2')): #Eg. -a 'd'makes no sense!
        error('a')
    if(option=='-a'and value[0]=='2'): #Expected Format for paired end sequencing is -a 2,200 where 200 is insert size
        try:
            testvar1=value[1]   #Testvar1 are variables for temporary use
            if(testvar1!=','):
                error('a')
            testvar2=int(value[2:])
        except ValueError:
            error('a')
        except IndexError:
            error('a')
    if(option=='-p' and value not in ('1', '2', '3', 1, 2, 3)): #Eg. -p 5 makes no sense!
        error('p') 
    if(option=='-g'  and   value not in ('l', 'c')): #Genome_Type : Circular/Linear
        error('g')
    if(option=='-v' and value not in ('True','False')):
        error('v')
    if(option=='-r'): #Readlength
        if(config['-p']!='3'):
            try:
                value=int(value)
            except ValueError:
                error('r')
        if(config['-p']=='3'): #Expected format -p mean, std
            try:
                m=value.split(', ')
                temp1=float(m[0])
                temp2=float(m[1])
            except AttributeError:
                error('r')
            except ValueError:
                error('r')
            except IndexError:
                error('r')
    if(option=='-s'):   #Stats  coverage and read numbers
        if(value[0] not in ('c', 'n')): #Expected format -s c50 or -s n5000000
            error('s')
        try:
            value=int(value[1:]) 
        except ValueError:
            error('s')





#Assigning values to variables

#Check for verbose first
ch=config['-v'] #Find user's choice for verbose
ch=bool_parse(ch) #Convert 'True' into True and 'False' into False
if ch==True:
    verbose=True
    print ('Verbose Mode On')
else:
    verbose=False


for option, value in config.items(): #A display of configuration
    if(option=='-a'):
        type1=value #SingleEnd/Paired
        if(type1[0]=='1'):
            sequencing_type=1
            if verbose:
                print('Sequencing Type selected : Single End Sequencing')
        if(type1[0]=='2'):
            insertsize=int(type1[2:])
            sequencing_type=2
            if verbose:
                print 'Sequencing Type Selected : Paired End Sequencing','\n','Insert Size :', insertsize

    if(option=='-p'):
        platform=int(value) #Illumina/SOlid/Pacb..
        if(platform==1):
            if verbose:
                print('Platform Selected: Illumina')
        if(platform==2):
            if verbose:
                print('Platform Selected: SOLiD')
        if(platform==3):
            if verbose:
                print('Platform Selected: PacBio')

    if(option=='-s'):
            if (value[0]=='c'):#First letter indicates if the user is entering coverage or read numbers required
                coverage=int(value[1:]) #Eg c30 means coverage of 30
                if verbose:
                    print 'Coverage Specified:', coverage
            if(value[0]=='n'):  #Total Numbers of Read to be generated
                readno=int(value[1:])
                if verbose:
                    print 'Read Numbers Specifed', readno

    if(option=='-r'):
        try:
            readlen=int(value)#ReadLength  except pacb
            if verbose:
                print 'ReadLength:', readlen
        except ValueError:    #For pacb
            mean_stddev=value.split(', ')
            mean=float(mean_stddev[0])
            stddev=float(mean_stddev[1])
            x= math.exp(mean+(stddev**2/2.0))  #This is mean readlength. It will be used for calculating coverage or No. of Reads
            readlen=math.exp(x)
            if verbose:
                print 'Readlength:', readlen

    if(option=='-i'):       #Input file path contains all chromosomes of a genome
        ipath=value
        try:
            inf=open(ipath, 'r') # inf: Input File Handle
            if verbose:
                print 'Input File Path:', ipath
        except IOError:     #If theres issue with input file
            print'Input File cannot be opened'
            sys.exit()

    if(option=='-x'):       #Output File path1 for Single End Sequencing
        opath1=value
        try:
            outf0=open(opath1, 'w')
            if verbose:
                print 'Output File Path:', opath1
        except:
            print 'Output File cannot be opened!'
    if(option=='-y'):
        opath2=value      #Output File path1 for paired end sequencing
        try:
            outf1=open(opath2, 'w')
            if verbose:
                print "Output File Path 1 3'-> 5'" , opath2
        except :
            print 'Output File cannot be opened'
    if(option=='-z'):   #OUtput file path2 for paired end sequencing
        opath3=value
        try:
            outf2=open(opath3, 'w')
            if verbose:
                print "Output File Path 2 5'--> 3'", opath3
        except:
            print 'Output File cannot be opened'
    if(option=='-g'):        #Genome Type - linear/circular
        genome_type=value
        if(genome_type=='l'):
            if verbose:
                print 'Genome Type', 'Linear'
        if(genome_type=='c'):
            if verbose:
                print 'Genome Type', 'Circular'

#Output thrown away
#outf1=open('/dev/null', 'w')
#outf0=open('/dev/null', 'w')
#outf2=open('/dev/null', 'w')

chromosome_data=fasta(inf) #Convert FASTA file into { idno:chromosome seq} dictionary format
if verbose:
    print('FASTA file Successfully Loaded!')

for idno, seq in chromosome_data.items():  #For every chromosome
    seq=seq.replace('\n', '')
    print len(seq)
    if verbose:
        print('Current Chromosome: ', idno)
    try: #If readNumbers are given generate coverage
        coverage=int(round(readno*readlen/float(len(seq))))
        if verbose:
            print('Coverage', coverage)
    except NameError: #Readno NOT given. Coverage is given
        readno=int(round((coverage*len(seq)/float(readlen))))
        if verbose:
            print('Total Number of reads', readno)

    if sequencing_type==1: #Single End Sequencing

        if verbose:
            counter=1 #This variable will be used to track how many reads have been generated
            if(10<=readno<100):
                cutoff=10    #After every 10 reads is written  a message will be displayed
            if(100<=readno<1000):
                cutoff=100  
            if(1000<=readno<10000):
                cutoff=1000
            if(10000<=readno<100000):
                cutoff=10000
            if(readno>=100000):
                cutoff=100000

        for i in range(readno):
            if(platform==1 or platform==2):
                read, loc = readg(seq, readlen, genome_type,sequencing_type) #Generates a random read from the chromosome sequence

            if(platform==1):
                read, qual, cigar_format=illumina_errorprofile(read, platform)

            if(platform==2): #SOLID
                cs=colorspace(read)
                read, qual, cigar_format=illumina_errorprofile(cs, platform) #cs sent for substitution errors
                read='A'+read  #Adaptor Base Added to the beginning

            if(platform==3): #Pacbio
                readlen=pacb_readlen(mean, stddev)
                while (readlen>50000):
                    readlen=pacb_readlen(mean, stddev)
                read, loc = readg(seq, int(round(readlen)), genome_type, sequencing_type)
                read, qual, cigar_format=pacb(read)  #Call the pacb function

            idn='@'+'Readno.'+str(i)+'Chrm'+str(idno)+'Loc'+str(loc)+str('_')+cigar_format+'\n'
            #First Line of FASTQ - i represents read no. and idno represent  the chromosome ID
            writer(outf0, idn, read, qual) #Writes everything on file

            if verbose:
                if(i==counter*cutoff):
                    print(i, ' reads have been generated and written successfully in FASTQ Format to the specified Directory!')
                    counter=counter+1


    if(sequencing_type==2): #Paired End Sequencing
        pairs=int(round(readno/2.0))
        if verbose:
            counter=1 #This variable will be used to track how many pairs have been generated
            if(10<=pairs<100):
               cutoff=10    #After every 10 reads is written  a message will be displayed
            if(100<=pairs<1000):
                cutoff=100
            if(1000<=pairs<10000):
                cutoff=1000
            if(10000<=pairs<100000):
                cutoff=10000
            if(pairs>=100000):
                cutoff=100000
        
        for i in range(pairs):

            if(platform==3):
                readlen=pacb_readlen(mean, stddev) #Generate a random readlength from normal distribution 
                while(readlen>50000):
                    readlen=pacb_readlen(mean, stddev)
                read3_5, read5_3, loc=readg(seq, readlen, genome_type, sequencing_type, insertsize)
                read3_5, qual35, cigar_format35=pacb(read3_5)  #Call the pacb function
                read5_3, qual53, cigar_format53=pacb(read5_3) 

            if( platform in (1,2)):
                read3_5, read5_3, loc=readg(seq,readlen, genome_type, sequencing_type, insertsize)
 
            if(platform==1):
                read3_5, qual35, cigar_format35=illumina_errorprofile(read3_5, platform)
                read5_3, qual53, cigar_format53=illumina_errorprofile(read5_3, platform)

            if(platform==2):
                cs35=colorspace(read3_5)
                cs53=colorspace(read5_3)
                read3_5, qual35, cigar_format35=illumina_errorprofile(cs35, platform)
                read5_3, qual53, cigar_format53=illumina_errorprofile(cs53, platform)

                read3_5='A'+read3_5   # 'A' is the adaptor base
                read5_3='A'+read5_3

            idn35='@'+'1.fq'+'Pair No.'+ str(i)+'chrm'+str(idno)+'Loc'+str(loc)+ cigar_format35+'\n'
            idn53='@'+'2.fq'+'Pair No.'+str(i)+'Chrm'+str(idno)+'Loc'+str(loc)+ cigar_format53+'\n'
            writer(outf1, idn35, read3_5, qual35) #outf1 : output file handle
            writer(outf2, idn53, read5_3, qual53)

            if verbose:
                if(i==counter*cutoff):
                    print i, 'pairs have been generated and written successfully in FASTQ Format to the specified Directory!'
                    counter=counter+1




if verbose:
    print('Total Reads Generated', readno)
    print('Operation Successful....')
inf.close()
if(sequencing_type==1):
    outf0.close()
if(sequencing_type==2):
    outf1.close()
    outf2.close()

