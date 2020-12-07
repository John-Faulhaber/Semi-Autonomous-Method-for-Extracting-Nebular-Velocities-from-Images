# Original use of this code: https://scholar.colorado.edu/concern/undergraduate_honors_theses/fx719n336


########################
### Preliminary code ###
########################


import numpy as np
import os


#################################################################################################################################################################################
#################################################################################################################################################################################
#################################################################################################################################################################################


#############################
### defining the function ###
#############################


def generateregions(oldclumps,newclumps):   # In the instance of the object of interest being a nebula, "oldclumps" and "newclumps" may represent noticible structure movement in the gas. 
    """ 
     Generates:numbered line region file
               velocity labeled line region file
               numbered vector region file
               400 year scaled numbered vector region file  
     
     Parameters
     ----------
         oldclumps: String with file name of older image region file with clumps specified. No extension.
         newclumps: String with file name of newer image region file with clumps specified. No extension.
         
    
    
    returns
    -------
        N/A
        
        (See description)
    
    
    
    example
    -------
        Older clumps region file: 1995P5.reg
        Newer clumps region file: 1995P5.reg
        
        --> generateregions('1995P5','2011P5')
    
    
    
    """
    
    import numpy as np   # redundancy
    import os



    my_file = oldclumps+'.reg'  # the function input only requires the file name (no extension) for convenience. Here, the file extention is added for use.
    base = os.path.splitext(my_file)[0]
    os.rename(my_file, base + '.txt')
    
    
#################################################################################################################################################################################
#################################################################################################################################################################################
#################################################################################################################################################################################

    
# the following code combs through the two region files and extracts the RA and DEC coordinates of each selected region [requires region file structure from SAOImage ds9 ver. 4.1]

    
    regionfile = []                             
    with open (oldclumps+'.txt', 'rt') as myfile: 
        for line in myfile:                
            regionfile.append(line)           
    
    
    regioncoordinates=regionfile[3:]
  
    
    
    individualregioncoords=[]
    for i in regioncoordinates:
        individualregioncoords.append(i[i.find("(")+1:i.find(")")])
    
    
    ras1=[]
    decs1=[]
    for i in individualregioncoords: 
        ras1.append(i.split(',')[0])
        decs1.append(i.split(',')[1])
      
    
    deliminatedras=[]
    deliminateddecs=[]
    for i in ras1:
        i=i.split(':')
        deliminatedras.append(i)
        
    
    for i in decs1:
        i=i.split(':')
        deliminateddecs.append(i)
      
    
    NinetyFiveRAs=deliminatedras   # The RA coordinates of the older data [original use of this code used data from 1995 and 2011].
    NiteyFiveDecs=deliminateddecs   # The DEC coordinates of the newer data.
    
    
    my_file = oldclumps+'.txt'
    base = os.path.splitext(my_file)[0]
    os.rename(my_file, base + '.reg')

    
    my_file = newclumps+'.reg'
    base = os.path.splitext(my_file)[0]
    os.rename(my_file, base + '.txt')

    
    regionfile = []                            
    with open (newclumps+'.txt', 'rt') as myfile:
        for line in myfile:         
            regionfile.append(line)                                   
    
    
    regioncoordinates=regionfile[3:]
    
    
    individualregioncoords=[]
    for i in regioncoordinates:
        individualregioncoords.append(i[i.find("(")+1:i.find(")")])
    
    
    ras2=[]
    decs2=[]
    for i in individualregioncoords: 
        ras2.append(i.split(',')[0])
        decs2.append(i.split(',')[1])
   
    
    deliminatedras=[]
    deliminateddecs=[]
    for i in ras2:
        i=i.split(':')
        deliminatedras.append(i)
        
    
    for i in decs2:
        i=i.split(':')
        deliminateddecs.append(i)
   
    
    ElevenRAs=deliminatedras   # The RA coordinates of the newer data [original use of this code used data from 1995 and 2011].   
    ElevenDecs=deliminateddecs   # The DEC coordinates of the newer data.
    
    
    my_file = newclumps+'.txt'
    base = os.path.splitext(my_file)[0]
    os.rename(my_file, base + '.reg')
    

#################################################################################################################################################################################
#################################################################################################################################################################################
#################################################################################################################################################################################

    
# The following code converts the extracted RA and DEC coordinates to degrees.
 
    
    NinetyFiveRADegrees=[]
    for i in NinetyFiveRAs:
        degreevalue=float(i[0])+(float(i[1])/60)+(float(i[2])/3600)
        NinetyFiveRADegrees.append(degreevalue)
    
    
    NinetyFiveDECDegrees=[]
    for i in NiteyFiveDecs:
        degreevalue=float(i[0])+(float(i[1])/60)+(float(i[2])/3600)
        NinetyFiveDECDegrees.append(degreevalue)

        
    ElevenRADegrees=[]
    for i in ElevenRAs:
        degreevalue=float(i[0])+(float(i[1])/60)+(float(i[2])/3600)
        ElevenRADegrees.append(degreevalue)

    
    ElevenDECDegrees=[]
    for i in ElevenDecs:
        degreevalue=float(i[0])+(float(i[1])/60)+(float(i[2])/3600)
        ElevenDECDegrees.append(degreevalue)
        
    
#################################################################################################################################################################################
#################################################################################################################################################################################
#################################################################################################################################################################################


# Calculate the distance between coordinate pairs and velocity of the movement.


    lengths=[]
    for i in range(len(NinetyFiveRADegrees)):
        Coslengtharcsec=206265*(2*np.arcsin(((np.sin((((float(ElevenDECDegrees[i])*np.pi)/180)-((float(NinetyFiveDECDegrees[i])*np.pi)/180))/2)**2+(np.cos(((float(NinetyFiveDECDegrees[i])*np.pi)/180))*np.cos(((float(ElevenDECDegrees[i])*np.pi)/180))*np.sin(((float(ElevenRADegrees[i])*np.pi)/180)-((float(NinetyFiveRADegrees[i])*np.pi)/180))**2))**0.5)))
        lengths.append(Coslengtharcsec)


##############        
#### NOTE ####
##############


# 1091 in following calculation is distance in pc to original structure this code was made for. Alter as necessary. Units must be pc.
# 15 in the following calculation accomodates for the ~15 year interval between the the original data this code was made for. Alter as necessary. Units must be years.

    tangentialvelovities=[]
    for i in range(len(lengths)):
        tangentialvelovities.append(round(((lengths[i]/15)*1091*4.74),1))
        
        
#################################################################################################################################################################################
#################################################################################################################################################################################
#################################################################################################################################################################################


# Generate the four additional vector files
# Output files will be in same file location as this code
        

    while True:
        answer=input("\nWould you like a numbered line region file? y/n.\n\n                                                 Answer: ")
        if answer=='n':  
            print( '\n  noted')
            break
        else:
            if answer=='y':
                numberedlineregionsfilename=input("\nName for numbered line region file?.\n\n                                                 Answer: ")
                f= open(numberedlineregionsfilename+'.txt',"w+")
                f.write('# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n')
                f.close()
                for i in range(len(lengths)):    
                    with open(numberedlineregionsfilename+'.txt', "a") as myfile:
                        myfile.write("line({a},{b},{c},{d}) # line=0 0 text={e}\n".format(a=ras1[i],b=decs1[i],c=ras2[i],d=decs2[i],e={i+1}))
                print('(☞ﾟヮﾟ)☞ you got it, mate!')
                fd=open(numberedlineregionsfilename+'.txt',"r")
                d=fd.read()
                fd.close()
                m=d.split("\n")
                s="\n".join(m[:-1])
                fd=open(numberedlineregionsfilename+'.txt',"w+")
                for i in range(len(s)):
                    fd.write(s[i])
                fd.close()
                my_file = numberedlineregionsfilename+'.txt'
                base = os.path.splitext(my_file)[0]
                os.rename(my_file, base + '.reg')       
                break
            else: 
                print( '\n\n ~ERROR~\n ~ERROR~\n ~ERROR~\n ~ERROR~\n\n\n    Unrecognized input ¯\_(ツ)_/¯ (try again)\n\n\n\n\n\n')
    
    
    while True:
        answer=input("\nWould you like a velocity labeled line region file? y/n.\n\n                                                 Answer: ")
        if answer=='n': 
            print( '\n  noted')
            break
        else:
            if answer=='y':
                velocitylabeledlineregionsfilename=input("\nName for velocity labeled line region file?.\n\n                                                 Answer: ")
                f= open(velocitylabeledlineregionsfilename+'.txt',"w+")
                f.write('# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n')
                f.close()
                for i in range(len(lengths)):
                    with open(velocitylabeledlineregionsfilename+'.txt', "a") as myfile:
                        myfile.write("line({a},{b},{c},{d}) # line=0 0 text={e}\n".format(a=ras1[i],b=decs1[i],c=ras2[i],d=decs2[i],e={tangentialvelovities[i]}))
                print('o(^▽^)o done!')
                fd=open(velocitylabeledlineregionsfilename+'.txt',"r")
                d=fd.read()
                fd.close()
                m=d.split("\n")
                s="\n".join(m[:-1])
                fd=open(velocitylabeledlineregionsfilename+'.txt',"w+")
                for i in range(len(s)):
                    fd.write(s[i])
                fd.close()
                my_file = velocitylabeledlineregionsfilename+'.txt'
                base = os.path.splitext(my_file)[0]
                os.rename(my_file, base + '.reg')                
                break
            else: 
                print( '\n\n ~ERROR~\n ~ERROR~\n ~ERROR~\n ~ERROR~\n\n\n    Unrecognized input ¯\_(ツ)_/¯ (try again)\n\n\n\n\n\n')
      
    
    while True:
            answer=input("\nWould you like a numbered vector region file? y/n.\n\n                                                 Answer: ")
            if answer=='n':
                print( '\n  noted')
                break
            else:
                if answer=='y':
                    numberedvectorregionsfilename=input("\nName for numbered vector region file?.\n\n                                                 Answer: ")    
                    f= open(numberedvectorregionsfilename+'.txt',"w+")
                    f.write('# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n')
                    f.close()
                    for i in range(len(lengths)):
                        with open(numberedvectorregionsfilename+'.txt', "a") as myfile:
                            myfile.write('# vector({a},{b},{c},{d}) vector=1 text={e}\n'.format(a=ras1[i],b=decs1[i],c=(((lengths[i]/206265)*180)/np.pi),d='',e={i+1}))
                    print('(⌐■_■) no problemo')
                    fd=open(numberedvectorregionsfilename+'.txt',"r")
                    d=fd.read()
                    fd.close()
                    m=d.split("\n")
                    s="\n".join(m[:-1])
                    fd=open(numberedvectorregionsfilename+'.txt',"w+")
                    for i in range(len(s)):
                        fd.write(s[i])
                    fd.close()
                    my_file = numberedvectorregionsfilename+'.txt'
                    base = os.path.splitext(my_file)[0]
                    os.rename(my_file, base + '.reg')    
                    break
                else: 
                    print( '\n\n ~ERROR~\n ~ERROR~\n ~ERROR~\n ~ERROR~\n\n\n    Unrecognized input ¯\_(ツ)_/¯ (try again)\n\n\n\n\n\n')
      
    
    while True:
            answer=input("\nWould you like a 400 year scaled numbered vector region file? y/n.\n\n                                                 Answer: ")
            if answer=='n':
                print( '\n  noted')
                break
            else:
                if answer=='y':
                    scalednumberedlineregionsfilename=input("\nName for 400 year scaled numbered vector region file?.\n\n                                                 Answer: ")
                    f= open(scalednumberedlineregionsfilename+'.txt',"w+")
                    f.write('# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n')
                    f.close()
                    for i in range(len(lengths)):
                        with open(scalednumberedlineregionsfilename+'.txt', "a") as myfile:
                            
##############
#### NOTE ####
##############

# The time difference between when the data was acquired directly corresponds to the distances between each corresponnding set of regions. In this original use case of this code, just under 16 years. As a result, to extend the vectors' lengths by 400 years movement in this particular instance, they needed to be multiplied by 26.6666... NOTE: This method of estimation assumes constant velocity. Alter as needed.


                            myfile.write('# vector({a},{b},{c},{d}) vector=1 text={e}\n'.format(a=ras1[i],b=decs1[i],c=(((lengths[i]/206265)*180)/np.pi)*26.66666666666666666666666666666666,d='',e={i+1}))   
                    print('ez pz  d=====(￣▽￣*)b')
                    fd=open(scalednumberedlineregionsfilename+'.txt',"r")
                    d=fd.read()
                    fd.close()
                    m=d.split("\n")
                    s="\n".join(m[:-1])
                    fd=open(scalednumberedlineregionsfilename+'.txt',"w+")
                    for i in range(len(s)):
                        fd.write(s[i])
                    fd.close()
                    my_file = scalednumberedlineregionsfilename+'.txt'
                    base = os.path.splitext(my_file)[0]
                    os.rename(my_file, base + '.reg')
                    break
                else:
                    print( '\n\n ~ERROR~\n ~ERROR~\n ~ERROR~\n ~ERROR~\n\n\n    Unrecognized input ¯\_(ツ)_/¯ (try again)\n\n\n\n\n\n')   

