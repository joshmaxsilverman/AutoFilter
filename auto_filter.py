import os
import string
import sets
import csv
import re
import numpy
path = '/Users/joshsilverman/Dropbox/Research/AutoFilter/gogat40c'

os.chdir(path)

######
#body#
######

#Get a list of filename stems

filelist = os.listdir(path)[2:]
stems = []

# for file in filelist:
#     for type in [".png",".contour", ".newcon", ".txt" , ".dat", ".fit"]:
#         file = string.replace(file, type, "") 
#     stems.append(file)
# stems = list(sets.Set(stems))

for file in filelist:
    stems.append(re.split(".[a-z]", file)[0])
stems = list(sets.Set(stems))

#Reading in the fit and dat files as numpy arrays

###########
#functions#
###########

def setup(stem):
    fit_file = open(stem + ".fit")
    global fit    
    fit = numpy.array([map(float, item) for item in csv.reader(fit_file)])
    global mzfit
    mzfit = fit.transpose()[0, 0:]
    
    dat_file = open(stem + ".dat")
    global data
    data = numpy.array([map(float, item) for item in csv.reader(dat_file)])
    global mzdat    
    mzdat = data.transpose()[0, 0:]

def myrange(input):
    return max(list(input)) - min (list(input))
    
def setuprange():
    global values
    values = fit.transpose()[1,0:]
    global valuesLength
    valuesLength = len(values)
    global punit
    punit = int(valuesLength/8)
    
    global partition #zip(*iter( )]*j) equivalent to Mathematica partition function
    partition = zip(*[iter(values)]*punit)[0:7]
    global ran    
    ran = map(myrange, partition)
    global rangeaccesshash
    rangeaccesshash = {}
    for i in range(len(ran)):
        rangeaccesshash[ran[i]] = i
    
    #Start expanding the "fitless" region from the center and endpoints
    #the "minimum range" is the region with the least spread in peak intensity
    global minrange    
    minrange = min(ran)
    global mymin
    mymin = rangeaccesshash[minrange]
    global baseline
    baseline = sum(partition[mymin])/len(partition[mymin])
    
    global leftborder, rightborder
    [leftborder, rightborder] = [mymin*punit + 1 - 1, (mymin + 1)*punit - 1]
    
    global peak1rightborder, peak2leftborder, peak2rightborder, peak1leftborder
    peak1rightborder = creepUp(leftborder , rightborder , 1000 , values, -1)
    peak2leftborder  = creepUp(leftborder , rightborder , 1000 , values, +1)
    peak2rightborder = creepUp(valuesLength - 1 , valuesLength - 1, 1000 , values, -1)
    peak1leftborder  = creepUp(0 , 0 , 1000 , values, +1)
    
    #minor disagreement between peak2rightborder in Mathematica and python versions, unresolved
    
    global fitlessrange
    fitlessrange = [peak1leftborder, peak1rightborder, peak2leftborder, peak2rightborder]
    
    #protects against the use of single isotope peaks in defining the datascale
    global FLR1, FLR2, datascale
    [FLR1, FLR2] = [[peak1leftborder, peak1rightborder], [peak2leftborder, peak2rightborder]]
    if cmp(FLR1[1]-FLR1[0], 0) < 0 or cmp(FLR2[1]-FLR2[0], 0) < 0:
        pass
    else:
        tempset= [ max( fit[FLR1[0]:FLR1[1]].transpose()[1]), max( fit[FLR2[0]:FLR2[1]].transpose()[1])]
        datascale = min(tempset)
    
def creepUp(LB, RB, EA, Val, direction):
    global TLB, TRB, TEA, temprange
    [TLB, TRB, TEA] = [LB, RB, EA]
    
    if direction < 0 :
        while TEA > 50 and TLB - TEA >= 0:
            temprange = Val[(TLB-TEA):TRB]
            if myrange(temprange) < 10 * minrange:
                TLB = TLB - TEA
            if myrange(temprange) > 10 * minrange:
                TEA = int(TEA/2)
        return TLB
    
    if direction > 0 :
        while TEA > 50 and TRB + TEA <= len(fit):
            temprange = Val[TLB:(TRB + TEA)]
            if myrange(temprange) < 10 * minrange:
                TRB = TRB + TEA
            if myrange(temprange) > 10 * minrange:
                TEA = int(TEA/2)
        return TRB

def xcoord(x):
    return fit[x][1]
    
#finds the peaks within the region fed to it [FLR1, FLR2]

def findpeaks(FLRlist):
    abovebelow = -1
    peaks = []
    global temppeaks
    temppeaks = []
    
    for i in range(FLRlist[0], FLRlist[1] + 1):
        pass      
        if cmp( fit[i][1] - 1.1*baseline,0) != abovebelow:
            temppeaks.append(i)
            abovebelow *= -1
        if len(temppeaks) == 2:
            tempfit = fit.transpose()[1]
            peakscale = max( tempfit[temppeaks[0]:temppeaks[1] + 1] )
            if peakscale > datascale / 25:
                peaks.append(temppeaks)
            temppeaks = []
        
    return peaks 

def nearest(array, value):
    index = (numpy.abs(array - value)).argmin()
    return array[index]

#makes pairs between the data and fit points within the region fed to it
def pairmaker(peaks):
    global dataendpointsmz
    peakendpoints = peaks
    
    def nearestdata (x):
        return nearest(mzdat,x)
    def nearestfit (x):
        return nearest(mzfit, x)
    
    expandedpeakendpoints = map( lambda x: nearestfit(nearestdata( mzfit[peakendpoints[x]] ) ), [0, 1])    
    expandedpeakendpoints = map( lambda y: [i for i, x in enumerate(mzfit) if x == expandedpeakendpoints[y]][0], [0,1])
    global tempregionfit    
    tempregionfit = ( lambda x: fit[x[0]:x[1]] )(expandedpeakendpoints)

    def nearesttempfit(x):
        return nearest(tempregionfit.transpose()[0], x)

    tempends = map( lambda y: tempregionfit[y][0], [1,-1] )
    dataendpointsmz = map(nearestdata, tempends)
    
    dataendpointsindex = map( lambda y: [i for i, x in enumerate(mzdat) if x == y][0], dataendpointsmz)
    tempregiondata = ( lambda y: data[ y[0]:y[1] + 1] )( dataendpointsindex ).transpose()[0:2].transpose()
    
    global datapairs
    datapairs = []
    
    for i in range(0, len( tempregiondata ) ):
        [dataposition, datavalue] = tempregiondata[i]
        fitposition = (lambda y: [j for j, x in enumerate(tempregionfit.transpose()[0] ) if x == y])( nearesttempfit(dataposition) )[0]       
        fitvalue = tempregionfit[fitposition][1]
        datapairs.append([datavalue, fitvalue])
    return numpy.array(datapairs)

def powersdifference (point, n):
    templist = map( lambda i: ( point[i][0] - point[i][1] ) ** n, range(0, len(point) )  )
    return sum( templist )
    
def absdifference (point):
    templist = map( lambda i: abs( point[i][0] - point[i][1] ) , range(0, len(point) ) )
    return sum( templist )

def score ( peakchoice, n ):
    datapairs = pairmaker( peakchoice )
    absdiff = absdifference( datapairs ) #powersdifference( datapairs, n )
    absarea = sum (datapairs.transpose()[1] )
    return absdiff / absarea

def scoren (peakschoice):
    return score(peakschoice, 2)

def doscoring ():
    global scorecollection
    peaks = map( findpeaks, [FLR1, FLR2] )
    scorecollection = map( lambda y: map( scoren, peaks[y]), [0,1] )
    tempscore = map( max, scorecollection )
    return max ( tempscore )    
    
import time
time1 = time.time()
valuedict = {}

for item in stems:
    setup(item)
    setuprange()
    while True:
        try:
            tempscore = doscoring()
            print "Score for %s is %0.4s" % (item, tempscore)
            valuedict[item] = tempscore
            break
        except NameError:
            break
        except ValueError:
            break
        except IndexError:
            break

print "process took " + str(time.time() - time1) + " seconds"

###################
#scored CSV output#
###################

csvpath = '/Users/joshsilverman/Dropbox/Research/AutoFilter/gogat_isocsv/gogat40_iso.csv'
scoredpath = '/Users/joshsilverman/Dropbox/Research/AutoFilter/gogat_isocsv/gogat40_iso_scored.csv'
datacsv = open(csvpath)
count = 0

scoredcsv = open(scoredpath, 'wr+')


for line in datacsv:
    if count > 0:
        tempkey = line.split(",")[-1].split("/")[1].split(".txt")[0]
        try:
            tempscore = str(valuedict[tempkey])
            scoredcsv.write( line.strip() + ", , " + tempscore + "\n" )
        except KeyError:
            pass
        count += 1
    else:
        count += 1
        scoredcsv.write( line.strip() + "," + "SCORE\n" )


scoredcsv.close()





