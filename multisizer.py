import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import re

#defining constants
countspervolt = 1/(4*298.02e-9)

def diameterToVolume(D):
    '''converts diameter to volume'''
    return (4.0/3.0) * np.pi *(D/2.0) **3.0
def volumeToDiameter(V):
    '''converts volume to diameter'''
    return np.cbrt(V * 0.75 / np.pi) * 2;


class coulterExperiment:
    def __init__(self,filename):
        self.filename = filename
        
        #open file, get data as a string
        with open(filename) as file:
            datastring = file.read()
            self.datastring = datastring
        
        linesList = datastring.split('\n')#break up data into an array with each line being an element of an array
        self.linesList = linesList
        
        
        
        #find sections 
        dataSections = re.findall('\[(.*)\]\n',datastring)
        self.dataSections = dataSections
        sectionIndices = [linesList.index('[' + section + ']') for section in dataSections] + [len(linesList)]#adding the last index so we can get the very end
        
        #construct dictionary with section data
        dataDict = {};
        
        #split datastring into the above dataSections, and put into a dictionary
        for i,section in enumerate(dataSections):
            dataDict[section] = linesList[sectionIndices[i]+1:sectionIndices[i+1]]
         
        #convert data that contains just numbers into floats or (if possible) ints
        for key,val in dataDict.iteritems():
            try:
                dataDict[key] = map(float,val)
                dataDict[key] = map(int,val)    
            except:
                pass

        #format the hex data into usable numbers    
        hex5Lists = [map(lambda x: int(x,16),subarray.split(',')) for subarray in dataDict['#Pulses5hex']]#converting from hexadecimal
        pulsesDataFrame = pd.DataFrame(hex5Lists,columns=('MaxHeight','MidHeight','Width','Area','Gain')) #convert from list of lists to useable dataframe with labeled columns
        self.pulsesDataFrame = pulsesDataFrame
        
        #convert parameter data into usable dictionaries
        for section in dataSections:
            sectionDict = {}
            for line in dataDict[section]: #going through all the data in this section, seeing if it fits the 'field=value' format

                try:
                    field,value = re.search('(.*)=(.*)',line).group(1,2) #looking for data of the form 'field=value'
                    sectionDict[field] = value #putting that data in the the section dictionary 
                    sectionDict[field] = float(value) #converting to float, if possible
                    sectionDict[field] = int(value) #converting to int, if possible
                except:
                    pass
                
            #replace the data in dataDict with the sectionDict, if it's of the proper format
            if len(sectionDict.keys()) > 0:
                dataDict[section] = sectionDict
            else:
                pass
        
        #put particular data of interest into object attributes for ease of reference
        self.instrumentData = dataDict['instrument']
        self.gain = self.instrumentData['Gain']
        self.current = self.instrumentData['Current']/1000.0 #converting to mA
        self.Kd = self.instrumentData['Kd']
        self.maxHeightCorr = self.instrumentData['MaxHtCorr']
        
        #calculating diameter from pulse data and other things
        self.pulsesDataFrame['height'] = self.pulsesDataFrame['MaxHeight'] + self.maxHeightCorr #correcting height
        self.pulsesDataFrame['diameter'] = self.Kd * np.cbrt((self.pulsesDataFrame.height) / (countspervolt * 25.0 * self.current * self.gain))
        self.pulsesDataFrame['volume'] = (4.0/3.0) * np.pi * (self.pulsesDataFrame['diameter']/2.0)**3.0
        
        
        self.dataDict = dataDict
        
        #calculating with default
        self.countCells()
    
    
    def countCells(self,boundType='diameter',minDiameter = 9.0, maxDiameter = 60.0, minVolume = diameterToVolume(9.0), maxVolume = diameterToVolume(60.0),dilution=100.0):
        '''sets the upper and lower size limit for what is considered a cell, and calculates the cell concentration in k/mL in that sample (subject to a dilution factor based on how the sample is prepared -- default is 100, since we usually put 100uL of cells into 10mL of solution for measurement.
        
        
        INPUTS:
        -------
        sortType : Determines whether to sort based on diameter or volume, and thereby whether to use minDiameter or minVolume for ounds
        minDiameter: in microns (um)
        maxDiameter: in microns (um)
        minVolume: in cubic microns (fL)
        max Volume: in cubic microns (fL). 
        
        CALCULATES:
        -----------
        count: how many cells, in k/mL
        meanDiameter: in microns (um)
        meanVolume: in cubic microns (fL)
        medianDiameter: in microns (um)
        medianVolume in cubic microns (fL)
        
        
        RETURNS:
        ----------
        
        
        '''
        if boundType.lower() in ['diameter','d','diam']:
            self.lowerBound,self.upperBound = minDiameter,maxDiameter
            self.boundType = 'diameter'
        elif boundType.lower() in ['volume','v','vol']:
            self.lowerBound,self.upperBound = minVolume,maxVolume
            self.boundType = 'volume'
        else:
            return None
        
        self.cellData = self.pulsesDataFrame[(self.pulsesDataFrame[self.boundType] >= self.lowerBound) & (self.pulsesDataFrame[self.boundType] <= self.upperBound)] #subsetting that data that represents just cells
        
        #calculating various data
        self.count = len(self.cellData) * dilution / 1000.0 #dividing by 1,000 to put it in k/mL
        self.meanDiameter = np.mean(self.cellData.diameter)
        self.meanVolume = np.mean(self.cellData.volume)
        self.medianDiameter = np.median(self.cellData.diameter)
        self.medianVolume = np.median(self.cellData.volume)
    

    
    def __str__(self):
        return '< Coulter counter experiment from \'{0}\' >'.format(self.filename)
    def __repr__(self):
        return self.__str__()

