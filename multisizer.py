import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import re
import os

#defining constants
countspervolt = 1/(4*298.02e-9)

def diameterToVolume(D):
    '''converts diameter to volume'''
    return (4.0/3.0) * np.pi *(D/2.0) **3.0
def volumeToDiameter(V):
    '''converts volume to diameter'''
    return np.cbrt(V * 0.75 / np.pi) * 2;


class coulterExperiment(object):
    def __init__(self,filename):
        self.filename = filename
        self.filenameWithoutPath = os.path.split(self.filename)[1]
        self.fileTitle = self.filenameWithoutPath[0:-4]
        
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
        self.summaryData = pd.DataFrame(columns=('count','meanDiameter','meanVolume','medianDiameter','medianVolume'))
        self.countCells()
    
    
    def countCells(self, boundType='diameter',minDiameter = 9.0, maxDiameter = 25.0, minVolume = diameterToVolume(9.0), maxVolume = diameterToVolume(25.0),dilution=100.0):
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
        
        self.summaryData.loc[self.filename,['count','meanDiameter','meanVolume','medianDiameter','medianVolume']] = self.count, self.meanDiameter, self.meanVolume, self.medianDiameter, self.medianVolume
    
    def histogram(self,dataType='diameter',bins=50,**kwargs):
        '''Plots a histogram of the dataType column in pulsesDataFrame. kwargs get fed into plot; bins goes into np.histogram'''
        #constructing and plotting the histogram
        freq,bins = np.histogram(self.pulsesDataFrame[dataType],bins=bins)
        plt.fill_between(bins[0:-1],freq,alpha=0.4,**kwargs)
        
        #drawing vertical lines where the cell bounds are
        plt.axvline(x=self.lowerBound,linestyle=':')
        plt.axvline(x=self.upperBound,linestyle=':')
        
        #setting the xlims appropriately
        plt.xlim(bins[0],max(bins[-1],self.upperBound*1.05)); #lower bound is the lowest size detected; upper bound is either the biggest sized cell detected or the upperBound we've chosen for the cellular size range -- whichever is bigger
        
    def regExColumn(self,pattern,columnName):
        '''Creates a new column in summaryData based on the results of a regular expression search of pattern, and names that column columnName.'''
        searchResult = re.search(pattern,self.filename)
        
        if searchResult is None:
            print 'No match found.'
            self.summaryData.loc[self.filename,columnName] = None;
        else:
            data = searchResult.group(1)
            self.summaryData.loc[self.filename,columnName] = data;
    
    def __str__(self):
        return '< Coulter counter experiment from \'{0}\' >'.format(self.filename)
    def __repr__(self):
        return self.__str__()


class batchExperiment:
    def __init__(self,source):
        '''Creates an object for a collection of coulter counter files . Source can be either a list of coulterExperiment objects, a folder path, or a list of filenames. Can also use the __add__ method to combine coulterExperiment objects.'''
        
        #figure out what the source type is
        if type(source) == str:#the source is a string for a folder path
            filenameList = [os.path.join(source,filename) for filename in os.listdir(source) if filename[-4:].lower() == '.#m4'] #getting all the .#m4 files in the folder
            self.experimentList = [coulterExperiment(source_object) for source_object in filenameList]
        else:#must be some sort of iterable
            assert hasattr(source,'__iter__'), 'Invalid source type. Source must be either a file path or a list (or similar).'
            if all([object_type == str for object_type in map(type,source)]):#if every item in the source array is a string:
                   self.experimentList = [coulterExperiment(filename) for filename in source]; #creating coulterExperiment from every file in list
                                          
            elif all([isinstance(source_object,coulterExperiment) for source_object in source]):#if every object is a coulterExperiment object
                   self.experimentList = list(source); #forcing into a list
                    
                    
        #gather summaryData dataframes from all the individual files, concatenate into one
        self.updateSummaryData()
        
        
    def updateSummaryData(self):
        self.summaryData = pd.concat([experiment.summaryData for experiment in self.experimentList]);
        
        
    def regExColumn(self,pattern,columnName):
        '''For each experiment in experimentList, creates a new column in summaryData based on the results of a regular expression search of pattern, names that column columnName, then updates the batch summaryData.'''
        
        for experiment in self:
            experiment.regExColumn(pattern,columnName)
            
        self.updateSummaryData()
        
    def histogramArray(self,subplotShape=None,**kwargs):
        '''Creates an array of subplots of histograms for each experiment.'''
        nExps = len(self)
        
        if subplotShape is None:
            nRows = np.floor(np.sqrt(nExps))
            nColumns = nExps / nRows + 1
            
        fig = plt.figure(figsize=(12,8))
        for i, experiment in enumerate(self):
            plt.subplot(nRows,nColumns,i+1)
            experiment.histogram(**kwargs)
            plt.title(experiment.fileTitle)
        
        plt.tight_layout()
        
    def __getitem__(self,key):
        return self.experimentList[key]
    
    def __len__(self):
        return len(self.experimentList)
            
        