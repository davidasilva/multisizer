'''A collection of functions and regular expression patterns that are broadly useful for extracting experiment data from filenames, via the regExColum method of coulterExperiment objects.'''

pmDict = {
    '+':True,
    '-':False}

def plusMinusToBool(pm):
    if pm not in ['+','-']: return None
    else: return pmDict[pm]
    
drugPattern = '([+-])drug'
experimentTimePattern = 't=([0-9]*\.*[0-9])[a-zA-Z]*__' #gets the float part of the time in, e.g., 't=1.4d'
experimentTimeUnitPattern = 't=[0-9]*\.*[0-9]([a-zA-Z]*)__' #gets the time unit part of, e.g. 't=1.4d' (would return 'd')
cellTypePattern = '[0-9]*-[0-9]*__(.*)__[+-]'
dateTimeFormat = '%H:%M:%S  %d %b %Y' #format for recording the experiment time in the .#m4 files.