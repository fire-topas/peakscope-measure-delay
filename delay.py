"""
This program uses the libary peakscope from Jonas Schäfers peakscope-project on GitHub (https://github.com/horazont/peakscope, 14th January 2023).
To run this you can find the needed file at https://github.com/horazont/peakscope/blob/main/peakscope.py as of the 14th January 2023.
Parts of the program that interact with the libary may contain portions of code from Jonas Schäfers project.
If so, they are marked with a "s1" in a comment.

All of the code used from the peakscope-projekt is under the MIT license:
Copyright 2016 Jonas Schäfer

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import peakscope

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import RangeSlider, Button, TextBox
import scipy
import tkinter

class Dataset():
    # loads one ore more binary file with data and writes it in numpy arrays
    def __init__(self, paths):
        peakscopeBin = peakscope.Bin(open(paths[0], "rb").read()) #s1

        self.xRange = 15*peakscopeBin.channels[0].timescale * 1e-9 #s1
        self.xData = np.linspace(-self.xRange/2, self.xRange/2, len(peakscopeBin.channels[0].data)) #s1
        self.yData = np.array([peakscopeBin.channels[0].data * peakscopeBin.channels[0].voltscale * 1e-3, 
                            peakscopeBin.channels[1].data * peakscopeBin.channels[0].voltscale * 1e-3]) #s1
        for path in paths[1:]:
            peakscopeBin = peakscope.Bin(open(path, "rb").read()) #s1

            self.yData += np.array([peakscopeBin.channels[0].data * peakscopeBin.channels[0].voltscale * 1e-3, 
                                peakscopeBin.channels[1].data * peakscopeBin.channels[0].voltscale * 1e-3]) #s1

    # converts seconds into a sample number
    def timeToSamples(self, time, offset=False):
        if(offset):
            time = time + self.xRange/2
        return int(time/self.xRange*len(self.xData))

    # converts a sample number into seconds
    def samplesToTime(self, samples):
        return samples/len(self.xData)*self.xRange

# a sinus function with parameters a, b, c, d
def sinus(x, a, b, c, d):
    return a * np.sin(b * x + c) + d

# derive an array
def derive(data):
    delta_y = (data - np.roll(data, 1))[1:]
    return delta_y

# find zeros of a function
def findX0s(data):
    lastValue = 0
    x0s = np.array([])
    for index, value in enumerate(data):
        if(value<0 and lastValue>0):
            x0s = np.append(x0s, index-0.5)
        elif(value>0 and lastValue<0):
            x0s = np.append(x0s, index-0.5)
        elif(value==0):
            x0s = np.append(x0s, index)               
        lastValue = value
    return x0s

# find proposal Parameters to fit the sinus to the data
def proposalParameters(data, dataset, extremumsRange, turning_pointsRange, proposalCAverageRange=10):
    propasalD = np.average(data) # the proposed d is the average of the data
    proposalA = (np.max(data)-np.min(data))/2 # proposed amplitude is calculated from the minumum and maximum values

    extremumDistances = (np.roll(np.array(extremumsRange), -1) - np.array(extremumsRange))[0:-1] # the distances between the extremums
    turning_pointsDistances = (np.roll(np.array(turning_pointsRange), -1) - np.array(turning_pointsRange))[0:-1] # the distances between the turning points

    # calculate a proposed B from the distances between extremums and turning points
    distances = np.append(extremumDistances, turning_pointsDistances)
    if(len(distances)!=0):
        distancesMin = np.min(distances); distancesMax = np.max(distances)
        distancesCenter = (distancesMax - distancesMin)/2
        T = np.average(distances[np.where(distances>distancesCenter)])
        proposalB = np.pi/dataset.samplesToTime(T)
    else:
        proposalB = 1

    # calculate a proposed C from all the other parameters and the beginning of the data
    y0 = np.average(data[0:proposalCAverageRange])
    x0 = proposalCAverageRange/2
    proposalC = (np.arcsin((y0-propasalD)/proposalA) - proposalB * x0)%2*np.pi

    return [proposalA, proposalB, proposalC, propasalD]

# filter data with a order 3 lowpass butterworth filter as shown at https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.lfilter.html as of the 14th January 2023
def filterData(data):
    b, a = scipy.signal.butter(3, 0.05)

    zi = scipy.signal.lfilter_zi(b, a)
    z, _ = scipy.signal.lfilter(b, a, data, zi=zi*data[0])

    z2, _ = scipy.signal.lfilter(b, a, z, zi=zi*z[0])

    return scipy.signal.filtfilt(b, a, data)

# fit a sinus through the data with scipy starting with propesed parameters
def fitSinus(xData, yData, sigma, dataset, extremumsRange, turning_pointsRange):
    p0 = proposalParameters(filterData(yData), dataset, extremumsRange, turning_pointsRange)
    return scipy.optimize.curve_fit(sinus, xData, yData, p0=p0, sigma=np.full((len(xData)), sigma), absolute_sigma=True), p0

class xRangeSlider():
    # init the class and create a placeholder slider
    def __init__(self, rect, title):
        self.rect = rect
        self.axes = plt.axes(rect)
        self.title = title
        self.slider = RangeSlider(self.axes, title, 0, 1)
    
    # create a new slider after data is loaded
    def onDataLoaded(self, maxRange, startRange, valStep, updateCallback):
        self.axes.remove()
        self.axes = plt.axes(self.rect)
        self.range = startRange
        self.slider = RangeSlider(self.axes, self.title, *maxRange, valinit=startRange, valstep=valStep)
        self.slider.on_changed(updateCallback)

class Plot():
    def __init__(self, dataset = None):
        # initialize the plots, widgets and tables
        self.fig, self.ax = plt.subplots(2, 3, gridspec_kw={'width_ratios': [3, 3, 1]}, constrained_layout = True)
        
        self.fig.tight_layout()
        self.fig.subplots_adjust(bottom=0.35)

        self.xRangeSliders = [xRangeSlider([0.1, 0.15, 0.8, 0.03], "Data 1 x-Range"), xRangeSlider([0.1, 0.1, 0.8, 0.03], "Data 2 x-Range")]

        openDatasetAxis = plt.axes([0.1, 0.25, 0.2, 0.05])
        openDatasetButton = Button(openDatasetAxis, "Open File")
        openDatasetButton.on_clicked(lambda value: self.getPath(value))

        toggleButtonAxis = plt.axes([0.4, 0.25, 0.2, 0.05])
        self.toggleButton = Button(toggleButtonAxis, "Toggle Measure Points")
        self.toggleButton.on_clicked(self.toggleImportantPoints)

        calculateButtonAxis = plt.axes([0.7, 0.25, 0.2, 0.05])
        calculateButton = Button(calculateButtonAxis, "Calculate")
        calculateButton.on_clicked(self.displayTimeDelay)

        errorAxis = plt.axes([0.1, 0.05, 0.2, 0.03])
        self.errorTextBox = TextBox(errorAxis, "error = ", "")

        setDefaultRangesAxis = plt.axes([0.4, 0.05, 0.2, 0.03])
        self.setDefaultRangesButton = Button(setDefaultRangesAxis, "Set Default Range")
        self.setDefaultRangesButton.on_clicked(self.setDefaultRanges)
        self.defaultRange = None
        self.defaultRangeSet = False

        timeDelayAxis = plt.axes([0.7, 0.05, 0.2, 0.03])
        self.timeDelayTextBox = TextBox(timeDelayAxis, "time delay = ", str(0))

        self.ranges = [[0, 0], [0, 0]]
        self.popts = [[0,0,0,0],[0,0,0,0]]

        self.parameterTables = []
        self.ax[0][2].set_axis_off() 
        self.parameterTables.append(self.ax[0][2].table(cellText = [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]], 
                                                                    rowLabels = ["a\n(Amplitude)", "b\n(Winkelgeschw.)", "c\n(x shift)", "d\n(y shift)"],  
                                                                    colLabels = ["proposed", "fitted", "delta"], loc='center'))
        self.parameterTables[0].auto_set_font_size(False)
        self.parameterTables[0].set_fontsize(5.5)
        self.parameterTables[0].scale(1.5, 1.5)
        

        self.ax[1][2].set_axis_off() 
        self.parameterTables.append(self.ax[1][2].table(cellText = [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]], 
                                                                    rowLabels = ["a\n(Amplitude)", "b\n(Winkelgeschw.)", "c\n(Winkelgeschw.)", "d\n(y shift)"],  
                                                                    colLabels = ["proposed", "fitted", "delta"], loc='center'))
        self.parameterTables[1].auto_set_font_size(False)
        self.parameterTables[1].set_fontsize(5.5)
        self.parameterTables[1].scale(1.5, 1.5)
        
        # load the data if possible
        if dataset != None:
            self.loadData(dataset)

        plt.show()

    # method called, when the "Set Default Range" button is clicked, that saves the current slider state
    def setDefaultRanges(self, value):
        self.defaultRange = np.copy(self.ranges)
        self.defaultRangeSet = True

    # loads the data into the plots
    def loadData(self, dataset):
        self.dataset = dataset
        self.ax[0][0].plot(dataset.xData, dataset.yData[0], alpha=0.5, color="blue"); self.ax[0][1].plot(dataset.xData, dataset.yData[0], alpha=0.5, color="blue")
        self.ax[0][0].plot(dataset.xData, filterData(dataset.yData[0]), color="blue"); self.ax[0][1].plot(dataset.xData, filterData(dataset.yData[0]), color="blue")
        self.ax[1][0].plot(dataset.xData, dataset.yData[1], alpha=0.5, color="blue"); self.ax[1][1].plot(dataset.xData, dataset.yData[1], alpha=0.5, color="blue")
        self.ax[1][0].plot(dataset.xData, filterData(dataset.yData[1]), color="blue"); self.ax[1][1].plot(dataset.xData,  filterData(dataset.yData[1]), color="blue")

        if(self.defaultRangeSet):
            self.ranges = np.copy(self.defaultRange)
        else:
            self.ranges = [[-dataset.xRange/2, dataset.xRange/2], [-dataset.xRange/2, dataset.xRange/2]]

        self.xRangeSliders[0].onDataLoaded((-dataset.xRange/2, dataset.xRange/2), self.ranges[0], dataset.xRange/(len(dataset.xData)-1), lambda value: self.xRangeSliderUpdate(value, 0))
        self.xRangeSliders[1].onDataLoaded((-dataset.xRange/2, dataset.xRange/2), self.ranges[1], dataset.xRange/(len(dataset.xData)-1), lambda value: self.xRangeSliderUpdate(value, 1))

        self.plotBorders = [[self.ax[0][0].axvline(x=-dataset.xRange/2, color='red', linestyle=':'), self.ax[0][0].axvline(x=dataset.xRange/2, color='red', linestyle=':')],
                            [self.ax[1][0].axvline(x=-dataset.xRange/2, color='red', linestyle=':'), self.ax[1][0].axvline(x=dataset.xRange/2, color='red', linestyle=':')]]
        
        self.toggleImportantPointsValue = False
        self.extremums = [findX0s(derive(filterData(dataset.yData[0]))), findX0s(derive(filterData(dataset.yData[1])))]
        self.extremumsPlot = [self.ax[0][1].plot(self.dataset.samplesToTime(self.extremums[0].astype(int)) - self.dataset.xRange/2, filterData(self.dataset.yData[0])[self.extremums[0].astype(int)], marker='.', linestyle='none', markersize=10, color='green'), 
                              self.ax[1][1].plot(self.dataset.samplesToTime(self.extremums[1].astype(int)) - self.dataset.xRange/2, filterData(self.dataset.yData[1])[self.extremums[1].astype(int)], marker='.', linestyle='none', markersize=10, color='green')]

        self.turning_points = [findX0s(derive(derive(filterData(dataset.yData[0])))), findX0s(derive(derive(filterData(dataset.yData[1]))))]
        self.turning_pointsPLot = [self.ax[0][1].plot(self.dataset.samplesToTime(self.turning_points[0].astype(int)) - self.dataset.xRange/2, filterData(self.dataset.yData[0])[self.turning_points[0].astype(int)], marker='.', linestyle='none', markersize=10, color='yellow'),
                                   self.ax[1][1].plot(self.dataset.samplesToTime(self.turning_points[1].astype(int)) - self.dataset.xRange/2, filterData(self.dataset.yData[1])[self.turning_points[1].astype(int)], marker='.', linestyle='none', markersize=10, color='yellow')]

        (popt1, pcov1), proposedParameters1 = fitSinus(dataset.xData, dataset.yData[0], 0.01, self.dataset, self.extremums[0], self.turning_points[0])
        (popt2, pcov2), proposedParameters2 = fitSinus(dataset.xData, dataset.yData[1], 0.01, self.dataset, self.extremums[1], self.turning_points[1])
        self.sinusPredictions = [self.ax[0][1].plot(dataset.xData, sinus(dataset.xData, *popt1), 'r:'), self.ax[1][1].plot(dataset.xData, sinus(dataset.xData, *popt2), 'r:')]
        self.updateTable(popt1, proposedParameters1, self.parameterTables[0])
        self.updateTable(popt2, proposedParameters2, self.parameterTables[1])
        self.popts = [popt1, popt2]

        x0s1 = ((np.linspace(-100, 100, 201) * np.pi) - (self.popts[0][2] % (2*np.pi)))/self.popts[0][1]
        x0s2 = ((np.linspace(-100, 100, 201) * np.pi) - (self.popts[1][2] % (2*np.pi)))/self.popts[1][1]
        
        self.selectablePoints = [self.ax[0][1].plot(x0s1, np.full((201), self.popts[0][3]), 'ro', linestyle='none', markersize=7.5, picker=True, pickradius=7.5),
                                    self.ax[1][1].plot(x0s2, np.full((201), self.popts[1][3]), 'ro', linestyle='none', markersize=7.5, picker=True, pickradius=7.5)]
        
        self.selectablePoints[0][0].id = 0
        self.selectablePoints[1][0].id = 1

        self.selectablePoints[0][0].set(alpha=0)
        self.selectablePoints[1][0].set(alpha=0)

        self.xValues = [0, 0]
        
        self.selectedPoints = [self.ax[0][1].plot(0, 0, marker='o', linestyle='none', color='green', markersize=7.5, zorder=2),
                            self.ax[1][1].plot(0, 0, marker='o', linestyle='none', color='green', markersize=7.5, zorder=2)]
        
        self.selectedPoints[0][0].set(alpha=0)
        self.selectedPoints[1][0].set(alpha=0)
        
        self.fig.canvas.callbacks.connect('pick_event', self.on_pick)

        self.xRangeSliderUpdate(self.ranges[0], 0)
        self.xRangeSliderUpdate(self.ranges[1], 1)

    # updates the plots, when the slider value changes
    def xRangeSliderUpdate(self, value, sliderId):
        self.ranges[sliderId] = value
        minIndex = self.dataset.timeToSamples(value[0], offset=True); maxIndex = self.dataset.timeToSamples(value[1], offset=True)
        xDataRange = self.dataset.xData[minIndex:maxIndex]
        yDataRange = self.dataset.yData[sliderId][minIndex:maxIndex]
        self.ax[sliderId][1].set_xlim(*value)
        self.ax[sliderId][1].set_ylim(np.min(yDataRange), np.max(yDataRange))
        self.plotBorders[sliderId][0].set_xdata(value[0])
        self.plotBorders[sliderId][1].set_xdata(value[1])

        currentRangeInSamples = (self.dataset.timeToSamples(self.ranges[sliderId][0], offset=True), self.dataset.timeToSamples(self.ranges[sliderId][1], offset=True))
        extremumsRange = self.extremums[sliderId][np.where(np.logical_and(self.extremums[sliderId]<currentRangeInSamples[1], self.extremums[sliderId]>currentRangeInSamples[0]))]
        turning_pointsRange = self.turning_points[sliderId][np.where(np.logical_and(self.turning_points[sliderId]<currentRangeInSamples[1], self.turning_points[sliderId]>currentRangeInSamples[0]))]

        (popt, pcov), proposedParameters = fitSinus(xDataRange, yDataRange, 0.01, self.dataset, extremumsRange, turning_pointsRange)
        self.updateTable(popt, proposedParameters, self.parameterTables[sliderId])

        self.sinusPredictions[sliderId][0].set_xdata(xDataRange)
        self.sinusPredictions[sliderId][0].set_ydata(sinus(xDataRange, *popt))
        self.popts[sliderId] = popt

    # update the values in a table
    def updateTable(self, popt, proposedParameters, table):
        table.get_celld()[(1, 0)].get_text().set_text("%.3E" % (popt[0]))
        table.get_celld()[(2, 0)].get_text().set_text("%.3E" % (popt[1]))
        table.get_celld()[(3, 0)].get_text().set_text("%.3E" % (popt[2]%(2*np.pi)))
        table.get_celld()[(4, 0)].get_text().set_text("%.3E" % (popt[3]))
        table.get_celld()[(1, 1)].get_text().set_text("%.3E" % (proposedParameters[0]))
        table.get_celld()[(2, 1)].get_text().set_text("%.3E" % (proposedParameters[1]))
        table.get_celld()[(3, 1)].get_text().set_text("%.3E" % (proposedParameters[2]%(2*np.pi)))
        table.get_celld()[(4, 1)].get_text().set_text("%.3E" % (proposedParameters[3]))
        table.get_celld()[(1, 2)].get_text().set_text("%.3E" % (popt[0] - proposedParameters[0]) + "\n{:.2%}".format(abs(popt[0] - proposedParameters[0])/popt[0]))
        table.get_celld()[(2, 2)].get_text().set_text("%.3E" % (popt[1] - proposedParameters[1]) + "\n{:.2%}".format(abs(popt[1] - proposedParameters[1])/popt[1]))
        table.get_celld()[(3, 2)].get_text().set_text("%.3E" % (popt[2] - proposedParameters[2]) + "\n{:.2%}".format(abs(popt[2]%(2*np.pi) - proposedParameters[2]%(2*np.pi))/popt[2]%(2*np.pi)))
        table.get_celld()[(4, 2)].get_text().set_text("%.3E" % (popt[3] - proposedParameters[3]) + "\n{:.2%}".format(abs(popt[3] - proposedParameters[3])/popt[3]))

    # method, that reacts to selecting a point
    def on_pick(self, event):
        x, y = event.artist.get_xdata(), event.artist.get_ydata()
        ind = event.ind

        self.xValues[event.artist.id] = x[ind[0]]
        self.selectedPoints[event.artist.id][0].set_data(x[ind[0]], y[ind[0]])
        self.selectedPoints[event.artist.id][0].set(alpha=1)
    
    # reveals the time delay
    def displayTimeDelay(self, value):
        if((self.xValues[1] != 0 and self.xValues[0] != 0)):
            self.timeDelayTextBox.set_val(abs(self.xValues[1] - self.xValues[0]))
        else:
            if(self.toggleImportantPointsValue == False):
                self.toggleImportantPoints(0)
            self.errorTextBox.set_val("Select two red points!")

    # toggles between extremums and turning points and the selectable red dots
    def toggleImportantPoints(self, value):
        self.toggleImportantPointsValue = not self.toggleImportantPointsValue

        if(self.toggleImportantPointsValue):
            self.toggleButton.label.set_text("Toggle Extremums/Turning Points")

            self.extremumsPlot[0][0].set(alpha=0)
            self.extremumsPlot[1][0].set(alpha=0)
            self.turning_pointsPLot[0][0].set(alpha=0)
            self.turning_pointsPLot[1][0].set(alpha=0)

            self.selectablePoints[0][0].set(alpha=1)
            self.selectablePoints[1][0].set(alpha=1)

            x0s1 = ((np.linspace(-100, 100, 201) * np.pi) - (self.popts[0][2] % (2*np.pi)))/self.popts[0][1]
            x0s2 = ((np.linspace(-100, 100, 201) * np.pi) - (self.popts[1][2] % (2*np.pi)))/self.popts[1][1]

            self.selectablePoints[0][0].set_data(x0s1, np.full((201), self.popts[0][3]))
            self.selectablePoints[1][0].set_data(x0s2, np.full((201), self.popts[1][3]))

            self.xValues = [0, 0]

            self.selectedPoints[0][0].set(alpha=0)
            self.selectedPoints[1][0].set(alpha=0)

        else:
            self.toggleButton.label.set_text("Toggle Measure Points")
            self.extremumsPlot[0][0].set(alpha=1)
            self.extremumsPlot[1][0].set(alpha=1)
            self.turning_pointsPLot[0][0].set(alpha=1)
            self.turning_pointsPLot[1][0].set(alpha=1)

            self.selectedPoints[0][0].set(alpha=0)
            self.selectedPoints[1][0].set(alpha=0)
            self.selectablePoints[0][0].set(alpha=0)
            self.selectablePoints[1][0].set(alpha=0)

    # clears the plots
    def clear(self):
        self.ax[0][0].lines.clear()
        self.ax[0][0].collections.clear()
        self.ax[0][1].lines.clear()
        self.ax[0][1].collections.clear()
        self.ax[1][0].lines.clear()
        self.ax[1][0].collections.clear()
        self.ax[1][1].lines.clear()
        self.ax[1][1].collections.clear()

    # opens a file dialog and loads the opened file
    def getPath(self, value):
        file_path = tkinter.filedialog.askopenfilenames()
        dataset = Dataset(file_path)
        self.clear()
        self.loadData(dataset)

if __name__ == "__main__":
    # initialize tkinter
    root = tkinter.Tk()
    root.withdraw()

    plot = Plot() # initialize the plot class