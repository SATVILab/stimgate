import os, sys, re
import numpy as np
import matplotlib as mpl
from fcm.graphics import bilinear_interpolate
from fcm.core.transforms import _logicle as logicle

def get_events(fileName,fileList,fcsList):
    '''
    convenience function to fetch events 
    '''
    fileInd = fileList.index(fileName)
    return fcsList[fileInd]

def calculate_fscore(neg_pdf,pos_pdf,beta=0.8,theta=2.0):
    '''
    The f-score is calculated as (precision*recall)/((beta^2*precision)+recall). 
    The TP, FP an FN are estimated from the overlayed probability density functions (pdf) from 
    the positive and negative event distributions.
    '''
    n = len(neg_pdf)
    fpos = np.where(pos_pdf > theta*neg_pdf, pos_pdf-neg_pdf, 0)
    tp = np.array([np.sum(fpos[i:]) for i in range(n)])
    fn = np.array([np.sum(fpos[:i]) for i in range(n)])
    fp = np.array([np.sum(neg_pdf[i:]) for i in range(n)])
    precision = tp/(tp+fp)
    precision[tp==0]=0
    recall = tp/(tp+fn)
    recall[recall==0]=0
    fscores = (1+beta*beta)*(precision*recall)/(beta*beta*precision + recall)
    fscores[np.where(np.isnan(fscores)==True)[0]]=0

    return fscores,precision,recall

def get_positivity_threshold(neg,pos,channelIndex,beta=0.8,theta=2.0, width=10, numBins=None):
    '''
    In order to calculate the f-score the pdfs are found using histogram representations of the 
    data. The number of bins numBins controls how smoothly the pdf fits the actual distribution 
    of events.
    '''

    def move_mean(x, window):
        xs = np.cumsum(x)
        x1 = xs[(window-1):]
        x2 = np.concatenate([[0], xs[:-window]])
        return np.concatenate([[np.nan]*(window-1), (x1-x2)/float(window)])

    neg,pos = neg[:,channelIndex].copy(),pos[:,channelIndex].copy()
    if numBins == None:
        numBins = int(np.sqrt(np.max([neg.shape[0],pos.shape[0]])))

    pdfNeg, bins = np.histogram(neg, bins=numBins, normed=True)
    pdfPos, bins = np.histogram(pos, bins=bins, normed=True)

    pdfNeg = move_mean(pdfNeg, window=width)
    pdfPos = move_mean(pdfPos, window=width)

    xs = (bins[:-1]+bins[1:])/2.0
    fscores,precision,recall = calculate_fscores(pdfNeg,pdfPos,beta=beta,theta=theta)
    fThreshold = xs[np.argmax(fscores)]

    return {'threshold':fThreshold, 'fscores':fscores,'pdfx':xs,'pdfpos':pdfPos,'pdfneg':pdfNeg,
            'precision':precision,'recall':recall}

def get_cytokine_positive_events(cytoIndex,fThreshold,fileList,fcsList):
    '''
    returns the percentages, counts and indices of cytokine positive events
    '''
    
    ## declare variables
    percentages = {}
    counts = {}
    idx = {}
    filterInds = np.array([])

    ## determine and save percentages, counts and indices 
    for fileName in fileList:
        events = get_events(fileName,fileList,fcsList)
        data = events[:,cytoIndex]
        positiveEventInds = np.where(data > fThreshold)[0]

        if events.shape[0] == 0 or len(positiveEventInds) == 0:
            percentages[fileName] = 0.0
        else:
            percentages[fileName] = (float(positiveEventInds.size)/float(events.shape[0])) * 100.0
        counts[fileName] = positiveEventInds.size
        idx[fileName] = positiveEventInds

    return percentages, counts, idx

def set_logicle_transformed_ticks(ax,axis='x',fontsize=10,fontname='sans'):
    '''
    to map an axis to a scale that immunologists are familar with
    '''

    if axis not in ['x','y','both']:
        print "ERROR set_logicle_transformed_ticks: invalid axis arg"
        return None
    
    ## setup scales
    scale = (10**5)*logicle(np.array([0, 100, 10**3, 10**4, 10**5]), 262144, 4.5, None, 0.5)
    tickPairs = [(1,9),(10,90),(100,900),(1000,9000),(10000,90000)]
    minorScale = [(10**5)*logicle(np.linspace(ab[0],ab[1],9),262144,4.5,None,0.5) for ab in tickPairs]
    labels = ['$0$', '$10^2$', '$10^3$', '$10^4$', '$10^5$']
    minorTicks = np.array([])

    for mTicks in minorScale:
        minorTicks = np.hstack([minorTicks,np.array(mTicks)])

    ## format the x axix
    if axis in ['x','both']:
        ax.set_xticks(scale)
        ax.set_xticks(minorTicks,minor=True)
        ax.set_xticklabels(labels,fontsize=fontsize-1,fontname=fontname)
        ax.xaxis.set_ticks_position('bottom')
        ax.set_xlim([0, 1e05])

    ## format the y axis
    if axis in ['y','both']:
        ax.set_yticks(scale)
        ax.set_yticks(minorTicks, minor=True)
        ax.set_yticklabels(labels,fontsize=fontsize-1,fontname=fontname)
        ax.yaxis.set_ticks_position('left')
        ax.set_ylim([0, 1e05])

def set_scatter_ticks(ax,axis,numTicks=6,fontsize=10,fontname='Arial'):
    '''
    formats a scatter axis ticks to the K format
    '''

    if numTicks == 6:
        tickVals = [250000,200000,150000,100000,50000,0]
        tickLabels = ['250K','200K','150K','100K','50K','0']
    elif numTicks == 4:
        tickVals = [250000,150000,50000,0]
        tickLabels = ['250K','150k','50K','0']
    else:
        print "ERROR: set_scatter_ticks: invalid number of ticks"
        return None

    ## format the x axix
    if axis in ['x','both']:
        ax.set_xticks(tickVals)
        ax.set_xticklabels(tickLabels,fontsize=fontsize-1,fontname=fontname)
        ax.xaxis.set_ticks_position('bottom')
        ax.set_xlim([0, 262144])

    ## format the y axis
    if axis in ['y','both']:
        ax.set_yticks(tickVals)
        ax.set_yticklabels(tickLabels,fontsize=fontsize-1,fontname=fontname)
        ax.yaxis.set_ticks_position('left')
        ax.set_ylim([0, 262144])

def create_plot(fig,ax,fileName,cytokine,channelDict,fileList,fcsList,fResults,cpe,counts,idx):

    print 'creating plot', cytokine

    def add_text(ax,txt,xAxLimit,yAxLimit,transform,isTitle=False):
        xPos = 0.5 * (xAxLimit[1] - xAxLimit[0])
        yPos = 0.94 * (yAxLimit[1] - yAxLimit[0])

        if isTitle == True:
            yPos = 1.97 * (yAxLimit[1] - yAxLimit[0])

        fColor = 'black'
        tColor = 'white'
                    
        ax.text(xPos,yPos,txt,color=tColor,fontsize=fontSize-1,
                ha="center", va="top",fontname=fontName,
                bbox = dict(boxstyle="round",facecolor=fColor,alpha=0.8)
                )

    def add_percent(ax,xPos,yPos,color,percent):
        percentToShow = "%2.2f"%round(percent,3)+"%"
        ax.text(xPos,yPos,percentToShow,color='white',fontsize=fontSize,
                ha="left", va="top",fontname=fontName,
                bbox = dict(boxstyle="round",facecolor=color,alpha=0.8)
                )
        
    ## variables
    fontSize = 10
    fontName = 'arial'
    bgEvents = get_events(fileName,fileList,fcsList)
    chan1 = cytokine
    chan2 = "SSC"
    myCmap = mpl.cm.gist_heat
    cytokineNames = {'TNFA':r"TNF-$\alpha$ Alexa 700",
                     'IL2':r"IL-2 PE",
                     'IFNG':r"IFN-$\gamma$ PE-Cy7",
                     }

    ## plot the background
    dataX,dataY = (bgEvents[:,channelDict[chan1]],bgEvents[:,channelDict[chan2]])
    borderEventsX1 = np.where(dataX == 0)[0]
    borderEventsX2 = np.where(dataY == dataX.max())[0]
    borderEventsY1 = np.where(dataY == 0)[0]
    borderEventsY2 = np.where(dataY == dataY.max())[0]
    borderEventsX = np.hstack([borderEventsX1,borderEventsX2])
    borderEventsY = np.hstack([borderEventsY1,borderEventsY2])
    borderEvents = np.hstack([borderEventsX,borderEventsY])
    nonBorderEvents = np.array(list(set(range(dataX.size)).difference(set(borderEvents))))
    colorList = bilinear_interpolate(dataX[nonBorderEvents],dataY[nonBorderEvents],bins=50)
    ax.scatter([dataX[nonBorderEvents]],[dataY[nonBorderEvents]],c=colorList,s=1,edgecolor='none',cmap=myCmap)
    
    ## add percentage box
    xAxLimit = [0,1e05]
    yAxLimit = [0,262144]
    xPos = 0.1 * (xAxLimit[1] - xAxLimit[0])
    yPos = 0.94 * (yAxLimit[1] - yAxLimit[0])
    add_percent(ax,xPos,yPos,'black',cpe[fileName])
    
    ## plot the fbeta threshold
    p = ax.plot(np.array(fResults['threshold']).repeat(50),np.linspace(yAxLimit[0],yAxLimit[1],50),color='black',linestyle='-',linewidth=1.0)

    ## axes and labels
    ax.set_title(re.sub("\.fcs","",fileName),fontsize=fontSize,fontname=fontName)
    set_logicle_transformed_ticks(ax,fontsize=fontSize,fontname=fontName)
    set_scatter_ticks(ax,'y',fontsize=fontSize,fontname=fontName)
    ax.set_xlabel(cytokineNames[cytokine.upper()],fontsize=fontSize-1,fontname=fontName)
    ax.set_ylabel(fcsList[0].channels[channelDict['SSC']],fontsize=fontSize-1,fontname=fontName)
    ax.set_aspect(1./ax.get_data_ratio())


def save_labels(labels,projectDir,fileName,labelsID):
    """
    saves clustering labels

    """

    if type(labels) == type([]):
        labels = np.array(labels)

    saveFilePath = os.path.join(projectDir,fileName+"_%s"%(labelsID)+".npy")
    np.save(saveFilePath,labels)

def load_labels(projectDir,fileName,labelsID):
    """
    loads clustering labels

    """

    saveFilePath = os.path.join(projectDir,fileName+"_%s"%(labelsID)+".npy")
    if not os.path.exists(saveFilePath):
        return None

    labels = np.load(saveFilePath)

    return labels