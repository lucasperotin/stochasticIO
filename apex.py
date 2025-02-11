import random
import math
import functools
import logging
from PIL import Image, ImageDraw
from itertools import groupby
from operator import itemgetter
import sys
import time
import argparse
import pickle
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from bitarray import bitarray, frozenbitarray

# This is the standard deviation of the gaussian noise we add to number of nodes and durations
STDEV=0.1
# This is the high pass filter for the noise. 0.8 means that we reject anything lower than 0.8 the average value
MAXLOW=0.95
# This is how far away we accept to be from the distribution of application type. 0.1 means 10%
MAXDELTA=0.1

FACTORS = {
    'nodes': 4,
    'cores_per_node': 10,
    'ram_per_node': 10,
    'node_bw': 1,
    'system_bw': 1
}

class AppClass:
    def __init__(self, **kwargs):
        for arg in ["name", "corehour_ratio", "walltime", "cores", "input_ratio", "output_ratio", "ckpt_ratio", "iterations", "ensemble", "ckpt_period"]:
            value = kwargs.pop(arg, None)
            if None is value:
                raise ValueError("Invalid initialization of AppClass: %s is not provided" % (arg,))
            setattr(self, arg, value)

        if kwargs:
            invalid_args = ", ".join(kwargs)
            raise ValueError("Invalid keyword argument(s) for AppClass creation: %s" % (invalid_args,))

    def __repr__(self):
        return f"AppClass(name={self.name}, corehour_ratio={self.corehour_ratio}, walltime={self.walltime}, cores={self.cores}, input/output/ckpt_ratio={self.input_ratio}/{self.output_ratio}/{self.ckpt_ratio}, iterations/ensemble={self.iterations}/{self.ensemble}, ckpt_period={self.ckpt_period})"

class Platform:
    def __init__(self, **kwargs):
        for arg in ["name", "nodes", "cores_per_node", "ram_per_node", "node_bw", "system_bw", "workload"]:
            value = kwargs.pop(arg, None)
            if None is value:
                raise ValueError("Invalid initialization of Platform: %s is not provided" % (arg,))
            setattr(self, arg, value)

        if kwargs:
            invalid_args = ", ".join(kwargs)
            raise ValueError("Invalid keyword argument(s) for Platform creation: %s" % (invalid_args,))

    def scale(self, d):
        for arg in ["nodes", "cores_per_node", "ram_per_node", "node_bw", "system_bw"]:
            value = d.pop(arg, None)
            if None is value:
                continue
            setattr(self, arg, int(math.ceil(getattr(self, arg) * float(value))))

    def __repr__(self):
        return f"Platform(name={self.name}, nodes={self.nodes}, cores_per_node={self.cores_per_node}, ram_per_node={self.ram_per_node}, node/system_bw={self.node_bw}/{self.system_bw}, workload={repr(self.workload)})"

debug = Platform(name = 'debug',
                 nodes = 10, 
                 cores_per_node = 1, 
                 ram_per_node = 1, 
                 node_bw = 1, 
                 system_bw = 2, 
                 workload = [AppClass(name = 'A', corehour_ratio = 0.5, walltime = 3600, cores = 2, input_ratio = 1, output_ratio = 1, ckpt_ratio = 0.3, iterations = 3, ensemble = 1, ckpt_period = 600),
                             AppClass(name = 'B', corehour_ratio = 0.5, walltime = 7200, cores = 4, input_ratio = 0.5, output_ratio = 0.5, ckpt_ratio = 0.3, iterations = 3, ensemble = 1, ckpt_period = 1200)])

trilab = Platform(name = 'trilab', nodes = 8894, cores_per_node = 16, ram_per_node = 32*1024*1024*1024, node_bw = 1.2*1024*1024*1024, system_bw = 200*1024*1024*1024, 
                  workload = [AppClass(name = "EAP", corehour_ratio = 0.2, walltime = 262.4*3600, cores = 16384, input_ratio = 0.03, output_ratio = 1.05, ckpt_ratio = 0.3, iterations = 30, ensemble = 1, ckpt_period = 3600),
                              AppClass(name = "LAP", corehour_ratio = 0.02, walltime = 64*3600, cores = 4096, input_ratio = 0.05, output_ratio = 3.2, ckpt_ratio = 0.75, iterations = 10, ensemble = 1, ckpt_period = 3600),
                              AppClass(name = "Silverton", corehour_ratio = 0.05, walltime = 128*3600, cores = 32768, input_ratio = 0, output_ratio = 0.43, ckpt_ratio = 2.1, iterations = 6, ensemble = 1, ckpt_period = 3600),
                              AppClass(name = "VPIC", corehour_ratio = 0.03, walltime = 157.2*3600, cores = 30000, input_ratio = 0.05, output_ratio = 2.7, ckpt_ratio = 0.1875, iterations = 4, ensemble = 1, ckpt_period = 3600),
                              AppClass(name = "DakotaA", corehour_ratio = 0.03, walltime = 100*3600, cores = 8192, input_ratio = 0.05, output_ratio = 0.085, ckpt_ratio = 0.002, iterations = 10, ensemble = 100, ckpt_period = 3600,),
                              AppClass(name = "DakotaS", corehour_ratio = 0.03, walltime = 100*3600, cores = 4096, input_ratio = 0.005, output_ratio = 0.0204, ckpt_ratio = 0.3, iterations = 30, ensemble = 300, ckpt_period = 3600)])

nersc = Platform(name = 'nersc', nodes = 8894, cores_per_node = 16, ram_per_node = 32*1024*1024*1024, node_bw = 1.2*1024*1024*1024, system_bw = 200*1024*1024*1024,
    workload = [AppClass(name = "ALS", corehour_ratio = 0.001, walltime = 8*0.5*3600, cores = 100, input_ratio = 0.3165, output_ratio = 1.0629, ckpt_ratio = 0, iterations = 10760/8, ensemble = 1, ckpt_period = 3600),
                AppClass(name = "CESM", corehour_ratio = 0.06, walltime = 240*3600, cores = 8000, input_ratio = 0, output_ratio = 10.2795, ckpt_ratio = 0.32, iterations = 8, ensemble = 1, ckpt_period = 3600),
                AppClass(name = "GTS", corehour_ratio = 0.06, walltime = 48*3600, cores = 16512, input_ratio = 0, output_ratio = 0.1478, ckpt_ratio = 0.67, iterations = 36, ensemble = 1, ckpt_period = 3600),
                AppClass(name = "HipMer", corehour_ratio = 0.001, walltime = 4*3600, cores = 960, input_ratio = 0.6583, output_ratio = 0.3434, ckpt_ratio = 0, iterations = 100, ensemble = 1, ckpt_period = 3600),
                AppClass(name = "Materials", corehour_ratio = 0.19, walltime = 41.7*3600, cores = 2400, input_ratio = 0, output_ratio = .2708, ckpt_ratio = 0.5, iterations = 100, ensemble = 1, ckpt_period = 3600),
                AppClass(name = "MILC", corehour_ratio = 0.11, walltime = 6*3600, cores = 4096, input_ratio = 0.7388, output_ratio = 0.7388, ckpt_ratio = 0, iterations = 1000, ensemble = 1, ckpt_period = 3600),
                AppClass(name = "SkySurvey", corehour_ratio = 0.06, walltime = 2*4*3600, cores = 16*24, input_ratio = 0.0879, output_ratio = 0.0269, ckpt_ratio = 0, iterations = 21000/16/2, ensemble = 1, ckpt_period = 3600)])

class  TimeEvent:
    def __init__(self, **kwargs):
        for arg in ["start", "duration", "freenodes"]:
            value = kwargs.pop(arg, None)
            if None is value:
                raise ValueError("Invalid initialization of TimeEvent: %s is not provided" % (arg,))
            setattr(self, arg, value)

        if kwargs:
            invalid_args = ", ".join(kwargs)
            raise ValueError("Invalid keyword argument(s) for TimeEvent creation: %s" % (invalid_args,))

    def __repr__(self):
        if self.freenodes.count() > 5:
            return f"(start={self.start}, duration={self.duration}, freenodes=|{self.freenodes.count()}|[...])"
        else:
            return f"(start={self.start}, duration={self.duration}, freenodes=|{self.freenodes.count()}|{repr(self.freenodes)})"

class Schedule:
    def __init__(self, **kwargs):
        for arg in ["maxEndTime", "minStartTime", "applist", "timeline", "workload_corehour", "tot_corehour", "app_corehour", "appnb"]:
            value = kwargs.pop(arg, None)
            if None is value:
                raise ValueError("Invalid initialization of Schedule: %s is not provided" % (arg,))
            setattr(self, arg, value)

        if kwargs:
            invalid_args = ", ".join(kwargs)
            raise ValueError("Invalid keyword argument(s) for Schedule creation: %s" % (invalid_args,))

    def __repr__(self):
        return f"Schedule(maxEnd/minStartTime = {self.maxEndTime}/{self.minStartTime}, workload/tot/app_corehour={self.workload_corehour}/{self.tot_corehour}/{self.app_corehour}, applist={repr(self.applist)}, timeline={repr(self.timeline)})"

class AppInstance:
    def __init__(self, **kwargs):
        for arg in ['idx', 'nbNodes', 'duration', 'start', 'nodes']:
            value = kwargs.pop(arg, None)
            if None is value:
                raise ValueError("Invalid initialization of AppInstance: %s is not provided" % (arg,))
            setattr(self, arg, value)

        if kwargs:
            invalid_args = ", ".join(kwargs)
            raise ValueError("Invalid keyword argument(s) for AppInstance creation: %s" % (invalid_args,))

    def __repr__(self):
        if self.nbNodes > 5:
            return f"App(idx = {self.idx}, nbNodes = {self.nbNodes}, start = {self.start}, duration = {self.duration}, nodes = [...])"
        return f"App(idx = {self.idx}, nbNodes = {self.nbNodes}, start = {self.start}, duration = {self.duration}, nodes = [{repr(self.nodes)}])"

def noisy(v):
    x = random.gauss(v, STDEV*v)
    if x <= MAXLOW*v:
        x = MAXLOW*v
    return x

def setToRanges(data):
    ranges =[]
    for k,g in groupby(enumerate(data),lambda x:x[0]-x[1]):
        group = (map(itemgetter(1),g))
        group = list(map(int,group))
        ranges.append((group[0],group[-1]))
    return ranges

def inputDuration(plat, app, nbNodes):
    ioBW = min(nbNodes*plat.node_bw, plat.system_bw)
    return app.input_ratio*nbNodes*plat.ram_per_node/ioBW

def outputDuration(plat, app, nbNodes):
    ioBW = min(nbNodes*plat.node_bw, plat.system_bw)
    return app.output_ratio*nbNodes*plat.ram_per_node/ioBW

def ckptDuration(plat, app, nbNodes):
    ioBW = min(nbNodes*plat.node_bw, plat.system_bw)
    return app.ckpt_ratio*nbNodes*plat.ram_per_node/ioBW

def nbCkpt(plat, app, nbNodes, walltime):
    inputD = inputDuration(plat, app, nbNodes)
    outputD = outputDuration(plat, app, nbNodes)
    if inputD + outputD > walltime:
        logging.error(f"Application {app.name} has a walltime of {walltime} which is less than it's I/O time ({inputD}+{outputD}+{ckptD}")
        raise Exception("Impossible to run this application in the given walltime")
    ckptD = ckptDuration(plat, app, nbNodes)
    if ckptD == 0:
        return 0
    if ckptD > app.ckpt_period:
        logging.error(f"Application {app.name} with {nbNodes} nodes and a walltime of {walltime} tries to checkpoint every {app.ckpt_period}, but its checkpoint time is {ckptD} which is longer than the period")
        raise Exception("Impossible to checkpoint with this period")
    nb = math.floor((walltime-inputD-outputD)/ckptD)
    if nb <= 1:
        return 0
    return nb

def compDuration(plat, app, nbNodes, walltime):
    inputD = inputDuration(plat, app, nbNodes)
    outputD = outputDuration(plat, app, nbNodes)
    if inputD + outputD >= walltime:
        logging.error(f"Application {app.name} has a walltime of {walltime} which is less than it's I/O time ({inputD}+{outputD}")
        raise Exception("Impossible to run this application in the given walltime")
    ckptD = ckptDuration(plat, app, nbNodes)
    nb = nbCkpt(plat, app, nbNodes, walltime)
#    if inputD + outputD + nb*ckptD >= walltime/2:
#        logging.warning(f"Application {app.name} has a walltime of {walltime} which is less than twice it's I/O time ({inputD}+{outputD}+{nb}*{ckptD})")
    return walltime - inputD - outputD - nb*ckptD

def drawRealisticSchedule(plat, sched, minTime, maxTime, width, height, filename):
    nodeHeight = height/plat.nodes
    if nodeHeight < 1:
        raise Exception(f"Image format does not give at least one pixel in height per node. Require at least {plat.nodes} in height")
    secondWidth = width/(maxTime-minTime)
    nbEmptyApps = 0
    im = Image.new(mode="RGB", size=(math.ceil(width), math.ceil(height)))
    im1 = ImageDraw.Draw(im)

    appColors = [ '#9e0142', '#d53e4f', '#f46d43', '#fdae61', '#fee08b', '#e6f598', '#abdda4', '#66c2a5', '#66c2a5', '#5e4fa2' ]

    for e in sched.timeline:
        e_start = math.floor( (e.start-minTime) * secondWidth)
        if e_start >= width:
            continue
        if e_start < 0:
            e_start = 0
        if e.duration == float('inf'):
            e_end = width
        else:
            e_end = math.floor( (e.start + e.duration - minTime) * secondWidth )
            if e_end > width:
                e_end = width
            if e_end < 0:
                continue
        for r in setToRanges(e.freenodes):
            node_start = math.floor(min(r)*nodeHeight)
            node_end = math.ceil((max(r)+1)*nodeHeight)-1
            shape=[(e_start, node_start), (e_end, node_end)]
            im1.rectangle(shape, fill = '#ffffff', outline = '#ffffff')

    for app in sched.applist:
        app_start = math.floor( (app.start-minTime) * secondWidth)
        if app_start >= width:
            continue
        if app_start < 0:
            app_start = 0
        app_end = math.floor( (app.start + app.duration - minTime) * secondWidth )
        if app_end > width:
            app_end = width
        if app_end < 0:
            continue
        if not (app.start <= sched.maxEndTime and app.start + app.duration <= sched.maxEndTime):
            logging.error(f"app.start = {app.start}, sched.maxEndTime = {sched.maxEndTime}, app.duration = {app.duration}")
            assert(False)
        if not (app_end <= width and app_start <= width):
            logging.error(f"app_end = {app_end}, width = {width}, app_start = {app_start}")
            assert(False)
        if app_end <= app_start + 1:
            nbEmptyApps += 1
        for r in setToRanges(app.nodes):
            node_start = math.floor(min(r)*nodeHeight)
            node_end = math.ceil((max(r)+1)*nodeHeight)-1
            shape=[(app_start, node_start), (app_end, node_end)]
            im1.rectangle(shape, fill = appColors[ app.idx ], outline = appColors[ app.idx ])
    
    im.save(filename)
    if nbEmptyApps > 0:
        logging.warning(f"There was {nbEmptyApps} applications that lasted at most a pixel long, they might be hidden by other apps in the gantt chart of {filename}")

def selectAppEx(plat, sched):
    exweights = [plat.weights[i] - sched.appnb[i] for i in range(len(plat.workload))]
    exweights = [0 if e < 0 else e for e in exweights]
    if sum(exweights) == 0:
        exweights = plat.weights
    sels = random.choices(range(len(plat.workload)), weights=exweights)
    sel = sels[0]
    return sel

def realisticSchedIsBalanced(plat, sched):
    appnb = sched.appnb
    weights = plat.weights

    for i in range(len(plat.workload)):
        app = plat.workload[i]
        #print(f"balance test: apps named {app.name} has {appnb[i]} instances, expect to have {weights[i]}")
        if abs( appnb[i] - weights[i] ) > MAXDELTA * appnb[i]:
            return False
    return True

def timelineToStr(timeline, maxset=0):
    s = "["
    for e in timeline:
        if e.freenodes.count() <= maxset:
            s += f"['start': {e.start}, 'duration': {e.duration}, 'freenodes': {repr(e.freenodes)}] "
        else:
            s += f"['start': {e.start}, 'duration': {e.duration}, 'freenodes': |{e.freenodes.count()}|] "
    s += "]"
    return s

def timelineIsConsistent(sched):
    for i in range(len(sched.timeline)):
        e = sched.timeline[i]
        if e.duration == float('inf') and i < len(sched.timeline)-1:
            logging.error(f"Timeline is inconsistent: event {e} at index {i} has infinite duration but is not last event of timeline {timelineToStr(sched.timeline)}")
            return False
        if i == 0:
            nt = e.start + e.duration
        else:
            if abs(e.start - nt) > 1:
                logging.error(f"Timeline is inconsistent: event {e} does not start at {nt} in {timelineToStr(sched.timeline)}")
                return False
            nt += e.duration
    if nt != float('inf'):
        logging.error(f"Timeline is inconsistent: does not end with an event of infinite duration {timelineToStr(sched.timeline)}")
        return False
    return True

def consolidateTimeline(sched, plat, maxTime):
    minStartTime = 0
    min_posStart = len(sched.timeline)-1
    nbtest = 0
    for app in sorted(plat.workload, key=lambda d: d.cores):
        nbtest += 1
        minNodes = math.ceil(MAXLOW*app.cores/plat.cores_per_node)
        minDuration = MAXLOW*app.walltime
        (canFit, pos_start, pos_end, all_freenodes) = appFitsFirstFit(sched, minNodes, minDuration, maxTime)
        if canFit:
            if pos_start < min_posStart:
                min_posStart = pos_start
        if min_posStart <= 1:
            break
    if min_posStart >= 1:
        del sched.timeline[0:min_posStart]
        logging.info(f"Optimization: removed {min_posStart} events")
        logging.info(f"Remaining {len(sched.timeline)} events in the timeline. First event start time is {sched.timeline[0].start}s")
        minStartTime = sched.timeline[0].start
        if logging.getLogger().isEnabledFor(logging.DEBUG):
            if not timelineIsConsistent(sched):
                sys.exit(1)
    else:
        logging.info(f"No merge (tested {nbtest} apps) -- remaining {len(sched.timeline)} events in the timeline. First event duration starts at {sched.timeline[0].start}")

    if False:
        mergedSome = True
        totMerged = 0
        before=len(sched.timeline)
        while mergedSome:
            mergedSome = False
            for p in range(len(sched.timeline)-1):
                if sched.timeline[p].freenodes == sched.timeline[p+1].freenodes:
                    sched.timeline[p].duration += sched.timeline[p+1].duration
                    del sched.timeline[p+1]
                    mergedSome = True
                    totMerged += 1
                    break
        mergedSome = True
        logging.info(f"Optimization: merged {totMerged} events during 2nd phase")

    return minStartTime

def appFitsFirstFit(sched, nbNodes, duration, maxTime):
    for pos_start in range(len(sched.timeline)):
        start = sched.timeline[pos_start].start
        if start > maxTime:
            return (False, None, None, None)
        if logging.getLogger().isEnabledFor(logging.DEBUG):
            logging.debug(f"Trying to start at phase {pos_start} ({start} s)")
        all_freenodes = sched.timeline[pos_start].freenodes.copy()
        tot_duration = 0
        for pos_end in range(pos_start, len(sched.timeline)):
            if logging.getLogger().isEnabledFor(logging.DEBUG):
                logging.debug(f"Considering {pos_start} to {pos_end} the latter starting at {sched.timeline[pos_end].start} and lasting {sched.timeline[pos_end].duration}")
            all_freenodes &= sched.timeline[pos_end].freenodes
            if logging.getLogger().isEnabledFor(logging.DEBUG):
                logging.debug(f"Freenodes in this range is of len {all_freenodes.count(1)}")
            if all_freenodes.count(1) < nbNodes:
                if logging.getLogger().isEnabledFor(logging.DEBUG):
                    logging.debug(f"It is smaller than {nbNodes}, need to consider another start")
                break
            tot_duration += sched.timeline[pos_end].duration
            if logging.getLogger().isEnabledFor(logging.DEBUG):
                logging.debug(f"Total duration is {tot_duration}")
            if tot_duration >= duration:
                return (True, pos_start, pos_end, all_freenodes)
    logging.error(f"There should always be a segment of infinite length at the end of the timeline. Can an application of {nbNodes} nodes fit on the platform?")
    raise Exception("Application does not fit on machine")

def insertApp(plat, sched, appIdx, nbNodes, duration, freenodes, pos_start, pos_end, maxEndTime, minStartTime, maxTime):
    start = sched.timeline[pos_start].start
    appFailedToInsertInRow = 0
    if logging.getLogger().isEnabledFor(logging.DEBUG):
        logging.debug(f"This interval is a good candidate")

    fnStart = 0
    appNodes = []
    revAppNodesArray = bitarray(plat.nodes)
    revAppNodesArray.setall(1)
    it = freenodes.itersearch(1)
    for i in range(nbNodes):
        n = next(it)
        appNodes.append(n)
        revAppNodesArray[n] = False
    if logging.getLogger().isEnabledFor(logging.DEBUG):
        logging.debug(f"Selected {len(appNodes)} nodes from free nodes")
            
    prev_duration = 0
    for p in range(pos_start, pos_end):
        plen = sched.timeline[p].freenodes.count(1)
        if logging.getLogger().isEnabledFor(logging.DEBUG):
            logging.debug(f"Removing {revAppNodesArray.count(0)} available nodes from interval {p} -- was of length {sched.timeline[p].freenodes.count(1)}")
        sched.timeline[p].freenodes &= revAppNodesArray
        if logging.getLogger().isEnabledFor(logging.DEBUG):
            logging.debug(f"Is now of length {sched.timeline[p].freenodes.count(1)}")
        prev_duration += sched.timeline[p].duration

    if logging.getLogger().isEnabledFor(logging.DEBUG):
        logging.debug(f"Adjusting last phase: new phase would last {sched.timeline[pos_end].duration-(duration-prev_duration)} s")
    end = start + duration
    if sched.timeline[pos_end].duration - (duration-prev_duration) > 0:
        new_duration = prev_duration + sched.timeline[pos_end].duration - duration
        if logging.getLogger().isEnabledFor(logging.DEBUG):
            logging.debug(f"Inserting new phase at time {end} of duration {new_duration}")
        sched.timeline.insert(pos_end+1, TimeEvent(start = end, duration = new_duration, freenodes = sched.timeline[pos_end].freenodes.copy()))
            
        if logging.getLogger().isEnabledFor(logging.DEBUG):
            logging.debug(f"Updating interval {pos_end}: new duration is {end-sched.timeline[pos_end].start}")
        sched.timeline[pos_end].duration = end-sched.timeline[pos_end].start
    elif sched.timeline[pos_end].duration - (duration-prev_duration) == 0:
        if logging.getLogger().isEnabledFor(logging.DEBUG):
            logging.debug(f"both phases end at exactly the same time, just remove nodes from the existing phase")
    else:
        logging.error(f"new phase would be negative! Duration of last event is {sched.timeline[pos_end].duration}, overall duration is {duration} and duration up to now is {prev_duration}")
        assert(False)

    if logging.getLogger().isEnabledFor(logging.DEBUG):
        logging.debug(f"Removing {len(appNodes)} available nodes from interval {pos_end} -- was of length {sched.timeline[pos_end].freenodes.count()}")
    sched.timeline[pos_end].freenodes &= revAppNodesArray
    if logging.getLogger().isEnabledFor(logging.DEBUG):
        logging.debug(f"Is now of length {sched.timeline[pos_end].freenodes.count()}")
            
    if logging.getLogger().isEnabledFor(logging.DEBUG):
        logging.debug(f"Application inserted at date {start}")
    sched.applist.append( AppInstance(idx = appIdx, nbNodes = nbNodes, duration = duration, start = start, nodes = appNodes) )
    if end > maxEndTime:
        maxEndTime = end
        sched.maxEndTime = maxEndTime
        sched.tot_corehour = end * plat.nodes * plat.cores_per_node
    sched.workload_corehour[appIdx] += nbNodes * plat.cores_per_node * duration
    sched.app_corehour += nbNodes * plat.cores_per_node * duration
    sched.appnb[appIdx] += 1
    return maxEndTime

statStartSelect = None
statStartLookup = None
statStartInsert = None
statDone        = None
stepIdx         = None

def stepRealisticSchedule(plat, sched, maxEndTime, minStartTime, maxTime):
    # Pick an app
    statStartSelect = time.time()
    appIdx = selectAppEx(plat, sched)
    app = plat.workload[appIdx]
    nbNodes = math.ceil(math.ceil(noisy(app.cores))/plat.cores_per_node)
    duration = math.ceil(noisy(app.walltime))
    # call compDuration to capture issues like walltime too short early
    compDuration(plat, app, nbNodes, duration)
    if logging.getLogger().isEnabledFor(logging.DEBUG):
        logging.debug(f"\n\nStep {stepIdx}: finding a spot for {app.name} of {nbNodes} nodes, {duration} s total")
    statStartLookup = time.time()
    # Find a slot
    (appInserted, pos_start, pos_end, freenodes) = appFitsFirstFit(sched, nbNodes, duration, maxTime)
    statStartInsert = time.time()
    if appInserted:
        maxEndTime = insertApp(plat, sched, appIdx, nbNodes, duration, freenodes, pos_start, pos_end, maxEndTime, minStartTime, maxTime)
    statDone = time.time()
    return (appInserted, app, maxEndTime)

def createRealisticSchedule(plat, minTime, maxTime):
    weights = [app.corehour_ratio*(plat.nodes*maxTime)/(app.walltime*math.ceil(app.cores/plat.cores_per_node)) for app in plat.workload]
    plat.weights = weights
    logging.info(f"Goal for number of applications over maxTime: {repr(weights)}")
    initialFreeNodes = bitarray(plat.nodes)
    initialFreeNodes.setall(1)
    e = TimeEvent(start = 0, duration = float('inf'), freenodes = initialFreeNodes)
    sched = Schedule(maxEndTime = 0, minStartTime = minTime, applist = [], timeline = [ e ],
                     workload_corehour = [0.0 for x in range(len(plat.workload))], tot_corehour =  0.0, app_corehour = 0.0, appnb = [0 for x in range(len(plat.workload))])
    maxEndTime = 0
    minStartTime = 0
    info_last = time.time()

    statStartSelect = 0
    statStartLookup = 0
    statStartInsert = 0
    statDone        = 0
    sumSelect  = 0
    sumLookup  = 0
    sumInsert  = 0
    nbInserted = 0
    nbFailedToInsert = 0

    stepIdx=0
    appFailedToInsertInRow = 0
    while appFailedToInsertInRow < 1000 and not (realisticSchedIsBalanced(plat, sched) and minStartTime >= minTime):
        if time.time() - info_last > 5:
            info_last = time.time()
            logging.info(f"")
            minStartTime = consolidateTimeline(sched, plat, maxTime)
            if sched.tot_corehour > 0:
                occ = sched.app_corehour / sched.tot_corehour
                logging.info(f"occupation: {occ} workload: {[ (sched.workload_corehour[x]/sched.app_corehour) for x in range(len(sched.workload_corehour))]}")
                goals = [f"{sched.appnb[i]}({sched.appnb[i]/plat.weights[i]}) " for i in range(len(plat.weights))]
                logging.info(f"#apps(goal): {goals}")
            if nbInserted+nbFailedToInsert > 0:
                logging.info(f"Timings ({nbInserted} tasks inserted, {nbFailedToInsert} tasks not inserted because of deadline): TaskSelection: {sumSelect/(nbInserted+nbFailedToInsert)}s FirstFit Lookup: {sumLookup/(nbInserted+nbFailedToInsert)}s TaskInsertion: {sumInsert/(nbInserted+nbFailedToInsert)}s")
            else:
                logging.info(f"No Timings: no task processed since last time")
            nbInserted = 0
            nbFailedToInsert = 0
            sumInsert = 0
            sumLookup = 0
            sumSelect = 0
            logging.info(f"maxEndTime is now {maxEndTime} ({maxEndTime/minTime} of minimum simulation, {maxEndTime/maxTime} of maximum simulation)")
            logging.info(f"minStartTime is now {minStartTime} ({minStartTime/minTime} of minimum simulation, {minStartTime/maxTime} of maximum simulation)")

        (appInserted, app, maxEndTime) = stepRealisticSchedule(plat, sched, maxEndTime, minStartTime, maxTime)

        if logging.getLogger().isEnabledFor(logging.DEBUG):
            assert(maxEndTime > 0)
            sched.maxEndTime = maxEndTime
            with open(f"debug-{stepIdx}.txt", 'wb') as f:
                pickle.dump([plat, sched], f)
            if not timelineIsConsistent(sched):
                sys.exit(1)

        stepIdx += 1
        sumSelect += statStartLookup - statStartSelect
        sumLookup += statStartInsert - statStartLookup
        sumInsert += statDone - statStartInsert
        if not appInserted:
            if logging.getLogger().isEnabledFor(logging.DEBUG):
                logging.debug(f"Not inserting application {app.name} ({appFailedToInsertInRow} in a row) because it would start after the maximal time required ({maxTime})")
            appFailedToInsertInRow += 1
            nbFailedToInsert += 1
            #raise Exception(f"Could not insert app {app.name} with {nbNodes} nodes and that last {duration} seconds, but that should be impossible as time is unbounded")
        else:
            nbInserted += 1

    if appFailedToInsertInRow >= 1000:
        logging.info(f"Stop scheduling: couldn't insert {appFailedToInsertInRow} applications in a row")
    if minStartTime >= minTime:
        logging.info(f"Stop scheduling: first possible task starts at {minStartTime/minTime} of the minimum time ({minTime})")
    if realisticSchedIsBalanced(plat, sched):
        logging.info(f"Done scheduling: scheduling is currently a realistic balance")

    sched.maxEndTime = maxEndTime
    return sched

def init_argparse() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        usage="%(prog)s [OPTION] [FILE]...",
        prog="apex.py",
        description="Create a schedule from the APEX data, and output it to a file, or read this file to\n"
                    "create windows in which no application leaves or join, and output a scenario for\n"
                    "the simulation of that window."

    )
    parser.add_argument("-v", "--version", action="version", version = f"{parser.prog} version 1.0.0")
    parser.add_argument("-d", "--debug", action="store_true", help="Enable debug mode")
    parser.add_argument("-m", "--mode", action="store", choices = ["schedule", "draw", "window"], help="Decide if we generate a schedule or scenarios from a schedule", default="schedule")

    parser.add_argument("-p", "--platform", action="store", choices = ["trilab", "nersc", "debug"], help="the platform to consider (only in scheduler mode)", default="trilab")
    parser.add_argument("--min", action="store", help="Minimum duration (in days) of the schedule plan (only in scheduler mode)", default=30*3)
    parser.add_argument("--max", action="store", help="Maximum duration (in days) of the schedule plan (only in scheduler mode)", default=30*6)

    parser.add_argument("-i", "--input", action="store", help="Input file (only in window and draw mode)", default="-")

    parser.add_argument("-W", "--width", action="store", help="Image width (only in draw mode)", default=2048)
    parser.add_argument("--minTime", action="store", help="Ignore events before this date (draw and window modes)", default=0)
    parser.add_argument("--maxTime", action="store", help="Ignore events after this date (draw and window modes) -- possible values include a time, minTime or maxTime", default="minTime")

    parser.add_argument("--window", action="store", help="Minimum duration of the window (only in window mode). Set 0 to generate an histogram of window sizes", default="0")

    parser.add_argument("-o", "--output", action="store", help="Output file", default="-")
    return parser

def windowEvents(t, IOName, workd, iod, windowStart, windowEnd, app, a):
    wStart = t
    wEnd = t+workd
    ioStart = wEnd
    ioEnd = ioStart+iod
    if ioStart <= windowEnd and ioEnd >= windowStart:
        #IO is in the window, we need to output [something, IOtime]
        winIOStart = max(ioStart, windowStart)
        winIOEnd = min(ioEnd, windowEnd)
        winIOD = winIOEnd - winIOStart
        if wStart <= windowEnd and wEnd >= windowStart:
            #Work is also in the window, easy case.
            winWStart = max(wStart, windowStart)
            winWEnd = min(wEnd, windowEnd)
            winWD = winWEnd - winWStart
        else:
            #(preceding) work is not in the window, we add a leading 0
            winWD = 0
        if logging.getLogger().isEnabledFor(logging.DEBUG):
            logging.debug(f"app {app.name} starting at {a.start} has Work at time {wStart} for {workd} then {IOName} at time {ioStart} for {iod}. Some overlapping the window [{windowStart}; {windowEnd}] for [{winWD},{winIOD}]")
        return [winWD, winIOD]
    #IO (which happens after the work) is not in the window. Maybe the work is
    if wStart <= windowEnd and wEnd >= windowStart:
        #There is only work, and then we're done.
        winWStart = max(wStart, windowStart)
        winWEnd = min(wEnd, windowEnd)
        winWD = winWEnd - winWStart
        if logging.getLogger().isEnabledFor(logging.DEBUG):
            logging.debug(f"app {app.name} starting at {a.start} has Work at time {wStart} for {workd} that overlaps the window [{windowStart}; {windowEnd}] for [{winWD}]")
        return [winWD]
    #Nothing is in the window
    return []

def cleanupWindow(seq):
    seqStartsWithWork = True
    p = 0
    while p < len(seq):
        if seq[p] < 1:
            # We truncate things that last less than 1s, and move them in the next phase of same type
            if p+2 < len(seq):
                seq[p+2] += seq[p]
                seq[p] = 0
            if p == 0:
                # If we are removing the first element, we need to remember that we don't start with some work
                seqStartsWithWork = not seqStartsWithWork
                del seq[0]
            elif p+1 < len(seq):
                # If we are in the middle of something, we merge both other ends
                seq[p] = seq[p-1]+seq[p+1]
                del seq[p+1]
                del seq[p-1]
            else:
                # If we are removing the last element, we try to give the time to the previous phase of the same type
                if p >= 2:
                    seq[p-2] += seq[p]
                elif p == 1:
                    # We give that negligible time to the previous phase (of wrong type) to keep the window full
                    seq[p-1] += seq[p]
                # Otherwithe the sequence has self-destructed anyway
                del seq[p]
            # The list changed, start from the beginning
            p = 0
        else:
            p += 1
    if seqStartsWithWork:
        return ('W', seq)
    else:
        return ('C', seq)

def generateWindow(plat, sched, minTime, maxTime, minDuration, outname):
    timeline = sorted( set([a.start for a in sched.applist]) )
    winIntervals = pd.DataFrame( {'windowDuration': [ timeline[i+1] - timeline[i] for i in range(len(timeline)-1)],
                                  'window': range(len(timeline)-1),
                                  'startDate': timeline[:-1]} )
    # Keep only the windows between minTime and maxTime
    winIntervals['endDate'] = winIntervals['startDate'] + winIntervals['windowDuration']
    winIntervals = winIntervals[ ((winIntervals.startDate <= maxTime) & (winIntervals.endDate >= minTime)) ]
    logging.info(f"There are {len(winIntervals)} windows between {minTime} and {maxTime}")

    if minDuration == 0:
        l = winIntervals.nlargest(20, 'windowDuration')
        llen = len(l)
        l20 = l.iloc[-1]['windowDuration']
        if outname.endswith('.png'):
            logging.info(f"Creating {outname} to store the histogram of the possible windows")
            sns.histplot(winIntervals, x='windowDuration')
            plt.savefig(outname)
            logging.info(f"Use --window {l20} to get the 20 longest windows or look at {outname} to decide what window size to use")
            print(l20)
            return None
        else:
            logging.info(f"Using window {l20} to get the 20 longest windows")
            minDuration = l20

    # Count how many there are that fit the goal
    windows = winIntervals[ winIntervals.windowDuration >= minDuration ]
    nbWindows = len(windows)
    if nbWindows == 0:
        logging.error(f"Could not find any window with duration at least {minDuration}. Try generating the histogram by passing --window 0")
        return None

    logging.info(f"Found {nbWindows} windows of {minDuration}, generating all")
    winIdx = 0
    for i, v in windows.iterrows():
        winIdx += 1
        start = timeline[v.window]
        end = timeline[v.window+1]
        logging.info(f"Window {i}: duration {v.windowDuration} starts at {start} ends at {end}")
        apps = []
        for a in sched.applist:
            if a.start <= end and a.start + a.duration >= start:
                apps.append(a)
        dotpos = outname.rfind('.')
        if dotpos > 0:
            fulloutname = outname[0:dotpos]+f"-{winIdx:03d}."+outname[dotpos+1:]
            runname = outname[0:dotpos]+f"-{winIdx:03d}.run"
        else:
            fulloutname = outname+f"-{winIdx:03d}"
            runname     = outname+f"-{winIdx:03d}.run"
        logging.info(f"Window {i} contains {len(apps)} applications -- creating {fulloutname}")
        with open(fulloutname, "w") as f:
            idx=0
            for a in apps:
                app = plat.workload[a.idx]
                nbNodes = a.nbNodes
                walltime = a.duration
                nbckpt = nbCkpt(plat, app, nbNodes, walltime)
                inputD = inputDuration(plat, app, nbNodes)
                outputD = outputDuration(plat, app, nbNodes)
                ckptD = ckptDuration(plat, app, nbNodes)
                totcompD = compDuration(plat, app, nbNodes, walltime)
                bwIO = min(nbNodes*plat.node_bw, plat.system_bw)
                t = a.start
                seq = windowEvents(t, 'input', 0, inputD, start, end, app, a)
                t += inputD
                workDone = 0
                for i in range(nbckpt):
                    seq += windowEvents(t, 'ckpt', app.ckpt_period, ckptD, start, end, app, a)
                    workDone += app.ckpt_period
                    t += app.ckpt_period
                    t += ckptD
                seq += windowEvents(t, 'output', totcompD-workDone, outputD, start, end, app, a)

                #Add some uniform noise
                for i in range(len(seq)-1):
                    v = min(seq[i], seq[i+1])
                    if v > 0:
                        n = random.gauss( 0, STDEV*v )
                        if seq[i+1]-n > 0 and seq[i]+n > 0:
                            seq[i] += n
                            seq[i+1] -= n

                #Some cleanup...
                seqDebug = seq.copy()
                (seqStartsWith, seq) = cleanupWindow(seq)

                if len(seq) == 0 or sum(seq) == 0:
                    # The algorithm will find events of duration 0 that start at the end of the window
                    continue
                if (seqStartsWith=='W') and (len(seq) == 1):
                    # Remove sequences with only some work at the beginning
                    continue

                #Convert communications back to bytes
                if seqStartsWith == 'C':
                    commOffset=0
                else:
                    commOffset=1
                for i in range(math.floor(len(seq)/2)):
                    seq[2*i+commOffset] = seq[2*i+commOffset]*bwIO

                nbPhases = len(seq)
                #Done: write the sequence with the rest of the information
                seqstr = " ".join(str(n) for n in seq)
                name = f"{app.name}{idx}"
                age = f"{start-a.start}"        # start of window - date of app start
                progress = f"{start-a.start}"   # actual amount of work done (comp time + effective I/O time) -- is the same as age for now...
                agePhase = f"0"                 # how long ago did the current phase start?
                wIter = "0"                     # wIter only for set10
                f.write(f"{name} {bwIO} {age} {progress} {agePhase} {nbPhases} {wIter} {seqStartsWith} {seqstr}\n")
                idx+=1
        logging.info(f"creating {runname}")
        with open(runname, "w") as f:
            for version in ["v1", "v2"]:
                for method in ["fairShare", "FCFS", "greedyYield", "greedyCom", "Set10Learn", "fixedWindow", "lookAheadGreedyYield", "nextEvLexMin"]:
                    f.write(f"{fulloutname} {plat.system_bw} 1 {method} {version} results/results{plat.name}/outw{winIdx}-{version}.txt\n")

def main() -> None:
    parser = init_argparse()
    args = parser.parse_args()
    
    plat = None
    if args.platform == "nersc":
        plat = nersc
    elif args.platform == "trilab":
        plat = trilab
    elif args.platform == "debug":
        plat = debug
    plat.scale(FACTORS)
    
    if args.debug:
        logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)
    else:
        logging.basicConfig(stream=sys.stderr, level=logging.INFO)

    if args.mode == "schedule":
        print("Use apex-sched.py to create a schedule")
        sys.exit(1)

    if args.mode == "draw":
        if args.output == "-":
            args.output = args.input
            suffix = args.output.rfind('.sched')
            if suffix > 0:
                args.output = args.output[:-len(".sched")] + '.png'
        with open(args.input, 'rb') if args.input != "-" else sys.stdin.buffer as f:
            [plat, sched] = pickle.load(f)
            if args.maxTime == "minTime":
                end = sched.minStartTime
            elif args.maxTime == "maxTime":
                end = sched.maxEndTime
            else:
                end = float(args.maxTime)
            if args.minTime == "minTime":
                start = 0
            elif args.minTime == "maxTime":
                start = sched.minStartTime
            else:
                start = float(args.minTime)
            if not args.output.endswith(".png"):
                pngfile = args.output+".png"
            else:
                pngfile = args.output
            logging.info(f"Generating {pngfile} between {int(start/24/3600)} days and {int(end/24/3600)} days ({args.maxTime})")
            drawRealisticSchedule(plat, sched, start, end, int(args.width), plat.nodes, pngfile)
            sys.exit(0)
        sys.exit(1)

    if args.mode == "window":
        if args.output == "-":
            args.output = args.input
            suffix = args.output.rfind('.sched')
            if suffix > 0:
                args.output = args.output[:-len(".sched")] + '.txt'
        with open(args.input, 'rb') if args.input != "-" else sys.stdin.buffer as f:
            [plat, sched] = pickle.load(f)
            if args.maxTime == "minTime":
                end = sched.minStartTime
            elif args.maxTime == "maxTime":
                end = sched.maxEndTime
            else:
                end = float(args.maxTime)
            if args.minTime == "minTime":
                start = 0
            elif args.minTime == "maxTime":
                start = sched.minStartTime
            else:
                start = float(args.minTime)
            generateWindow(plat, sched, start, end, float(args.window), args.output)
            sys.exit(0)
        sys.exit(1)

main()
