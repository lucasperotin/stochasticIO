#!env python3

import sys
import glob
import math

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
        for arg in ["nodes", "cores_per_node", "ram_per_node", "node_bw", "system_bw", "workload"]:
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
        return f"Platform(nodes={self.nodes}, cores_per_node={self.cores_per_node}, ram_per_node={self.ram_per_node}, node/system_bw={self.node_bw}/{self.system_bw}, workload={repr(self.workload)})"

debug = Platform(nodes = 10, 
                 cores_per_node = 1, 
                 ram_per_node = 1, 
                 node_bw = 1, 
                 system_bw = 2, 
                 workload = [AppClass(name = 'A', corehour_ratio = 0.5, walltime = 3600, cores = 2, input_ratio = 1, output_ratio = 1, ckpt_ratio = 0.3, iterations = 3, ensemble = 1, ckpt_period = 600),
                             AppClass(name = 'B', corehour_ratio = 0.5, walltime = 7200, cores = 4, input_ratio = 0.5, output_ratio = 0.5, ckpt_ratio = 0.3, iterations = 3, ensemble = 1, ckpt_period = 1200)])

trilab = Platform(nodes = 8894, cores_per_node = 16, ram_per_node = 32*1024*1024*1024, node_bw = 1.2*1024*1024*1024, system_bw = 200*1024*1024*1024, 
                  workload = [AppClass(name = "EAP", corehour_ratio = 0.2, walltime = 262.4*3600, cores = 16384, input_ratio = 0.03, output_ratio = 1.05, ckpt_ratio = 0.3, iterations = 30, ensemble = 1, ckpt_period = 3600),
                              AppClass(name = "LAP", corehour_ratio = 0.02, walltime = 64*3600, cores = 4096, input_ratio = 0.05, output_ratio = 3.2, ckpt_ratio = 0.75, iterations = 10, ensemble = 1, ckpt_period = 3600),
                              AppClass(name = "Silverton", corehour_ratio = 0.05, walltime = 128*3600, cores = 32768, input_ratio = 0, output_ratio = 0.43, ckpt_ratio = 2.1, iterations = 6, ensemble = 1, ckpt_period = 3600),
                              AppClass(name = "VPIC", corehour_ratio = 0.03, walltime = 157.2*3600, cores = 30000, input_ratio = 0.05, output_ratio = 2.7, ckpt_ratio = 0.1875, iterations = 4, ensemble = 1, ckpt_period = 3600),
                              AppClass(name = "DakotaA", corehour_ratio = 0.03, walltime = 100*3600, cores = 8192, input_ratio = 0.05, output_ratio = 0.085, ckpt_ratio = 0.002, iterations = 10, ensemble = 100, ckpt_period = 3600,),
                              AppClass(name = "DakotaS", corehour_ratio = 0.03, walltime = 100*3600, cores = 4096, input_ratio = 0.005, output_ratio = 0.0204, ckpt_ratio = 0.3, iterations = 30, ensemble = 300, ckpt_period = 3600)])

nersc = Platform(nodes = 8894, cores_per_node = 16, ram_per_node = 32*1024*1024*1024, node_bw = 1.2*1024*1024*1024, system_bw = 200*1024*1024*1024,
    workload = [AppClass(name = "ALS", corehour_ratio = 0.001, walltime = 8*0.5*3600, cores = 100, input_ratio = 0.3165, output_ratio = 1.0629, ckpt_ratio = 0, iterations = 10760/8, ensemble = 1, ckpt_period = 3600),
                AppClass(name = "CESM", corehour_ratio = 0.06, walltime = 240*3600, cores = 8000, input_ratio = 0, output_ratio = 10.2795, ckpt_ratio = 0.32, iterations = 8, ensemble = 1, ckpt_period = 3600),
                AppClass(name = "GTS", corehour_ratio = 0.06, walltime = 48*3600, cores = 16512, input_ratio = 0, output_ratio = 0.1478, ckpt_ratio = 0.67, iterations = 36, ensemble = 1, ckpt_period = 3600),
                AppClass(name = "HipMer", corehour_ratio = 0.001, walltime = 4*3600, cores = 960, input_ratio = 0.6583, output_ratio = 0.3434, ckpt_ratio = 0, iterations = 100, ensemble = 1, ckpt_period = 3600),
                AppClass(name = "Materials", corehour_ratio = 0.19, walltime = 41.7*3600, cores = 2400, input_ratio = 0, output_ratio = .2708, ckpt_ratio = 0.5, iterations = 100, ensemble = 1, ckpt_period = 3600),
                AppClass(name = "MILC", corehour_ratio = 0.11, walltime = 6*3600, cores = 4096, input_ratio = 0.7388, output_ratio = 0.7388, ckpt_ratio = 0, iterations = 1000, ensemble = 1, ckpt_period = 3600),
                AppClass(name = "SkySurvey", corehour_ratio = 0.06, walltime = 2*4*3600, cores = 16*24, input_ratio = 0.0879, output_ratio = 0.0269, ckpt_ratio = 0, iterations = 21000/16/2, ensemble = 1, ckpt_period = 3600)])

PLATFORMS = {
    'debug': debug,
    'nersc': nersc,
    'trilab': trilab
}

def windowStats(filename):
    dirs = filename.split('/')
    bnum = 0
    elts = dirs[-1].split('-')
    platname = elts[0]
    if len(dirs) >= 2:
        elts2 = dirs[-2].split('-')
        if len(elts2) == 2:
            bnum = elts2[1]
            platname = elts2[0]
    plat = PLATFORMS[platname]
    fnode = float(elts[1])
    fcore = float(elts[2])
    fram  = float(elts[3])
    fnodebw = float(elts[4])
    fsysbw  = float(elts[5])
    win = []
    tot_nodes = 0
    with open(filename, 'r') as f:
        for line in f.readlines():
            a = line.split(' ')
            name = a[0]
            bw = float(a[1])
            n = int(a[5])
            comm_next = (a[7] == 'C')
            phases = []
            tot_dur = 0
            nodes = int(math.ceil(bw/(plat.node_bw*fnodebw)))
            tot_nodes += nodes
            for i in range(8, 8+n):
                if comm_next:
                    d = float(a[i])/bw
                else:
                    d = float(a[i])
                tot_dur += d
                phases.append( {'comm': comm_next, 'value': a[i], 'dur': d} )
                comm_next = not comm_next
            win.append({'name': name, 'bw': bw, 'nphases': n, 'comm_first': phases[0]['comm'], 'phases': phases, 'tot_dur': tot_dur, 'nodes': nodes})
        min_dur = win[0]['tot_dur']
        max_dur = win[0]['tot_dur']
        sum_dur = win[0]['tot_dur']
        for i in range(1, len(win)):
            if win[i]['tot_dur'] < min_dur:
                min_dur = win[i]['tot_dur']
            if win[i]['tot_dur'] > max_dur:
                max_dur = win[i]['tot_dur']
            sum_dur += win[i]['tot_dur']
        print(f"{platname},{fnode},{len(win)},{tot_nodes},{tot_nodes/(fnode*plat.nodes)},{sum_dur/len(win)}")

print("plat,fnode,nbApps,nbNodes,usage,windur")
for g in sys.argv[1:]:
    for f in glob.glob(g):
        windowStats(f)
