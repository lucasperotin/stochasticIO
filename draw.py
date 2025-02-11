import glob
import numpy as np
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import os

objectives=["minYield","Util","UtilGopi"]
objectives2=["MINYIELD","EFFICIENCY","UTILIZATION"]

draw=["Gopi","Gopi2"]
labels=["fairShare", "FCFS","Set10Learn", "greedyYield","greedyCom","lookAheadGreedyYield","fixedWindow","nextEvLexMin"]
labels2=["FairShare","FCFS","Set-10","GreedyYield","GreedyCom","LookAheadGreedyYield","PeriodicGreedyYield","BestNextEvent"]
#colors = ["red", "blue", "green", "purple", "orange","brown","black","yellow"]
#colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray']

# Obtention de l'instance du colormap Set1
cmap = cm.get_cmap('Set1')

# Obtention de 8 couleurs du colormap Set1
colors = cmap.colors[:8]

# Use the Seaborn default color palette
colors = [(0.12156862745098039, 0.4666666666666667, 0.7058823529411765), (1.0, 0.4980392156862745, 0.054901960784313725), (0.8392156862745098, 0.15294117647058825, 0.1568627450980392),  (0.5490196078431373, 0.33725490196078434, 0.29411764705882354), (0.17254901960784313, 0.6274509803921569, 0.17254901960784313), (0.8901960784313725, 0.4666666666666667, 0.7607843137254902), (0.5803921568627451, 0.403921568627451, 0.7411764705882353), (0.4980392156862745, 0.4980392156862745, 0.4980392156862745)][:8]


if ("Gopi" in draw):
# List of possible values for w
    ws = [[0.2, 0.5, 0.8, 0.9, 1.0, 1.1],[0,0.25,0.5,0.75,1],[0,0.25,0.5,0.75,1],[10,20,30]]
    Letters=["w","b","s","nH"]
    Letters2=["$W^{GOAL}$","$\\nu$","$\sigma$","$n_{small}$"]
    Values=[0,0,0,0]
    for x in range(4):
        for obj in objectives:
            for i in range(1,3):
                data = {name: {w: [] for w in ws[x]} for name in labels}
                if(x==0):
                    files = glob.glob("./results/resultsGopi/outw*-nH20-b0.5-s0.5-v"+(str)(i)+".txt")
                if(x==1):
                    files = glob.glob("./results/resultsGopi/outw0.8-nH20-b*-s0.5-v"+(str)(i)+".txt")
                if(x==2):
                    files = glob.glob("./results/resultsGopi/outw0.8-nH20-b0.5-s*-v"+(str)(i)+".txt")
                if(x==3):
                    files=glob.glob("./results/resultsGopi/outw0.8-nH*-b0.5-s0.5-v"+(str)(i)+".txt")
                for file in files:
                    with open(file, "r") as f:
                        lines = f.readlines()
                        for line in lines:
                            name, minYield, utilization, utilGopi = line.strip().split()
                            if (name in labels):
                                if(Letters[x]=="s"):
                                    w = float(file.split(Letters[x])[5].split("-")[0])  # Extract the value of w from the file name
                                else:
                                    w = float(file.split(Letters[x])[1].split("-")[0])  # Extract the value of w from the file name
                                    
                                if (obj=="minYield"):
                                    data[name][w].append(float(minYield))
                                elif (obj=="Util"):
                                    data[name][w].append(float(utilization))
                                else:
                                    data[name][w].append(float(utilGopi))
                                    
                fig, ax = plt.subplots()
                artists = []
                avg_values = {name: [] for name in labels}
                for z,w in enumerate(ws[x]):
                    for j, name in enumerate(labels):
                        position = z + j/(len(labels)+1)
                        box = ax.boxplot(data[name][w], positions=[position], widths = 1/(len(labels)+1), labels=[w], showfliers=False, patch_artist=True, boxprops=dict(facecolor=colors[j], edgecolor=colors[j]),whiskerprops=dict(color=colors[j]),capprops=dict(color=colors[j]), medianprops=dict(color='None'),whis=(10,90))
                        if (z==0):
                            artists.append(matplotlib.patches.Patch(color=colors[j],label=labels2[labels.index(name)]))
                        #avg_values[name].append(np.nanmean(data[name][w]))
                        
                        
                        #avg_values[name].append(np.nanmean(data[name][w]))
                        ax.plot(position, np.nanmean(data[name][w]), 'x', color='black', markersize=2)
                #         if (name=="nextEvLexMin" and x==0 and w==0.9 and i==1 and obj=="minYield"):
                #             data[name][w].sort()
                #             print(data[name][w])
    
                # for j, name in enumerate(labels):
                #     ax.plot(np.arange(len(ws[x]))+j*(1/(len(labels)+1)), avg_values[name], '-', color=colors[j],label=labels2[labels.index(name)], linewidth=0.5)
                    
                #leg=ax.legend(loc = "upper center", bbox_to_anchor=(1,0),prop={'size': 6}, ncol=8)
                #for legobj in leg.legendHandles:
                 #   legobj.set_linewidth(3.0)
                ax.set_xlabel(Letters2[x])
                ax.set_ylabel(objectives2[objectives.index(obj)])
                ax.set_xticks(np.arange(len(ws[x]))+0.5,minor=False)
                ax.set_xticklabels(ws[x],  ha='center', rotation=0)
                plt.savefig("figures/figuresGopi/"+obj+Letters[x]+"v"+(str)(i)+".pdf",bbox_inches='tight')
                plt.clf()

    
    for x in range(1,4):
        for obj in objectives:
            for i in range(1,3):
                data = {name: {w: [] for w in ws[x]} for name in labels}
                if(x==1):
                    files = glob.glob("./results/resultsGopi/outw1.1-nH20-b*-s0.5-v"+(str)(i)+".txt")
                if(x==2):
                    files = glob.glob("./results/resultsGopi/outw1.1-nH20-b0.5-s*-v"+(str)(i)+".txt")
                if(x==3):
                    files=glob.glob("./results/resultsGopi/outw1.1-nH*-b0.5-s0.5-v"+(str)(i)+".txt")
                for file in files:
                    with open(file, "r") as f:
                        lines = f.readlines()
                        for line in lines:
                            name, minYield, utilization, utilGopi = line.strip().split()
                            if (name in labels):
                                if(Letters[x]=="s"):
                                    w = float(file.split(Letters[x])[5].split("-")[0])  # Extract the value of w from the file name
                                else:
                                    w = float(file.split(Letters[x])[1].split("-")[0])  # Extract the value of w from the file name
                                    
                                if (obj=="minYield"):
                                    data[name][w].append(float(minYield))
                                elif (obj=="Util"):
                                    data[name][w].append(float(utilization))
                                else:
                                    data[name][w].append(float(utilGopi))
                                    
                fig, ax = plt.subplots()
                artists = []
                avg_values = {name: [] for name in labels}
                for z,w in enumerate(ws[x]):
                    for j, name in enumerate(labels):
                        position = z + j/(len(labels)+1)
                        box = ax.boxplot(data[name][w], positions=[position], widths = 1/(len(labels)+1), labels=[w], showfliers=False, patch_artist=True, boxprops=dict(facecolor=colors[j], edgecolor=colors[j]),whiskerprops=dict(color=colors[j]),capprops=dict(color=colors[j]), medianprops=dict(color='None'),whis=(10,90))
                        if z==0:
                            artists.append(matplotlib.patches.Patch(color=colors[j],label=labels2[labels.index(name)]))
                        #avg_values[name].append(np.nanmean(data[name][w]))
                        ax.plot(position, np.nanmean(data[name][w]), 'x', color='black', markersize=2)
                            
                        
    
                #for j, name in enumerate(labels):
                    
                 #   ax.plot(np.arange(len(ws[x]))+j*(1/(len(labels)+1)), avg_values[name], '-', color=colors[j],label=labels2[labels.index(name)], linewidth=0.5)
                    
             #   leg=ax.legend(loc = "lower left", bbox_to_anchor=(1,0),prop={'size': 10},ncol=1)
              #  for legobj in leg.legendHandles:
               #     legobj.set_linewidth(3.0)
                ax.set_xlabel(Letters2[x])
                ax.set_ylabel(objectives2[objectives.index(obj)])
                ax.set_xticks(np.arange(len(ws[x]))+0.5,minor=False)
                ax.set_xticklabels(ws[x],  ha='center', rotation=0)
                plt.savefig("figures/figuresGopi/"+obj+Letters[x]+"v"+(str)(i)+"DENSE.pdf",bbox_inches='tight')
                plt.clf()

                 

if ("Apex" in draw):
# List of possible values for w
    ws = [0.2, 0.5, 0.8, 0.9, 1, 1.1]
    Values=[0,0,0,0]
    stris=["Nersc","Trilab"]
    for stri in stris:
        if (stri=="Nersc"):
            stri2="nersc"
        else:
            stri2="trilab"
        for obj in objectives:
            for i in range(1,3):
                data = {name: {w: [] for w in ws} for name in labels}
                files = glob.glob("./results/results"+stri2+"/outw*-v"+(str)(i)+".txt")
                
                for file in files:
                    w = float(file.split("w")[1].split("-")[0])  # Extract the value of w from the file name
                    with open(file, "r") as f:
                        lines = f.readlines()
                        for line in lines:
                            #print(line)
                            name, minYield, utilization, utilGopi = line.strip().split()
                            if (name in labels):
                                if (obj=="minYield"):
                                    data[name][w].append(float(minYield))
                                elif (obj=="Util"):
                                    data[name][w].append(float(utilization))
                                else:
                                    data[name][w].append(float(utilGopi))
        
                avg_data = {name: {w: np.mean(values) for w, values in data[name].items()} for name in data}
        
                fig, ax = plt.subplots()
                artists = []
                avg_values = {name: [] for name in labels}
                for z,w in enumerate(ws):
                    for j, name in enumerate(labels):
                        position = z + j/(len(labels)+1)
                        box = ax.boxplot(data[name][w], positions=[position], widths = 1/(len(labels)+1), labels=[w], showfliers=False, patch_artist=True, boxprops=dict(facecolor=colors[j], edgecolor=colors[j]),whiskerprops=dict(color=colors[j]),capprops=dict(color=colors[j]), medianprops=dict(color='None'),whis=(10,90))
                        if z==0:
                            artists.append(matplotlib.patches.Patch(color=colors[j],label=labels2[labels.index(name)]))
                        avg_values[name].append(np.nanmean(data[name][w]))
                        
        
                for j, name in enumerate(labels):
                    ax.plot(np.arange(len(ws))+j*(1/6), avg_values[name], '-', color=colors[j],label=labels2[labels.index(name)], linewidth=0.5)
                    
                leg = plt.legend(loc='lower left', bbox_to_anchor=(0.5, 1.15), ncol=1, prop={'size': 6})
                ax.set_xlabel("W")
                ax.set_ylabel(objectives2[objectives.index(obj)])
                ax.set_xticks(np.arange(len(ws))+0.5,minor=False)
                ax.set_xticklabels(ws,  ha='center', rotation=0)
                plt.savefig("figures/figures"+stri+"/"+obj+"wv"+(str)(i)+".pdf")
                plt.clf()

if("Gopi2" in draw):
    # Define the names of the possible labels
    # Create a dictionary to store the data
    
    # Iterate over all files in the "detailed" folder
    for file_name in os.listdir("results/detailedGopi"):
        if file_name.endswith(".txt"):
            data = {name: {i: [] for i in range(60)} for name in labels}
            # Open the file
            with open(os.path.join("results/detailedGopi", file_name)) as file:
                # Iterate over each line of the file
                for line in file:
                    line = line.strip()
                    values = line.split(" ")
                    label = values[0]
                    if label in labels:
                        for i in range(60):
                            data[label][i].append(float(values[i+1]))
            # Iterate over the labels and plot the averages
            for j,name in enumerate(labels):
                averages = [sum(data[name][i])/len(data[name][i]) for i in range(60)]
                plt.plot(range(1, 61), averages, label=labels2[labels.index(name)],color=colors[j])
                
            # Add the plot labels and legend
            plt.xlabel("Index")
            plt.ylabel("Average Value")
            #plt.legend()
            #plt.legend(loc='lower left', bbox_to_anchor=(0.5, 1.15), ncol=1, prop={'size': 6})
            # for legobj in leg.legendHandles:
            #     legobj.set_linewidth(3.0)
            
        # Save the plot to a file with the same name as the data file
        plt.savefig("figures/figuresGopi/"+os.path.splitext(file_name)[0].replace('.', ',', 1) + ".pdf")
        
        # Clear the plot for the next iteration
        plt.clf()
                
if ("Apex2" in draw):
    
    # Define the names of the possible labels
    # Create a dictionary to store the data
    
    stris=["Nersc","Trilab"]
    for stri in stris:
        if (stri=="Nersc"):
            NB=60
        else:
            NB=15
    # Iterate over all files in the "detailed" folder
        for file_name in os.listdir("results/detailed"+stri):
            if file_name.endswith(".txt"):
                data = {name: {i: [] for i in range(NB)} for name in labels}
                # Open the file
                with open(os.path.join("results/detailed"+stri, file_name)) as file:
                    # Iterate over each line of the file
                    for line in file:
                        line = line.strip()
                        values = line.split(" ")
                        label = values[0]
                        for i in range(NB):
                            data[label][i].append(float(values[i+1]))
                # Iterate over the labels and plot the averages
                for j,name in enumerate(labels):
                    averages = [sum(data[name][i])/len(data[name][i]) for i in range(NB)]
                    plt.plot(range(1, NB+1), averages, label=labels2[labels.index(name)],color=colors[j])
                    
                # Add the plot labels and legend
                plt.xlabel("Index")
                plt.ylabel("Average Value")
                plt.legend()
                
                # Save the plot to a file with the same name as the data file
                plt.savefig("figures/figures"+stri+"/"+os.path.splitext(file_name)[0].replace('.', ',', 1) + ".pdf")
                
                # Clear the plot for the next iteration
                plt.clf()
