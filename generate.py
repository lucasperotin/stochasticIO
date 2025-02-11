import random
import math

def generate(nH, h, w, b, sigma, filename):
  # Generate the table AIO
  AIO = [random.uniform(0, 1) for _ in range(60)]
  # Calculate AIOtemp
  AIOtemp = sum(AIO)
  # Scale the values in AIO by w/AIOtemp
  AIO = [a * (w / AIOtemp) for a in AIO]

  # Open the file in write mode
  with open(filename, 'w') as f:
    # Write the first nH lines
    for i in range(nH):
      # Generate wIter from normal distribution
      wIter=0
      while(wIter<=2):
          wIter = random.normalvariate(1000, 1000*sigma)
      # Calculate n and write it to the file
      n = math.ceil(h/wIter)
      m=2*n+1
      f.write("App"+(str)(i)+" 1 0 0 0 "+str(m)+" "+str(wIter)+" W ")
      # Generate the next element and write it to the file
      f.write(str(random.uniform(0, wIter))+" ")
      TCPU=wIter*(1-AIO[i]) 
      TIO= wIter*AIO[i]
      for j in range(n):
          x=random.uniform(-b,b)
          temp=(1+x)*TIO
          f.write(str(temp)+" ")
          x=random.uniform(-b,b)
          temp=(1+x)*TCPU
          f.write(str(temp)+" ")
      f.write('\n')

    # Write the next 20 lines
    for i in range(nH, 20+nH):
      # Generate wIter from normal distribution
      wIter=0
      while(wIter<=2):
          wIter = random.normalvariate(10000, 10000*sigma)
      # Calculate n and write it to the file
      n = math.ceil(h/wIter)
      m=2*n+1
      f.write("App"+(str)(i)+" 1 0 0 0 "+str(m)+" "+str(wIter)+" W ")
      # Generate the next element and write it to the file
      f.write(str(random.uniform(0, wIter))+" ")
      TCPU=wIter*(1-AIO[i]) 
      TIO= wIter*AIO[i]
      for j in range(n):
          x=random.uniform(-b,b)
          temp=(1+x)*TIO
          f.write(str(temp)+" ")
          x=random.uniform(-b,b)
          temp=(1+x)*TCPU
          f.write(str(temp)+" ")
      f.write('\n')

    # Write the remaining lines
    for i in range(20+nH, 60):
      # Generate wIter from normal distribution
      wIter=0
      while(wIter<=2):
          wIter = random.normalvariate(100000, 100000*sigma)
      # Calculate n and write it to the file
      n = math.ceil(h/wIter)
      m=2*n+1
      f.write("App"+(str)(i)+" 1 0 0 0 "+str(m)+" "+str(wIter)+" W ")
      # Generate the next element and write it to the file
      f.write(str(random.uniform(0, wIter))+" ")
      TCPU=wIter*(1-AIO[i]) 
      TIO= wIter*AIO[i]
      for j in range(n):
          x=random.uniform(-b,b)
          temp=(1+x)*TIO
          f.write(str(temp)+" ")
          x=random.uniform(-b,b)
          temp=(1+x)*TCPU
          f.write(str(temp)+" ")
      f.write('\n')


folder="files/gopi/"
sigmas=[0,0.25,0.75,1]
bs=[0,0.25,0.75,1]
ws=[0.2,0.5,0.8,0.9,1.0,1.1]
h=2000000
nHs=[10,30]

b=0.5
sigma=0.5
w1=0.8
w2=1.1
nH=20

nbSamples=10

nH=20
if(True):
    for sigma in sigmas:
        for i in range(nbSamples):
            generate(nH,h,w1,b,sigma,folder+"w"+str(w1)+"-nH"+str(nH)+"-b"+str(b)+"-s"+str(sigma)+"_"+str(i)+".txt")
            generate(nH,h,w2,b,sigma,folder+"w"+str(w2)+"-nH"+str(nH)+"-b"+str(b)+"-s"+str(sigma)+"_"+str(i)+".txt")      
    sigma=0.5
    
    
    for b in bs:
        for i in range(nbSamples):
            generate(nH,h,w1,b,sigma,folder+"w"+str(w1)+"-nH"+str(nH)+"-b"+str(b)+"-s"+str(sigma)+"_"+str(i)+".txt") 
            generate(nH,h,w2,b,sigma,folder+"w"+str(w2)+"-nH"+str(nH)+"-b"+str(b)+"-s"+str(sigma)+"_"+str(i)+".txt")       
    b=0.5 
    
    for w in ws:
        for i in range(nbSamples):
            generate(nH,h,w,b,sigma,folder+"w"+str(w)+"-nH"+str(nH)+"-b"+str(b)+"-s"+str(sigma)+"_"+str(i)+".txt")      
    w=1.1
               
    for nH in nHs:
        for i in range(nbSamples):
            generate(nH,h,w1,b,sigma,folder+"w"+str(w1)+"-nH"+str(nH)+"-b"+str(b)+"-s"+str(sigma)+"_"+str(i)+".txt") 
            generate(nH,h,w2,b,sigma,folder+"w"+str(w2)+"-nH"+str(nH)+"-b"+str(b)+"-s"+str(sigma)+"_"+str(i)+".txt")       
    nH=20
                    

heuristics=["fairShare", "FCFS", "greedyYield", "greedyCom", "Set10Learn", "fixedWindow","lookAheadGreedyYield","nextEvLexMin"]
with open("args.txt", "w") as out_file:
    for sigma in sigmas:
        for i in range(nbSamples):
            filename = "w{}-nH{}-b{}-s{}_{}.txt".format(w1, nH, b, sigma, i)
            outname="w{}-nH{}-b{}-s{}".format(w1, nH, b, sigma)
            for j in heuristics:
                out_file.write("./"+folder+"{} 1 1 {} v1 results/resultsGopi/out{}-v1.txt\n".format(filename, j ,outname))
                out_file.write("./"+folder+"{} 1 1 {} v2 results/resultsGopi/out{}-v2.txt\n".format(filename, j, outname))
            filename = "w{}-nH{}-b{}-s{}_{}.txt".format(w2, nH, b, sigma, i)
            outname="w{}-nH{}-b{}-s{}".format(w2, nH, b, sigma)
            for j in heuristics:
                out_file.write("./"+folder+"{} 1 1 {} v1 results/resultsGopi/out{}-v1.txt\n".format(filename, j ,outname))
                out_file.write("./"+folder+"{} 1 1 {} v2 results/resultsGopi/out{}-v2.txt\n".format(filename, j, outname))
    sigma=0.5
    
    for b in bs:
        for i in range(nbSamples):
            filename = "w{}-nH{}-b{}-s{}_{}.txt".format(w1, nH, b, sigma, i)
            outname="w{}-nH{}-b{}-s{}".format(w1, nH, b, sigma)
            for j in heuristics:
                out_file.write("./"+folder+"{} 1 1 {} v1 results/resultsGopi/out{}-v1.txt\n".format(filename, j ,outname))
                out_file.write("./"+folder+"{} 1 1 {} v2 results/resultsGopi/out{}-v2.txt\n".format(filename, j, outname))
            filename = "w{}-nH{}-b{}-s{}_{}.txt".format(w2, nH, b, sigma, i)
            outname="w{}-nH{}-b{}-s{}".format(w2, nH, b, sigma)
            for j in heuristics:
                out_file.write("./"+folder+"{} 1 1 {} v1 results/resultsGopi/out{}-v1.txt\n".format(filename, j ,outname))
                out_file.write("./"+folder+"{} 1 1 {} v2 results/resultsGopi/out{}-v2.txt\n".format(filename, j, outname))
    b=0.5 
    
    for w in ws:
        for i in range(nbSamples):
            filename = "w{}-nH{}-b{}-s{}_{}.txt".format(w, nH, b, sigma, i)
            outname="w{}-nH{}-b{}-s{}".format(w, nH, b, sigma)
            for j in heuristics:
                out_file.write("./"+folder+"{} 1 1 {} v1 results/resultsGopi/out{}-v1.txt\n".format(filename, j ,outname))
                out_file.write("./"+folder+"{} 1 1 {} v2 results/resultsGopi/out{}-v2.txt\n".format(filename, j, outname))
 
    w=0.9
                
    for nH in nHs:
        for i in range(nbSamples):
            
            filename = "w{}-nH{}-b{}-s{}_{}.txt".format(w1, nH, b, sigma, i)
            outname="w{}-nH{}-b{}-s{}".format(w1, nH, b, sigma)
            for j in heuristics:
                out_file.write("./"+folder+"{} 1 1 {} v1 results/resultsGopi/out{}-v1.txt\n".format(filename, j ,outname))
                out_file.write("./"+folder+"{} 1 1 {} v2 results/resultsGopi/out{}-v2.txt\n".format(filename, j, outname))
            
            filename = "w{}-nH{}-b{}-s{}_{}.txt".format(w2, nH, b, sigma, i)
            outname="w{}-nH{}-b{}-s{}".format(w2, nH, b, sigma)
            for j in heuristics:
                out_file.write("./"+folder+"{} 1 1 {} v1 results/resultsGopi/out{}-v1.txt\n".format(filename, j ,outname))
                out_file.write("./"+folder+"{} 1 1 {} v2 results/resultsGopi/out{}-v2.txt\n".format(filename, j, outname))
 
    nH=20
