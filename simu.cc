#include <iostream>
#include <iomanip>
#include <vector>
#include <queue>
#include <string>
#include <fstream>
#include <utility>
#include <algorithm>
#include <math.h>
#include <cassert>
#include <assert.h>
#include <cmath>
#include <sstream>

#ifdef _OPENMP
   #include <omp.h>
#else
   #define omp_get_num_threads() 0
#endif

using namespace std;

#define ERROR 0.0000000001
#define ERROR2 0.000000001
#define INV_PRECISION 10000
#define SINGULARITY_PRECISION 0.001
#define PRINTCOMPUTENEXTEVENT 0
#define DEBUGNELEXMIN 0
#define ALLDEBUGNELEXMIN 0
#define EXTRAALLDEBUGNELEXMIN 0
#define HACKSHORTCOMM 1
#define CATCHERROR false
#define AUTOSCALE false
#define LIMITRECURSIVEDEPTH 10
#define CHECKOPTIMISATION 0

//**************GLOBAL VARIABLES AND CLASSES******

class Application
{
public:
    Application(string n, double r, double i,vector<double> v, double k,double rW,double rR) : name(n), release(r), iniProgress(i), Volumes(v), knownWiter(k), realWiter(rW), realRatio(rR) {}
    void setTProgress(double p)  // Set & Get the time required by the application to progress as much as it would do if it was running solo.
    {
        totProgress=p;
    }
    double getTProgress()
    {
        return totProgress;
    }
    double getYield(double t) /// Convention, yield is 0 if the application is new. You may not set the yield
    {
        if(t>release)
        {
            return (totProgress/(t-release));
        }
        else
        {
            return 0;
        }
    }
    void setRemVol(double p) // Set & get the remaining volume of the current phase
    {
        curRemVol=p;
    }
    double getRemVol()
    {
        return curRemVol;
    }
    void setPhaseId(long v) // Set & get the Id of the current Phase (likely useless for you)?
    {
        phaseId=v;
    }
    long getPhaseId()
    {
        return phaseId;
    }
    double getCurAct()  // Get the total volume of the current phase (might be useless?).
    {
        return Volumes[phaseId];
    }
    void setType(bool b) // Set & get the type of the current type of phase : Communication (true) or Work (false)
    {
        curType=b;
    }
    bool getType()
    {
        return curType;
    }
    void setBW(double b) // Set & get the current PassBand allocated to the application
    {
        allocatedBandwidth=b;
    }
    double getBW()
    {
        return allocatedBandwidth;
    }
    void setReleasePhase(double t) // Set & get the current PassBand allocated to the application
    {
        releasePhase=t;
    }
    double getReleasePhase()
    {
        return releasePhase;
    }
    void setWiter(double w)
    {
        wIter=w;
    }
    double getWiter(){
        return wIter;
    }
    void setMaxBandwidth(double w)
    {
        maxBandwidth=w;
    }
    double getMaxBandwidth(){
        return maxBandwidth;
    }
    void setCountIter(long w)
    {
        countIter=w;
    }
    double getCountIter(){
        return countIter;
    }
    void setUtilGopi(long w)
    {
        UtilGopi=w;
    }
    double getUtilGopi(){
        return UtilGopi;
    }
    void setNew(bool w)
    {
        isNew=w;
    }
    bool getNew(){
        return isNew;
    }
    // These may not be modified and are constant
    string name;
    double release;         // Release date
    double iniProgress;     // Progress at the beginning of the window
    vector<double> Volumes; // All volumes
    bool knownWiter;
    double realWiter;
    double realRatio;

private:
    double maxBandwidth;        // Maximum bandwidth that may be allocated
    double totProgress;
    double curRemVol;
    long phaseId;
    double allocatedBandwidth;
    double releasePhase;
    bool curType;
    double wIter;
    long countIter;
    double UtilGopi;
    bool isNew;
};

//**************PRINT AND DEBUG FUNCTIONS ****************

template <int debuga>
void printApplication(unsigned long ind, double cTime, std::vector<Application>& apps)
{
    unsigned long i;
    vector<double> vols=apps[ind].Volumes;
    string sttype;
    if (apps[ind].getType()==true)
    {
        sttype="Com (true)";
    }
    else
    {

    	sttype="Work (false)";
    }
    std::cout<<"Print Application " << apps[ind].name << "\n\n";
    cout << "Release date : "<< apps[ind].release << ", (max) bandwidth : " << apps[ind].getMaxBandwidth() << "\n";
    if constexpr (debuga>1){
    	cout << "Volumes : ";
    	for (i=0; i<vols.size(); i++)
    	{
    		cout << vols[i] << ' ';
    	}
    }
    cout<<"\n\n";
    cout << "knownWiter ? "<< apps[ind].knownWiter << " wIter " << apps[ind].getWiter()<<"\n";
    cout << "Total progress: "<< apps[ind].getTProgress() << ", (in window : " << apps[ind].getTProgress()-apps[ind].iniProgress << "), yield : "<<apps[ind].getYield(cTime)<<"\n";
    cout << "Current allocated bandwidth: " << apps[ind].getBW() << "\n";
    cout << "Current activity: Index " << apps[ind].getPhaseId() << ", Volume " << apps[ind].getCurAct()  << ", Type "<< sttype << ", out of which "<< apps[ind].getRemVol() << " remains.\n\n\n";
}

template <int debuga>
void printApplications(double cTime, std::vector<Application> &apps)
{
    unsigned long i;
    for (i=0; i<apps.size(); i++)
    {
        if constexpr (debuga>2) {
            printApplication<debuga>(i,cTime,apps);
        } else {
            if (apps[i].getType()){
                printApplication<debuga>(i,cTime,apps);
            }
        }
    }
}

template <int debuga>
void printIni(double cTime, std::vector<Application> &apps)
{
    cout<<"At time " << cTime << ", The state of applications is as follows : \n\n";
    printApplications<debuga>(cTime,apps);
}

void printStepA(vector<unsigned long> appInCom, std::vector<Application> &apps)
{
    unsigned long m=appInCom.size();
    unsigned long ind;
    unsigned long i;
    cout<<"\n### STEP A : bandwidth decision ###\n\n";

    cout<<"The Heuristic made the following choices :\n\n";
    for(i=0; i<m; i++)
    {
        ind=appInCom[i];
        cout << "Application "<<apps[ind].name << " gets allocated bandwidth "<<apps[ind].getBW()<<" out of " << apps[ind].getMaxBandwidth() << "\n";
    }
}

//************** INPUT MANAGEMENT ****************

void addApps2(ifstream& appStr, std::vector<Application> &apps, double &endWin,double &totBand, string &heuristic)
{
    unsigned long i;

    string nam;
    double r;
    double inip;
    double bw;
    unsigned long nAct;
    vector<double> Volumes;
    double volume;
    double fiVol;
    double progressMax;
    double relPhase;
    double wIter;
    bool knownWiter=false;
    int type;
    char initype;
    bool found=false;
    double progressMax2=0;
    double progressTot=0;
    double countVol=0;
    double countApp=0;
    
    progressMax=0;
    double progressCom=0;
    double progressWork=0;
    double tempprogressCom=0;
    double tempprogressWork=0;
    while(appStr>>nam)
    {
        countApp+=1;
        found=true;
        progressMax=0;
        Volumes= {};
        appStr >> bw >> r >> inip >> relPhase >> nAct  >> wIter >> initype;
        if (wIter>0 and heuristic!="Set10Learn"){
            knownWiter=true;
        }
        else{
            knownWiter=false;
        }
        if (initype=='C')
        {
            type=1;
        }
        else
        {
            type=0;
        }
        appStr >> fiVol;
        countVol++;
        if(type){
            fiVol=fiVol/totBand;
        }
        Volumes.push_back(fiVol);
        if(bw>totBand){
            bw=totBand;
        }
        bw=bw/totBand;
        tempprogressCom=0;
        tempprogressWork=0;
        progressMax+=fiVol*type/bw+fiVol*(1-type);
        //cout << progressMax << "\n";
        tempprogressCom+=fiVol*type;
        tempprogressWork+=fiVol*(1-type);
        for (i=1; i<nAct; i++)
        {
            type=1-type;
            appStr >> volume;
            countVol++;
            if(type){
                volume=volume/totBand;
            }
            progressMax+=volume*type/bw+volume*(1-type);
            
            //cout << progressMax << "\n";
            tempprogressCom+=volume*type;
            tempprogressWork+=volume*(1-type);
            Volumes.push_back(volume);
        }
        //cout<<fiVol<<"aa\n";
        Application cApp=Application(nam,-r,inip,Volumes,knownWiter,(tempprogressCom+tempprogressWork)/nAct*2,tempprogressCom/tempprogressWork);
        cApp.setMaxBandwidth(bw);
        cApp.setTProgress(inip);
        cApp.setRemVol(fiVol);
        cApp.setPhaseId(0);
        cApp.setBW(0);
        cApp.setType(initype=='C');
        cApp.setReleasePhase(-relPhase);
        cApp.setWiter(wIter);
        cApp.setCountIter(0);
        cApp.setUtilGopi(0);
        cApp.setNew(true);
        progressTot+=progressMax;
        progressCom+=tempprogressCom;
        progressWork+=tempprogressWork;
        if(progressMax<progressMax2 or progressMax2==0){
            progressMax2=progressMax;
        }
        //cout << nam << " " << progressMax << "\n";
    
    }
    if (not found){
        throw std::runtime_error("FATAL ERROR : File not found");
    }
    endWin=endWin*progressMax2;
    progressWork/=countApp;
    //cout << progressCom << "\n";
    //cout << progressWork<<"\n";
    //cout << "Saturation : " << progressCom/progressWork<<"\n";
}

pair<double,double> addApps(ifstream& appStr, std::vector<Application> &apps, double &endWin,double &totBand, string &heuristic)
{
    unsigned long i;

    string nam;
    double r;
    double inip;
    double bw;
    unsigned long nAct;
    vector<double> Volumes;
    double volume;
    double fiVol;
    double progressMax;
    double relPhase;
    double wIter;
    bool knownWiter=false;
    int type;
    char initype;
    bool found=false;
    double progressMax2=0;
    double progressTot=0;
    double countVol=0;
    double countApp=0;
    bool itCounts;
    
    progressMax=0;
    double progressCom=0;
    double progressWork=0;
    double tempprogressCom=0;
    double tempprogressWork=0;
    while(appStr>>nam)
    {
        countApp+=1;
        found=true;
        progressMax=0;
        Volumes= {};
        appStr >> bw >> r >> inip >> relPhase >> nAct  >> wIter >> initype;
        if (wIter>0 and heuristic!="Set10Learn"){
            knownWiter=true;
        }
        else{
            knownWiter=false;
        }
        if (initype=='C')
        {
            type=1;
        }
        else
        {
            type=0;
        }
        appStr >> fiVol;
        countVol++;
        if(type){
            fiVol=fiVol/totBand;
        }
        if(bw>totBand){
            bw=totBand;
        }
        bw=bw/totBand;
        tempprogressCom=0;
        tempprogressWork=0;
        
        itCounts=true;
        if(fiVol*type/bw+fiVol*(1-type)>endWin){
            itCounts=false;
            if (not type){
                fiVol=endWin;
            }
            else{
                fiVol=bw*endWin;
            }
        }
        assert(fiVol>0);
        Volumes.push_back(fiVol);
        progressMax+=fiVol*type/bw+fiVol*(1-type);
        tempprogressCom+=fiVol*type;
        tempprogressWork+=fiVol*(1-type);
        
        for (i=1; i<nAct; i++)
        {
            type=1-type;
            appStr >> volume;
            if(type){
                volume=volume/totBand;
            }
            if (progressMax+volume*type/bw+volume*(1-type)>endWin and itCounts){
                if (not type){
                    volume=endWin-progressMax;
                }
                else{
                    volume=bw*(endWin-progressMax);
                }
                progressMax+=volume*type/bw+volume*(1-type);
                tempprogressCom+=volume*type;
                tempprogressWork+=volume*(1-type);
                assert(volume>0);
                Volumes.push_back(volume);
                itCounts=false;
            }
            if(itCounts){
                countVol++;
                progressMax+=volume*type/bw+volume*(1-type);
                tempprogressCom+=volume*type;
                tempprogressWork+=volume*(1-type);
                assert(volume>0);
                Volumes.push_back(volume);
            }
        }
        assert(progressMax>endWin*0.999);
        assert(progressMax<endWin*1.001);
        assert(tempprogressCom>=0);        //cout<<fiVol<<"aa\n";
        Application cApp=Application(nam,-r,inip,Volumes,knownWiter,(tempprogressCom+tempprogressWork)/nAct*2,tempprogressCom/tempprogressWork);
        cApp.setMaxBandwidth(bw);
        cApp.setTProgress(inip);
        cApp.setRemVol(fiVol);
        cApp.setPhaseId(0);
        cApp.setBW(0);
        cApp.setType(initype=='C');
        cApp.setReleasePhase(-relPhase);
        cApp.setWiter(wIter);
        cApp.setCountIter(0);
        cApp.setUtilGopi(0);
        cApp.setNew(true);
        apps.push_back(cApp);
        progressTot+=progressMax;
        progressCom+=tempprogressCom;
        progressWork+=tempprogressWork;
        if(progressMax<progressMax2 or progressMax2==0){
            progressMax2=progressMax;
        }
    
    }
    if (not found){
        throw std::runtime_error("FATAL ERROR : File not found");
    }
    progressWork/=countApp;
    //cout << progressCom << "\n";
    //cout << progressWork<<"\n";
    //cout << "Saturation : " << progressCom/progressWork<<"\n";
    return (make_pair(progressTot/countVol/countApp,progressCom/endWin));
}

void autoScale(std::vector<Application> &apps, double &scaleFactor){
    
    for (unsigned int i=0;i<apps.size();i++){
        apps[i].setMaxBandwidth(apps[i].getMaxBandwidth()/scaleFactor);
    }
}
//**************WRITE HEURISTICS HERE****************

// Computes the yield of application ``appli'' at date cTime+delta if it is allocated the bandwidth allocatedBW
long double fYield(unsigned long appli, double delta, double cTime, double allocatedBW, double * releaseDate, double * progress, double * maxBandwidth,
		bool withAssert){

	if (allocatedBW > maxBandwidth[appli]+ERROR){
		cout << "application: " << appli << endl;
		cout << "allocatedBW = " << allocatedBW << endl;
		cout << "maxBandwidth= " << maxBandwidth[appli] << endl;
	}

	// We check that the assigned bandwidth is achievable
	if (withAssert) assert(allocatedBW <= maxBandwidth[appli]+ERROR);

	if ((delta == 0) && (cTime <= releaseDate[appli] + ERROR)) return 0;
	long double yield = (progress[appli]+allocatedBW*delta/maxBandwidth[appli])/(cTime+delta-releaseDate[appli]);
//	cout << "\t\t**** " << appli << " " << delta << " " << allocatedBW << " " << yield << endl;
	return yield;
}
long double fYield(unsigned long appli, double delta, double cTime, double allocatedBW, double * releaseDate, double * progress, double * maxBandwidth){
	return fYield(appli, delta, cTime, allocatedBW, releaseDate, progress, maxBandwidth, true);
}

// Computes the yield of application ``appli'' if it completes its communication at date cTime+delta
long double fYieldAtCompletion(unsigned long appli, double delta, double cTime, double * remCommVol, double * releaseDate, double * progress, double * maxBandwidth){
	assert(delta > 0);
	// The application should be able to complete its communication at time delta
	assert(delta >= remCommVol[appli]/maxBandwidth[appli]-ERROR);
	long double yield = (progress[appli]+remCommVol[appli]/maxBandwidth[appli])/(cTime+delta-releaseDate[appli]);
	return yield;
}

// Computes the maximum yield achievable at date cTime+delta by the application ``appli''
long double fYieldMax(unsigned long appli, double delta, double cTime, double * remCommVol, double * releaseDate, double * progress, double * maxBandwidth){
	long double yield = 0;

	if (EXTRAALLDEBUGNELEXMIN) cout << "\t\t\t\tfYieldMax: " << delta << " ?<=? " <<  remCommVol[appli] << "/" << maxBandwidth[appli] << " = " << remCommVol[appli]/maxBandwidth[appli] << endl;

	if (delta <= remCommVol[appli]/maxBandwidth[appli])
	{
		yield = fYield(appli, delta, cTime, maxBandwidth[appli], releaseDate, progress, maxBandwidth);
		if (ALLDEBUGNELEXMIN) cout << "\t\t\t\t\t" << appli << " fYield= " << yield << endl;
	}
	else
	{
		yield = fYieldAtCompletion(appli, delta, cTime, remCommVol, releaseDate, progress, maxBandwidth);
		if (ALLDEBUGNELEXMIN) cout << "\t\t\t\t\t"<< appli << " fYieldAtCompletion= " << yield << endl;
	}

	return yield;
}

// Check whether the new solution is better than the current best one
void keepBestSolution(unsigned long nApps,
		double * tempYields, double * tempBandwidths,
		double * bestTempYields, double * bestTempBandwidths){


//	cout << "*** keepBestSolution" << endl;
//	for (unsigned long i=0; i<nApps; i++)
//	{
//		cout << " i: " << bestTempYields[i] << " vs. " << tempYields[i] << endl;
//	}

//	// We sort the array tempYields
//	sort(tempYields, tempYields + nApps);

	bool solutionOrdered = false;
	bool foundBetterSolution = false;
	unsigned long j = 0;
	while((!solutionOrdered) && (j < nApps))
	{
		// If the yields at the considered rank are equal we just move on
		// Therefore, we check this is not the case

		// If the yield of the new solution at the current rank is greater than the current best
		// the new solution is better and we are done
		if (tempYields[j] > bestTempYields[j]+ERROR)
		{
			solutionOrdered = true;
			foundBetterSolution = true;
		}
		else
		{
			// If the yield of the new solution at the current rank is smaller than the current best
			// the new solution is worse and we are done
			if(tempYields[j] < bestTempYields[j]-ERROR)
			{
				solutionOrdered = true;
				foundBetterSolution = false;
			}
		}

		j++;
	}

	long double usedBW = 0;

    if (CHECKOPTIMISATION){
        if (foundBetterSolution){
        	cout << "\tSolution was updated from: "<< endl << "\t\t";
        }
        else{
        	cout << "\tSolution was UNCHANGED: "<< endl << "\t\t";
        }
    	for (unsigned long i=0; i<nApps; i++)
    	{
    		cout << bestTempYields[i] << " ";
    	}
    	cout << endl << "\t  vs. " << endl << "\t\t";
    	for (unsigned long i=0; i<nApps; i++)
    	{
    		cout << tempYields[i] << " ";
    	}
    	cout << endl;
    }


    if (foundBetterSolution){
//        if (CHECKOPTIMISATION){
//        	cout << "\tSolution was updated from: "<< endl;
//        	for (unsigned long i=0; i<nApps; i++)
//        	{
//        		cout << bestTempYields[i] << " ";
//        	}
//        	cout << endl << " to: " << endl;
//        	for (unsigned long i=0; i<nApps; i++)
//        	{
//        		cout << tempYields[i] << " ";
//        	}
//        	cout << endl;
//        }


//    	cout << "BEST BW allocation: ";
    	for (unsigned long i=0; i<nApps; i++)
    	{
    		bestTempBandwidths[i] = tempBandwidths[i];
    		bestTempYields[i] = tempYields[i];

    		usedBW += bestTempBandwidths[i];

//    		cout << bestTempBandwidths[i] << " ";
    	}
//    	cout << endl;
        if (ALLDEBUGNELEXMIN) cout << "\tkeepBestSolution: new solution: BW = " << usedBW << endl;
    }

}


// Compute the bandwidth required by application theApp to achieve the yield targetYield at date cTime+delta
long double requiredBandwidth(int theApp, long double targetYield, double delta, double cTime,
		double * remCommVol, double * releaseDate, double * progress, double * maxBandwidth)
{

//	cout << "\t\t\t\t\trequiredBandwidth for " << theApp << " BW=0 y= " << fYield(theApp, delta, cTime, 0, releaseDate, progress, maxBandwidth);
//	cout << "\tBW=max= " << maxBandwidth[theApp] << " y= " << fYield(theApp, delta, cTime, maxBandwidth[theApp], releaseDate, progress, maxBandwidth) << endl;
	long double BW = 0;
	if(targetYield - fYield(theApp, delta, cTime, 0, releaseDate, progress, maxBandwidth) > ERROR)
	{
		BW = (maxBandwidth[theApp]/delta)*(targetYield*(cTime+delta-releaseDate[theApp])-progress[theApp]);
		assert(BW>=0);
	}

//	cout << "\t\t\t\t\t\tBW= " << BW << endl;


	return BW;
}



// Compute the total bandwidth required to achieve the yield targetYield at date cTime+delta
long double totalRequiredBandwidth(unsigned long nApps, long double targetYield, double delta, double cTime,
		double * remCommVol, double * releaseDate, double * progress, double * maxBandwidth)
{

	if (EXTRAALLDEBUGNELEXMIN){
		cout << "\t\t\t\ttotalRequiredBandwidth for a yield of  " << targetYield << endl;
	}
	long double totalBW = 0;
	for(unsigned long i = 0; i<nApps; i++)
	{
		long double localBW = requiredBandwidth(i,  targetYield,  delta,  cTime, remCommVol, releaseDate, progress, maxBandwidth);
		totalBW += localBW;
		if (EXTRAALLDEBUGNELEXMIN){
			cout << "\t\t\t\tAppli " << i << " BW: " << localBW << endl;
		}
	}
	if (ALLDEBUGNELEXMIN){
		cout << "\t\t\t\ttotalRequiredBandwidth for a yield of  " << targetYield << " : " << totalBW << endl;
	}


	return totalBW;
}

// Compute the total bandwidth required to achieve the yield targetYield at date cTime+delta
// by all applications with the exception of constrainingApp
long double totalRequiredBandwidth(unsigned long nApps, long double targetYield, double delta, double cTime,
		double * remCommVol, double * releaseDate, double * progress, double * maxBandwidth,
		int constrainingApp)
{

	if (ALLDEBUGNELEXMIN){
		cout << "\t\t\t\ttotalRequiredBandwidth for a yield of  " << targetYield << endl;
	}
	long double totalBW = 0;
	for(unsigned long i = 0; i<nApps; i++)
	{
		if(i!=constrainingApp){
			long double localBW = requiredBandwidth(i,  targetYield,  delta,  cTime, remCommVol, releaseDate, progress, maxBandwidth);
			totalBW += localBW;
			if (ALLDEBUGNELEXMIN){
				cout << "\t\t\t\tAppli " << i << " BW: " << localBW << endl;
			}
		}
	}

	return totalBW;
}


// Compute the total bandwidth required for the free applications to achieve the yield targetYield at date cTime+delta
long double totalRequiredBandwidth(unsigned long nApps,
		bool * freeApps,
		double targetYield, double delta, double cTime,
		double * remCommVol, double * releaseDate, double * progress, double * maxBandwidth){
	// Compute the number of free applications
	int nbFreeApps = 0;
	for(unsigned long i=0; i<nApps; i++)
	{
		if (freeApps[i])
		{
			nbFreeApps++;
		}
	}

	// We allocate the data structure for the subset of applis we focus on
	double actualVolRem[nbFreeApps];
	double actualReleaseDates[nbFreeApps];
	double actualProgresses[nbFreeApps];
	double actualMaxBandwidths[nbFreeApps];

	// We copy the data
	int index = 0;
	for(unsigned long i=0; i<nApps; i++)
	{
		if (freeApps[i])
		{
			actualVolRem[index]     = remCommVol[i];
			actualReleaseDates[index]   = releaseDate[i];
			actualProgresses[index] = progress[i];
			actualMaxBandwidths[index] = maxBandwidth[i];

			index++;
		}
	}

	return totalRequiredBandwidth(nbFreeApps, targetYield, delta, cTime,
			actualVolRem, actualReleaseDates, actualProgresses, actualMaxBandwidths);
}



// Computes the best mean yield that can be achieved at date cTime+delta
long double MinYield(double delta, double cTime,
		unsigned long nApps, double * remCommVol, double * releaseDate, double * progress, double * maxBandwidth,
		double availableBanwidth, double &totBand,
		bool * bwIsExhausted, bool * constrainingApps){

	if (ALLDEBUGNELEXMIN) cout << "\t\t\t\tMinYield computation delta= " << delta << " avail bw= " << availableBanwidth  << endl;

	if(ALLDEBUGNELEXMIN){
		for(unsigned long i = 0; i<nApps; i++){
			cout <<"\t\t\t\t\tAppli " << i << " max bw: " << maxBandwidth[i] << "\tRemaining vol: " << remCommVol[i] << endl;
		}
	}

	// Initialization of the data structure defining the type of solution found
	*bwIsExhausted = false;
	for(unsigned long i = 0; i<nApps; i++){
		constrainingApps[i] = false;
	}
	long double yields[nApps];


	// Compute an upper bound on the minimum yield
	long double upperBound = -1;
	for(unsigned long i = 0; i<nApps; i++)
	{
		yields[i] = fYieldMax(i, delta, cTime, remCommVol, releaseDate, progress, maxBandwidth);
		if ((i == 0) || (yields[i] < upperBound))
		{
			upperBound = yields[i];
		}
		if (ALLDEBUGNELEXMIN) cout << "\t\t\t\t\t " << i << "\tyield: " << yields[i] << "\tUpperBound: " << upperBound << endl;
	}

	if (ALLDEBUGNELEXMIN) cout << "\t\t\t\tMinYield: upperBound= " << upperBound << endl;


	if (ALLDEBUGNELEXMIN){
		long double tmpRequiredBW = totalRequiredBandwidth(nApps, upperBound, delta, cTime, remCommVol, releaseDate, progress, maxBandwidth);
		cout << "\t\t\t\tIs UpperBoundFeasible? " << tmpRequiredBW << "?<=?" <<  availableBanwidth+ERROR << endl;
	}

	// Check whether the upper bound defines a feasible solution
	if (totalRequiredBandwidth(nApps, upperBound, delta, cTime, remCommVol, releaseDate, progress, maxBandwidth) <= availableBanwidth+ERROR)
	{
		if (ALLDEBUGNELEXMIN) cout << "\t\t\t\tMinYield: upperBound is REACHABLE " << endl;
		*bwIsExhausted = false;
		for(unsigned long i = 0; i<nApps; i++){
			if(yields[i] <= upperBound+ERROR){
				constrainingApps[i] = true;
			}
		}
		return upperBound;
	}

	// Build and sort the set of yields achieved when no bandwidth is allocated
	vector<long double> yieldsWOBandwidth;
	for(unsigned long i = 0; i<nApps; i++)
	{
		yieldsWOBandwidth.push_back(fYield(i, delta, cTime, 0, releaseDate, progress, maxBandwidth));
	}
    sort(yieldsWOBandwidth.begin(), yieldsWOBandwidth.end());

	// Identification of applications which need bandwidth to reach the target yield
	bool applisNeedingBW[nApps];
    // Identification of the interval of ``yields without any bandwidth allocated'' including the optimal solution
    unsigned long indexMin = 0;
    unsigned long indexMax = yieldsWOBandwidth.size()-1;
	if (totalRequiredBandwidth(nApps, yieldsWOBandwidth[indexMax], delta, cTime, remCommVol, releaseDate, progress, maxBandwidth) <= availableBanwidth+ERROR)
	{
		indexMin = indexMax;
		for(unsigned long i = 0; i<nApps; i++)
		{
			applisNeedingBW[i] = true;
		}
	}
	else
	{
		while (indexMax-indexMin>1)
		{
			unsigned long indexMiddle = (indexMin+indexMax)/2;
			if (totalRequiredBandwidth(nApps, yieldsWOBandwidth[indexMiddle], delta, cTime, remCommVol, releaseDate, progress, maxBandwidth) > availableBanwidth-ERROR)
			{
				indexMax = indexMiddle;
			}
			else
			{
				indexMin = indexMiddle;
			}
		}
		for(unsigned long i = 0; i<nApps; i++)
		{
			if (fYield(i, delta, cTime, 0, releaseDate, progress, maxBandwidth) < yieldsWOBandwidth[indexMin]+ERROR)
			{
				applisNeedingBW[i] = true;
			}
			else
			{
				applisNeedingBW[i] = false;
			}
		}
	}

	// Computation of the optimal yield
	long double sumInNumerator = 0;
	long double sumInDenominator = 0;
	for(unsigned long i = 0; i<nApps; i++)
	{
//		cout << "\t\t\t\tAppli " << i;
//		if (applisNeedingBW[i]){
//			cout << " needs bandwidth." << endl;
//		}
//		else{
//			cout << " does not need bandwidth." << endl;
//
//		}

		if (applisNeedingBW[i])
		{
			sumInNumerator += progress[i]*maxBandwidth[i];
			sumInDenominator += (cTime+delta-releaseDate[i])*maxBandwidth[i];
		}
	}

	if (nApps==1){
		long double ywo = (sumInNumerator)/sumInDenominator;
		long double ywith = (delta*availableBanwidth+sumInNumerator)/sumInDenominator;
		cout << "Yield without bandwidth: " << ywo << endl;
		cout << "Yield with bandwidth: " << ywith << endl;
		cout << "Difference: " << ywith-ywo << endl;
	}


	long double yield = (delta*availableBanwidth+sumInNumerator)/sumInDenominator;

	*bwIsExhausted = true;

return yield;
}


// Computes the best mean yield that can be achieved at date cTime+delta for a subset of applications
long double MinYield(double delta, double cTime,
		unsigned long nApps,
		bool * freeApps,
		double * remCommVol, double * releaseDate, double * progress, double * maxBandwidth,
		double availableBanwidth, double &totBand,
		bool * bwIsExhausted, bool * constrainingApps){
	// Compute the number of free apps
	int nbFreeApps = 0;
	for(unsigned long i=0; i<nApps; i++)
	{
		if (freeApps[i])
		{
			nbFreeApps++;
		}
	}

	// We allocate the data structure for the subset of applications we focus on
	double actualVolRem[nbFreeApps];
	double actualReleaseDate[nbFreeApps];
	double actualProgress[nbFreeApps];
	double actualMaxBandwidth[nbFreeApps];
	bool tmpConstrainingApps[nbFreeApps];

	// We copy the data
	int index = 0;
	for(unsigned long i=0; i<nApps; i++)
	{
		if (freeApps[i])
		{
			actualVolRem[index]     = remCommVol[i];
			actualReleaseDate[index]   = releaseDate[i];
			actualProgress[index] = progress[i];
			actualMaxBandwidth[index] = maxBandwidth[i];

			index++;
		}
		constrainingApps[i] = false;
	}

	long double theYield =  MinYield(delta, cTime,
			nbFreeApps, actualVolRem, actualReleaseDate, actualProgress, actualMaxBandwidth,
			availableBanwidth, totBand,
			bwIsExhausted, tmpConstrainingApps);

	// We copy back the data
	index = 0;
	for(unsigned long i=0; i<nApps; i++)
	{
		if (freeApps[i])
		{
			constrainingApps[i] = tmpConstrainingApps[index];
			index++;
		}
	}

	return theYield;
}


bool LexMinYield(double delta, double cTime,
		unsigned long eventDefiningApp,
		unsigned long nApps, double * remCommVol, double * releaseDate, double * progress, double * maxBandwidth,
		double * allocatedBandwidths, double * yields, double &totBand,
		bool fullDEBUG){


	if (DEBUGNELEXMIN || fullDEBUG) {
		cout << endl << "\t\t***Entering LexMinYield*** with constraining application: " << eventDefiningApp;
		cout << " at time " << cTime << " + " << delta << " with debug " << fullDEBUG << endl;
	}

//	if(cTime > 5200){
//		fullDEBUG = true;
//	}

	assert(delta >= 0);

	// Amount of bandwidth that remains to be distributed
	double remainingBandwidth = totBand;
	// Data structure to hold the QUALITY of the solution
	vector<double> achievedYields;
	// Data structure to record for which applications the bandwidth (and thus the yield) has not already been fixed
	bool freeApplis[nApps];
	for(unsigned long i=0; i<nApps; i++)
	{
		freeApplis[i] = true;
	}
	int nbFreeApps = nApps;
	bool constrainingApps[nApps];

	// We fix the solution for the event defining application
	allocatedBandwidths[eventDefiningApp] = remCommVol[eventDefiningApp]/delta;
	remainingBandwidth -= allocatedBandwidths[eventDefiningApp];
	achievedYields.push_back(fYield(eventDefiningApp, delta, cTime, allocatedBandwidths[eventDefiningApp], releaseDate, progress, maxBandwidth));
	freeApplis[eventDefiningApp] = false;
	nbFreeApps--;

	if(fullDEBUG){
		cout << "Number of recorded yields before while: " << achievedYields.size() << endl;
	}

	// While there are free applications and bandwidth to be distributed
	while((nbFreeApps>0) && (remainingBandwidth > ERROR))
	{
		int formerNBFreeApps = nbFreeApps;

		if (DEBUGNELEXMIN || fullDEBUG)
			cout << "\t\t   remainingBandwidth= " << remainingBandwidth << " for " << nbFreeApps << " free applications" << endl;

		// If there is a simple free application the solution is obvious
		if (nbFreeApps==1){
			for(unsigned long i=0; i<nApps; i++)
			{
				if (freeApplis[i])
				{
					if (remainingBandwidth<maxBandwidth[i]){
						allocatedBandwidths[i] = remainingBandwidth;
					}
					else{
						allocatedBandwidths[i] = maxBandwidth[i];
					}
					if(remCommVol[i]/allocatedBandwidths[i] < delta){
						allocatedBandwidths[i] = remCommVol[i]/delta;
					}
					freeApplis[i] = false;
					nbFreeApps--;
					// We record its yield
					achievedYields.push_back(fYield(i, delta, cTime, allocatedBandwidths[i], releaseDate, progress, maxBandwidth));
				}
			}
		}
		else{

			// Compute the best achievable minimum yield for the free applications
			bool bwIsExhausted;
			long double yield = MinYield(delta, cTime,
					nApps,
					freeApplis,
					remCommVol,  releaseDate,  progress,  maxBandwidth,
					remainingBandwidth, totBand,
					&bwIsExhausted, constrainingApps);

			if (DEBUGNELEXMIN || fullDEBUG) cout << "\t\t\tTarget yield= " << yield << endl;

			// Compute the needed bandwidth
			double requiredBW = totalRequiredBandwidth(nApps,
					freeApplis,
					yield, delta, cTime,
					remCommVol, releaseDate, progress, maxBandwidth);

			//		// We should always use some of the bandwidth
			//		assert(requiredBW>0);

			// Sanity check
			if (DEBUGNELEXMIN || fullDEBUG) cout << "\t\t\trequiredBW= " << requiredBW << "\tremainingBandwidth= " << remainingBandwidth << "\tBandwidth left: " << remainingBandwidth-requiredBW<< endl;
			//		cout <<"\t\t\tError margin: " << nApps*ERROR << endl;
			//		assert(requiredBW <= remainingBandwidth+nApps*ERROR);
			// If the whole remaining bandwidth is used we are done
			//		if (requiredBW >= remainingBandwidth - nbFreeApps*ERROR)
			if (bwIsExhausted||(requiredBW >= remainingBandwidth - nbFreeApps*ERROR)||(remainingBandwidth -requiredBW < remainingBandwidth/10000))
			{
				if (DEBUGNELEXMIN || fullDEBUG){
					cout << "\t\t\tWhole bandwidth is used" << endl;
					if (bwIsExhausted & (!(requiredBW >= remainingBandwidth - nbFreeApps*ERROR))){
						cout << "\t\tPOTENTIAL PROBLEM: " << requiredBW << " ?>=? " << remainingBandwidth - nbFreeApps*ERROR << endl;
					}
				}
				long double scalingFactor = remainingBandwidth/requiredBW;
				for(unsigned long i=0; i<nApps; i++)
				{
					if (freeApplis[i])
					{
						// We set the bandwidth allocated to application i
						if (yield > fYield(i, delta, cTime, 0, releaseDate, progress, maxBandwidth))
						{
							allocatedBandwidths[i] = (maxBandwidth[i]/delta)*(yield*(cTime+delta-releaseDate[i])-progress[i])*scalingFactor;
							// To avoid rounding problems
							if (allocatedBandwidths[i] > maxBandwidth[i]){
								allocatedBandwidths[i] = maxBandwidth[i];
							}
							if (allocatedBandwidths[i] > remainingBandwidth){
								allocatedBandwidths[i] = remainingBandwidth;
							}
						}
						else
						{
							allocatedBandwidths[i] = 0;
						}
						if (DEBUGNELEXMIN || fullDEBUG) cout << "\t\t\t\tBW allocated to application "<<i<<": " << allocatedBandwidths[i] << endl;
						remainingBandwidth -= allocatedBandwidths[i];
						freeApplis[i] = false;
						nbFreeApps--;
						// We record its yield
						achievedYields.push_back(fYield(i, delta, cTime, allocatedBandwidths[i], releaseDate, progress, maxBandwidth));
					}
				}
				//			remainingBandwidth = 0;
			}
			// Otherwise we identify the applications whose yield is equal to the optimal minimum yield
			else
			{
				// If requiredBW is equal to 0, this means that the remaining bandwidth is so small that allocating it
				// does not change the minimum yield. The processing of this case is a hack to cope with rounding problems
				if(requiredBW==0){
					// We allocate the remaining bandwidth equally among the free apps
					int nbAllocatedApps = 0;
					for(unsigned long i=0; i<nApps; i++)
					{
						if (freeApplis[i])
						{
							allocatedBandwidths[i] = remainingBandwidth/nbFreeApps;
							freeApplis[i] = false;
							nbAllocatedApps++;
							// We record its yield
							achievedYields.push_back(fYield(i, delta, cTime, allocatedBandwidths[i], releaseDate, progress, maxBandwidth));
						}
					}
					assert(nbAllocatedApps==nbFreeApps);
					remainingBandwidth = 0;
					nbFreeApps=0;
				}
				else{
					if (DEBUGNELEXMIN || fullDEBUG) cout << "\t\t\tSome bandwidth will be left" << endl;
					for(unsigned long i=0; i<nApps; i++)
					{
						if (freeApplis[i])
						{
							// Compute the bandwidth needed by application i to achieve the target yield
							double bwForAppli = 0;
							if (yield > fYield(i, delta, cTime, 0, releaseDate, progress, maxBandwidth))
							{
								bwForAppli = (maxBandwidth[i]/delta)*(yield*(cTime+delta-releaseDate[i])-progress[i]);
								if(bwForAppli > maxBandwidth[i]){
									bwForAppli = maxBandwidth[i];
								}
								if(bwForAppli > remainingBandwidth){
									bwForAppli = remainingBandwidth;
								}
							}
							else
							{
								bwForAppli = 0;
							}
							// Sanity check
							assert(bwForAppli < maxBandwidth[i]+ERROR);

							if (DEBUGNELEXMIN || fullDEBUG)
							{
								cout << "\t\t\t\tApplication " << i << " BW: " << bwForAppli << " (<= max: " << maxBandwidth[i] << ") ";
								cout << "\tYield: " << fYield(i, delta, cTime, bwForAppli, releaseDate, progress, maxBandwidth);
								cout << "\tRemaining comm volume: " << remCommVol[i] - delta*bwForAppli << endl;
							}

							//					cout << "**MYDEBUG**:\t" << bwForAppli << "?>=?" << maxBandwidth[i]-ERROR << " \tOR\t " << bwForAppli*delta << "?>=?" <<  remCommVol[i]-ERROR << endl;


							// Check whether application i is constraining
							//					if ((bwForAppli >= maxBandwidth[i]-ERROR) || (bwForAppli*delta >= remCommVol[i]-ERROR) || (bwForAppli == 0))
							//						cout << "\t\t\t\t" << bwForAppli << " ?>=? " << maxBandwidth[i]-ERROR << "   OR   " << bwForAppli << "*" << delta << " = " << bwForAppli*delta << " ?>=? " << remCommVol[i]-ERROR << endl;
							if ((bwForAppli >= maxBandwidth[i]-ERROR) || (bwForAppli*delta >= remCommVol[i]-ERROR) || constrainingApps[i])
							{

								if (DEBUGNELEXMIN || fullDEBUG) cout << "\t\t\t\t\tAppli " << i << " is tight" << endl;
								// We set the bandwidth allocated to application i
								allocatedBandwidths[i] = bwForAppli;
								if (DEBUGNELEXMIN || fullDEBUG) cout << "\t\t\t\t\tAllocated bandwidth: " << allocatedBandwidths[i] << endl;
								remainingBandwidth -= allocatedBandwidths[i];
								freeApplis[i] = false;
								nbFreeApps--;
								// We record its yield
								achievedYields.push_back(fYield(i, delta, cTime, allocatedBandwidths[i], releaseDate, progress, maxBandwidth));
							}
						}
					}
				}
			}
		}
		//		cout << formerNBFreeApps << " > " << nbFreeApps << " OR " << remainingBandwidth << " < " << nbFreeApps*ERROR << endl;
		if(CATCHERROR){
			if ((formerNBFreeApps<=nbFreeApps)&(remainingBandwidth>nbFreeApps*ERROR)){
				return false;
			}
		}
		else{
			assert((formerNBFreeApps>nbFreeApps) || (remainingBandwidth<nbFreeApps*ERROR));
		}

	}

	if(fullDEBUG){
		cout << "Number of recorded yields after while: " << achievedYields.size() << endl;
	}

	// If there remain some free applications (because the bandwidth was exhausted and they did not need any)
	// we complete the solution
	if(nbFreeApps>0){
		for(unsigned long i=0; i<nApps; i++){
			if (freeApplis[i]){
				allocatedBandwidths[i] = 0;
				achievedYields.push_back(fYield(i, delta, cTime, allocatedBandwidths[i], releaseDate, progress, maxBandwidth));
			}
		}
	}

	long double totAllocatedBW = 0;
	if (DEBUGNELEXMIN || fullDEBUG){
		cout << "\t\tAllocated bandwidths: ";
		for(unsigned long i=0; i<nApps; i++){
			cout << i << ": " << allocatedBandwidths[i] << " (<=" << maxBandwidth[i] << ")  ";
			totAllocatedBW += allocatedBandwidths[i];
		}
		cout << endl;
		cout << "\t\tTotal allocation: " << totAllocatedBW << " (<=" << totBand << ")"<< endl;
		assert(totAllocatedBW < totBand + ERROR);
		//		assert(totAllocatedBW >= totBand -ERROR);
	}

	if (DEBUGNELEXMIN || fullDEBUG){
		cout << "\t\tAchieved yields: ";
		for(unsigned long i=0; i<nApps; i++){
			cout << i << ": " << fYield(i, delta, cTime, allocatedBandwidths[i], releaseDate, progress, maxBandwidth) << " ";
		}
		cout << "  [total allocated BW: " << totAllocatedBW << "] with event defining application: " << eventDefiningApp << endl;

		cout << "\t\tThe recorded achieved yields are :";
		for(unsigned long i=0; i<nApps; i++)
		{
			cout << achievedYields[i] << " ";
		}
		cout << endl;

	}

	// We sort the achieved yields
	sort(achievedYields.begin(), achievedYields.end());
	// We store the lexicographically sorted vector of achieved yields
	for(unsigned long i=0; i<nApps; i++)
	{
		yields[i] = achievedYields[i];
	}


	// We check whether the maxBandwidth constraints are satisfied
	bool problemDetected = false;
	for(unsigned long i=0; i<nApps; i++){
		if(allocatedBandwidths[i] > maxBandwidth[i] + ERROR){
			problemDetected = true;
		}
	}
	if(problemDetected){
		cout << "***** PROBLEM DETECTED *****" << endl;
	}
	if(problemDetected && (!fullDEBUG)){

		return LexMinYield(delta, cTime,
				eventDefiningApp,
				nApps, remCommVol, releaseDate, progress, maxBandwidth,
				allocatedBandwidths, yields, totBand, true);
	}

	return true;
}


bool LexMinYield(double delta, double cTime,
		unsigned long eventDefiningApp,
		unsigned long nApps, double * remCommVol, double * releaseDate, double * progress, double * maxBandwidth,
		double * allocatedBandwidths, double * yields, double &totBand){

	bool PECULIAR_CASE_DEBUG = false;
//	if ((cTime > 315) && (cTime<350)){
//		PECULIAR_CASE_DEBUG = true;
//	}


return LexMinYield(delta, cTime,
		eventDefiningApp,
		nApps, remCommVol, releaseDate, progress, maxBandwidth,
		allocatedBandwidths, yields, totBand, PECULIAR_CASE_DEBUG);
}


bool firstIntersection(double deltaMin, double deltaMax,
		long double firstNumLin, long double firstNumCnt,
		long double firstDenomLin, long double firstDenomCnt,
		long double secondNumLin, long double secondNumCnt,
		long double secondDenomLin, long double secondDenomCnt,
		long double * theFirstIntersection
		){
	// We compute the coefficients of the second degree polynomial
	long double polynA = firstNumLin*secondDenomLin - firstDenomLin*secondNumLin;
	long double polynB = firstNumLin*secondDenomCnt+firstNumCnt*secondDenomLin-firstDenomLin*secondNumCnt-firstDenomCnt*secondNumLin;
	long double polynC = firstNumCnt*secondDenomCnt-firstDenomCnt*secondNumCnt;

	// We check whether the second degree polynomial is not degenerated
	if (polynA!=0){
		long double discriminant = polynB*polynB-4*polynA*polynC;
		if (discriminant >= 0){
			// We compute and sort the two roots
			long double sqrtDiscriminant = sqrt(discriminant);
			long double firstRoot = (-polynB-sqrtDiscriminant)/(2*polynA);
			long double secondRoot = (-polynB+sqrtDiscriminant)/(2*polynA);
			if (firstRoot>secondRoot){
				long double tmpRoot = firstRoot;
				firstRoot = secondRoot;
				secondRoot = tmpRoot;
			}
			// We check whether the smallest root is in the target interval
			if ((deltaMin <= firstRoot) && (firstRoot <= deltaMax)){
				// It is in the interval this is THE solution
				*theFirstIntersection = firstRoot;
				return true;
			}
			else{
				if ((deltaMin <= secondRoot) && (secondRoot <= deltaMax)){
					// It is in the interval this is THE solution
					*theFirstIntersection = secondRoot;
					return true;
				}
				else{
					return false;
				}
			}
		}
		else{
			// There is no solution
			return false;
		}
	}
	else{
		// The polynomial is thus a first degree polynomial
		// We check whether the first degree polynomial is degenerate
		if (polynB!=0){
			long double singleRoot = -polynC/polynB;
			// We check whether the root is in the target interval
			if ((deltaMin <= singleRoot) && (singleRoot <= deltaMax)){
				// It is in the interval this is THE solution
				*theFirstIntersection = singleRoot;
				return true;
			}
			else{
				return false;
			}
		}
		else{
			// The difference is thus a constant
			if (polynC==0){
				// The two functions are identical
				*theFirstIntersection = deltaMax;
				return true;
			}
			else{
				return false;
			}
		}
	}
}

bool secondIntersection(double deltaMin, double deltaMax,
		long double firstNumLin, long double firstNumCnt,
		long double firstDenomLin, long double firstDenomCnt,
		long double secondNumLin, long double secondNumCnt,
		long double secondDenomLin, long double secondDenomCnt,
		long double * theFirstIntersection,
		bool * firstIsLarger
){
	// We compute the coefficients of the second degree polynomial
	long double polynA = firstNumLin*secondDenomLin - firstDenomLin*secondNumLin;
	long double polynB = firstNumLin*secondDenomCnt+firstNumCnt*secondDenomLin-firstDenomLin*secondNumCnt-firstDenomCnt*secondNumLin;
	long double polynC = firstNumCnt*secondDenomCnt-firstDenomCnt*secondNumCnt;

	// We check that the second degree polynomial is not degenerate
	if (polynA!=0){

		long double discriminant = polynB*polynB-4*polynA*polynC;
		// We know that deltaMin is a root. Hence, the discriminant cannot be negative
		assert(discriminant >= 0);
		assert(polynA != 0);

		// We compute and sort the two roots
		long double sqrtDiscriminant = sqrt(discriminant);
		long double firstRoot = (-polynB-sqrtDiscriminant)/(2*polynA);
		long double secondRoot = (-polynB+sqrtDiscriminant)/(2*polynA);
		if (firstRoot>secondRoot){
			long double tmpRoot = firstRoot;
			firstRoot = secondRoot;
			secondRoot = tmpRoot;
		}

		// We check whether the first root is equal to deltaMin (or almost)
		double distFirstRootToDeltaMin = firstRoot-deltaMin;
		if (distFirstRootToDeltaMin < 0 ){
			distFirstRootToDeltaMin = - distFirstRootToDeltaMin;
		}
		// We check whether the second root is equal to deltaMin (or almost)
		double distSecondRootToDeltaMin = secondRoot-deltaMin;
		if (distSecondRootToDeltaMin < 0 ){
			distSecondRootToDeltaMin = - distSecondRootToDeltaMin;
		}

		// One of the two roots should be (almost) equal to deltaMin
		if(DEBUGNELEXMIN){
			cout << "\t\tIn secondIntersection for interval [" << deltaMin << ", " << deltaMax << "].\n\tThe roots are:" << endl;
			cout << "\t\t\t" << firstRoot  << " (Hence, a distance to deltaMin of " << distFirstRootToDeltaMin << ")" << endl;
			cout << "\t\t\t" << secondRoot  << " (Hence, a distance to deltaMin of " << distSecondRootToDeltaMin << ")" << endl;
			cout << "\t\t\tpolynA = " << polynA << endl;
			cout << "\t\t\tpolynB = " << polynB << endl;
			cout << "\t\t\tpolynC = " << polynC << endl;
		}

		assert( (distFirstRootToDeltaMin <= ERROR) || (distFirstRootToDeltaMin <= SINGULARITY_PRECISION)
				||  (distSecondRootToDeltaMin <= ERROR) || (distSecondRootToDeltaMin <= SINGULARITY_PRECISION));

		// We check which of the two fraction/function is larger in the interval [deltaMin=firstRoot, min(secondRoot, deltaMax)]
		// If polynA is positive the difference is positive outside the roots and negative between them
		// IF THERE ARE TWO ROOTS
		if(discriminant > 0){
			// If the firstRoot is deltaMin
			if (distFirstRootToDeltaMin < distSecondRootToDeltaMin){
				if (polynA > 0){
					*firstIsLarger = false;
				}
				else{
					if (polynA < 0){
						*firstIsLarger = true;
					}
				}
			}
			// If the secondRoot is deltaMin
			else{
				if (polynA > 0){
					*firstIsLarger = true;
				}
				else{
					if (polynA < 0){
						*firstIsLarger = false;
					}
				}
			}
		}
		// If there is a double root
		else{
			if (polynA > 0){
				*firstIsLarger = true;
			}
			else{
				if (polynA < 0){
					*firstIsLarger = false;
				}
			}
		}


		*theFirstIntersection = deltaMax;
		// We check whether the second root is in the target interval
		// If deltaMin is the second root this is not the case
		if (distSecondRootToDeltaMin < distFirstRootToDeltaMin){
			return false;
		}
		else{
			if ((deltaMin <= secondRoot) && (secondRoot <= deltaMax)){
				// It is in the interval this is THE solution
				*theFirstIntersection = secondRoot;
				return true;
			}
			else{
				return false;
			}
		}
	}
	else{
		// The equation is a first degree polynomial
		// We know it has at least one root in deltaMin
		// We check whether the first degree polynomial is degenerate
		if (polynB!=0){
			// Then the only solution is in deltaMin
			return false;
		}
		else{
			// Then the two functions must be identical
			assert(polynC==0);
			*firstIsLarger = true; // Both assertions are true
			*theFirstIntersection = deltaMax;
			return true;
		}
	}
}


bool isYieldIncreasing(double cTime, int constrainingApp,
		unsigned long nApps, double * remCommVol, double * releaseDate, double * progress, double * maxBandwidth,
		double availableBandwidth, bool * applisNeedingBW){
	// We have identified the set of applications to which bandwidth should be allocated if delta = deltaMin
	// We check whether the common yield achieved by this function is increasing or decreasing
	long double firstSum = 0;
	long double secondSum = 0;
	long double thirdSum = 0;
	for(unsigned long i = 0; i<nApps; i++)
	{
		if (applisNeedingBW[i] && (i!=constrainingApp))
		{
			firstSum += maxBandwidth[i];
			secondSum += maxBandwidth[i]*progress[i];
			thirdSum += maxBandwidth[i]*(availableBandwidth*(cTime-releaseDate[i])+remCommVol[constrainingApp]);
		}
	}
	long double trendIndicator = firstSum*secondSum-thirdSum;
	// If the yield is decreasing along the interval the best solution is found at deltaMin
	if (trendIndicator>=0){
		return false;
	}
	else{
		return true;
	}
}

void computeYields(double delta, double cTime, int constrainingApp,
		unsigned long nApps, double * remCommVol, double * releaseDate, double * progress, double * maxBandwidth,
		double availableBandwidth, bool * applisNeedingBW,
		double * allocatedBandwidth, double * yields){
	// Initializations
	applisNeedingBW[constrainingApp] = false;
	for(unsigned long i = 0; i<nApps; i++){
		allocatedBandwidth[i] = 0;
	}

	// Variable to enable the debug of a peculiar case
	bool PECULIAR_CASE_DEBUG = false;
//	if (cTime > 26541){
//		PECULIAR_CASE_DEBUG = true;
//	}

	if(PECULIAR_CASE_DEBUG){
		cout << "IN computeYields for time " << delta << endl;
	}

	if (PECULIAR_CASE_DEBUG) cout << "\t\t\tallocatedBandwidth[constrainingApp] = " << allocatedBandwidth[constrainingApp] <<  " (computeYields)" << endl;

	allocatedBandwidth[constrainingApp] = remCommVol[constrainingApp]/delta;
	long double remainingBandwidth = availableBandwidth - allocatedBandwidth[constrainingApp];

//	cout << "\t\t\tallocatedBandwidth[constrainingApp] = " << allocatedBandwidth[constrainingApp] <<  " (computeYields)" << endl;


	// Computation of the optimal yield
	long double sumInNumerator = 0;
	long double sumInDenominator = 0;
	for(unsigned long i = 0; i<nApps; i++)
	{
		if (applisNeedingBW[i])
		{
			sumInNumerator += progress[i]*maxBandwidth[i];
			sumInDenominator += (cTime+delta-releaseDate[i])*maxBandwidth[i];
		}
	}
	long double yield = (delta*remainingBandwidth+sumInNumerator)/sumInDenominator;

	int nbAllocatedApps = 0;
	int lastAllocatedApp = -1;
	long double totalBWAllocated = 0;

	for(unsigned long i = 0; i<nApps; i++)
	{
		if (PECULIAR_CASE_DEBUG) cout << "Application " << i <<  " needing BW " << applisNeedingBW[i] << " allocated " << allocatedBandwidth[i]<< endl;

		if (ALLDEBUGNELEXMIN || PECULIAR_CASE_DEBUG)
		{
			cout << "\t\t\tAbout application: " << i << " yield WO BW = " << fYield(i, delta, cTime, 0, releaseDate, progress, maxBandwidth) << endl;
			double maxBW = remainingBandwidth;
			if (maxBW > maxBandwidth[i]){
				maxBW = maxBandwidth[i];
			}
			cout << "\t\t\t\t yield with all the bandwidth = " << fYield(i, delta, cTime, maxBW, releaseDate, progress, maxBandwidth) << endl;
		}
		if ((remainingBandwidth > ERROR) && (applisNeedingBW[i]))
		{
			if (ALLDEBUGNELEXMIN || PECULIAR_CASE_DEBUG)
				cout << "\t\t\t\t need BW" << endl;
			allocatedBandwidth[i] = requiredBandwidth(i, yield, delta, cTime, remCommVol, releaseDate, progress, maxBandwidth);
			if(allocatedBandwidth[i] > maxBandwidth[i]){
				allocatedBandwidth[i] = maxBandwidth[i];
			}
			totalBWAllocated += allocatedBandwidth[i];
			if (ALLDEBUGNELEXMIN || PECULIAR_CASE_DEBUG)
				cout << "\t\t\t\t allocated: " << allocatedBandwidth[i] << endl;
			yields[i] = yield;
			nbAllocatedApps++;
			lastAllocatedApp = i;
		}
		else{
//			cout <<"***before fYield" << allocatedBandwidth[i] << endl;
			if (ALLDEBUGNELEXMIN || PECULIAR_CASE_DEBUG)
				cout << "\t\t\t\t does not need BW" << endl;
			yields[i] = fYield(i, delta, cTime, allocatedBandwidth[i], releaseDate, progress, maxBandwidth);
		}
	}

	if (ALLDEBUGNELEXMIN || PECULIAR_CASE_DEBUG)
	{
		cout << "\t\t\tnbAllocatedApps = " << nbAllocatedApps << " totalBWAllocated=  " << totalBWAllocated << " for " << remainingBandwidth << endl;
	}

	if(remainingBandwidth>ERROR){
		assert(nbAllocatedApps>0);
		if (nbAllocatedApps==1){
			allocatedBandwidth[lastAllocatedApp] = remainingBandwidth;
			if(allocatedBandwidth[lastAllocatedApp] > maxBandwidth[lastAllocatedApp]){
				allocatedBandwidth[lastAllocatedApp] = maxBandwidth[lastAllocatedApp];
			}
		}
		else{
			if (totalBWAllocated==0){
				for(unsigned long i = 0; i<nApps; i++){
					if (applisNeedingBW[i])
					{
						allocatedBandwidth[i] = remainingBandwidth/nbAllocatedApps;
						if(allocatedBandwidth[i] > maxBandwidth[i]){
							allocatedBandwidth[i] = maxBandwidth[i];
						}
						yields[i] = fYield(i, delta, cTime, allocatedBandwidth[i], releaseDate, progress, maxBandwidth);
					}
				}
			}
			else{
				if (totalBWAllocated > remainingBandwidth){
					long double correctionRatio = remainingBandwidth/totalBWAllocated;
					for(unsigned long i = 0; i<nApps; i++){
						if (applisNeedingBW[i])
						{
							allocatedBandwidth[i] = allocatedBandwidth[i]*correctionRatio;
							if(allocatedBandwidth[i] > maxBandwidth[i]){
								allocatedBandwidth[i] = maxBandwidth[i];
							}
							yields[i] = fYield(i, delta, cTime, allocatedBandwidth[i], releaseDate, progress, maxBandwidth);
						}
					}
				}
				else{
					if ((remainingBandwidth > ERROR) && (remainingBandwidth-totalBWAllocated > remainingBandwidth/INV_PRECISION)){
						long double correctionRatio = remainingBandwidth/totalBWAllocated;
						for(unsigned long i = 0; i<nApps; i++){
							if (applisNeedingBW[i])
							{
								allocatedBandwidth[i] = allocatedBandwidth[i]*correctionRatio;
								if(allocatedBandwidth[i] > maxBandwidth[i]){
									allocatedBandwidth[i] = maxBandwidth[i];
								}
								yields[i] = fYield(i, delta, cTime, allocatedBandwidth[i], releaseDate, progress, maxBandwidth);
							}
						}
					}
				}
			}
		}
	}
}

// Computes the best mean yield that can be achieved between the dates cTime+firstBound and cTime+secondBound if application ``constraining''
// is the one defining the event
 void MinYieldAtAppEvent(double delta, double cTime,
		int constrainingApp,
		unsigned long nApps, double * remCommVol, double * releaseDate, double * progress, double * maxBandwidth,
		double availableBandwidth,
		double * allocatedBandwidth, double * yields){

	if (ALLDEBUGNELEXMIN) cout << "\t\t\t\tMinYieldAtAppEvent computation at time [ " << delta << "] avail bw= " << availableBandwidth  << endl;

	if(ALLDEBUGNELEXMIN){
		for(unsigned long i = 0; i<nApps; i++){
			cout <<"\t\t\t\t\tAppli " << i << " max bw: " << maxBandwidth[i] << "\tRemaining vol: " << remCommVol[i] << endl;
		}
	}



	// We initialize the data structure which will hold the data allocation
	for(int i=0; i<nApps; i++){
		allocatedBandwidth[i] = 0;
	}
	allocatedBandwidth[constrainingApp] = remCommVol[constrainingApp]/delta;

	double remainingBandwidth = availableBandwidth - allocatedBandwidth[constrainingApp];
	if (remainingBandwidth < ERROR){
		remainingBandwidth = 0;
	}

	// Build and sort the set of yields achieved when no bandwidth is allocated
	vector<long double> yieldsWOBandwidth;
	for(unsigned long i = 0; i<nApps; i++)
	{
		if (i != constrainingApp) yieldsWOBandwidth.push_back(fYield(i, delta, cTime, 0, releaseDate, progress, maxBandwidth));
	}
    sort(yieldsWOBandwidth.begin(), yieldsWOBandwidth.end());

	// Identification of applications which need bandwidth to reach the target yield
	bool applisNeedingBW[nApps];
    // Identification of the interval of ``yields without any bandwidth allocated'' including the optimal solution
    unsigned long indexMin = 0;
    unsigned long indexMax = yieldsWOBandwidth.size()-1;

//    cout << "before if" << endl;

	if (totalRequiredBandwidth(nApps, yieldsWOBandwidth[indexMax], delta, cTime, remCommVol, releaseDate, progress, maxBandwidth, constrainingApp) <= remainingBandwidth+ERROR)
	{
//	    cout << "in then " << endl;


		indexMin = indexMax;
		for(unsigned long i = 0; i<nApps; i++)
		{
			applisNeedingBW[i] = true;
		}
		applisNeedingBW[constrainingApp] = false;
	}
	else
	{
//	    cout << "in else " << endl;
//
//		cout << "\tindexMin= " << indexMin << "\tindexMax= " << indexMax << endl;

		while (indexMax-indexMin>1)
		{

//			cout << "indexMin= " << indexMin << "\tindexMax= " << indexMax << endl;

			unsigned long indexMiddle = (indexMin+indexMax)/2;
			if (totalRequiredBandwidth(nApps, yieldsWOBandwidth[indexMiddle], delta, cTime, remCommVol, releaseDate, progress, maxBandwidth,constrainingApp) > remainingBandwidth-ERROR)
			{
				indexMax = indexMiddle;
			}
			else
			{
				indexMin = indexMiddle;
			}
		}
		for(unsigned long i = 0; i<nApps; i++)
		{
			if(i!=constrainingApp){
				if (fYield(i, delta, cTime, 0, releaseDate, progress, maxBandwidth) < yieldsWOBandwidth[indexMin]+ERROR)
				{
					applisNeedingBW[i] = true;
				}
				else
				{
					applisNeedingBW[i] = false;
				}
			}
		}
	}

	// Sanity check
	applisNeedingBW[constrainingApp] = false;

	// Computation of the optimal yield
	long double sumInNumerator = 0;
	long double sumInDenominator = 0;
	for(unsigned long i = 0; i<nApps; i++)
	{
		if (applisNeedingBW[i])
		{
			sumInNumerator += progress[i]*maxBandwidth[i];
			sumInDenominator += (cTime+delta-releaseDate[i])*maxBandwidth[i];
		}
	}
	long double yield = (delta*remainingBandwidth+sumInNumerator)/sumInDenominator;

//cout << "\t\tYield = " << yield << endl;

	for(unsigned long i = 0; i<nApps; i++)
	{
//		cout << "About application: " << i << endl;
		if (applisNeedingBW[i])
		{
//			cout << "\t need BW" << endl;
			allocatedBandwidth[i] = requiredBandwidth(i, yield, delta, cTime, remCommVol, releaseDate, progress, maxBandwidth);
			if(allocatedBandwidth[i] > maxBandwidth[i]){
				allocatedBandwidth[i] = maxBandwidth[i];
			}
//			cout << "\t allocated: " << allocatedBandwidth[i] << endl;
			yields[i] = yield;
		}
		else{
//			cout << "\t does not need BW" << endl;
			yields[i] = fYield(i, delta, cTime, allocatedBandwidth[i], releaseDate, progress, maxBandwidth);
		}
	}
//	cout << "Everything is OK..." << endl;
}

void MinYieldAtAppEventInInterval(double deltaMin, double deltaMax, double cTime,
		int constrainingApp,
		unsigned long nApps, double * remCommVol, double * releaseDate, double * progress, double * maxBandwidth,
		double availableBandwidth,
		double * allocatedBandwidth, double * yields, double * bestDelta, int recurDepth){

	// Boolean variable to enable the debugging of a peculiar case
	bool PECULIAR_CASE_DEBUG = false;
//	if ((cTime > 315) && (cTime<350)){
//				PECULIAR_CASE_DEBUG = true;
//	}

//	if ((cTime>5210) && (deltaMin > 33800) && (constrainingApp==57)){
//	if ((cTime>10880)){
//		PECULIAR_CASE_DEBUG = true;
//	}
//	if (cTime > 26541){
//		PECULIAR_CASE_DEBUG = true;
//	}


	if ((DEBUGNELEXMIN) || PECULIAR_CASE_DEBUG || CHECKOPTIMISATION)
	{
		cout << "\tMinYieldAtAppEventInInterval computation in interval [ " << deltaMin << ", " << deltaMax << "] for app " << constrainingApp;
		cout << " avail bw= " << availableBandwidth  << " recursive level: " << recurDepth;
		cout << endl;
	}

	if(DEBUGNELEXMIN){
		for(unsigned long i = 0; i<nApps; i++){
			cout <<"\t\tAppli " << i << " max bw: " << maxBandwidth[i] << "\tRemaining vol: " << remCommVol[i] << endl;
		}
	}

	// We initialize the data structure which will hold the data allocation
	for(int i=0; i<nApps; i++){
		allocatedBandwidth[i] = 0;
	}
	allocatedBandwidth[constrainingApp] = remCommVol[constrainingApp]/deltaMin;

//	if (PECULIAR_CASE_DEBUG){
//		cout << "BW at step 01 : " << allocatedBandwidth[constrainingApp] << endl;
//	}
	if (allocatedBandwidth[constrainingApp] > maxBandwidth[constrainingApp] + ERROR){
		cout << "***ALLOCATION PROBLEM " << allocatedBandwidth[constrainingApp] << "vs. " << maxBandwidth[constrainingApp] << endl;
	}
//	cout << "\t\t\tallocatedBandwidth[constrainingApp] = " << allocatedBandwidth[constrainingApp] << endl;

	double remainingBandwidth = availableBandwidth - allocatedBandwidth[constrainingApp];
	if (remainingBandwidth < ERROR){
		remainingBandwidth = 0;
	}
	if(DEBUGNELEXMIN || PECULIAR_CASE_DEBUG) cout << "\t\tRemaining bandwidth at time deltaMin: " << remainingBandwidth << endl;

	// Data structure to identify applications which need bandwidth to reach the optimal yield
	bool applisNeedingBW[nApps];
	int nbApplisNeedingBW = 0;
	for(unsigned long i = 0; i<nApps; i++)
	{
		applisNeedingBW[i] = false;
	}

	if(remCommVol[constrainingApp]<0.001){
		allocatedBandwidth[constrainingApp] = remCommVol[constrainingApp]/deltaMin;
    	double totRemainingRequestedBW = 0;
		for(unsigned long i = 0; i<nApps; i++){
			if(i!=constrainingApp){
				totRemainingRequestedBW += maxBandwidth[i];
			}
		}

		if(totRemainingRequestedBW>remainingBandwidth){
			for(unsigned long i = 0; i<nApps; i++){
				if(i!=constrainingApp){
					allocatedBandwidth[i] = remainingBandwidth*maxBandwidth[i]/totRemainingRequestedBW;
					applisNeedingBW[i] = true;
				}
			}
		}
		else{
			for(unsigned long i = 0; i<nApps; i++){
				if(i!=constrainingApp){
					allocatedBandwidth[i] = maxBandwidth[i];
					applisNeedingBW[i] = true;
				}
			}
		}

//		cout << "BW at step 02 : " << allocatedBandwidth[constrainingApp] << endl;

		*bestDelta = deltaMin;
		computeYields(*bestDelta, cTime, constrainingApp,
				nApps, remCommVol, releaseDate, progress, maxBandwidth,
				availableBandwidth, applisNeedingBW,
				allocatedBandwidth, yields);

//		cout << "BW at step 03 : " << allocatedBandwidth[constrainingApp] << endl;

		return;
	}



	// Build and sort the set of yields achieved when no bandwidth is allocated
	vector<pair<long double, int> > yieldsWOBandwidth;
	for(unsigned long i = 0; i<nApps; i++)
	{
		if (i != constrainingApp) yieldsWOBandwidth.push_back(make_pair(fYield(i, deltaMin, cTime, 0, releaseDate, progress, maxBandwidth), i));
	}
    sort(yieldsWOBandwidth.begin(), yieldsWOBandwidth.end());

    if (ALLDEBUGNELEXMIN || PECULIAR_CASE_DEBUG){
    	for(int i = 0; i<yieldsWOBandwidth.size(); i++){
    		cout << "\t\tApplication " << get<1>(yieldsWOBandwidth[i]) << " as a base yield of " << get<0>(yieldsWOBandwidth[i]) << endl;
    	}
    }

    // Computation of an upper bound on the maximum achievable yield
    double upperBoundYield = -1;
    double currentYield = fYieldAtCompletion(constrainingApp, deltaMin, cTime, remCommVol, releaseDate, progress, maxBandwidth);
    for(unsigned long i = 0; i<nApps; i++)
    {
    	if(i!=constrainingApp){
    		currentYield = fYieldMax(i, deltaMin, cTime, remCommVol, releaseDate, progress, maxBandwidth);
    	}
    	if ((i==0) || (currentYield < upperBoundYield)){
    		upperBoundYield = currentYield;
    	}
    }
    if (DEBUGNELEXMIN || PECULIAR_CASE_DEBUG) cout << "\t\tUpper bound yield: " << upperBoundYield << endl;

    // We check whether this upper bound is achievable
	long double requiredBWForUpperBound = totalRequiredBandwidth(nApps,
			upperBoundYield, deltaMin, cTime,
			remCommVol, releaseDate, progress, maxBandwidth, constrainingApp);
    if (requiredBWForUpperBound -remainingBandwidth> remainingBandwidth/INV_PRECISION)
    {// The upper bound is not achievable

    	if (DEBUGNELEXMIN || PECULIAR_CASE_DEBUG){
    		cout << "\t\tThe upper bound is not achievable: required " << requiredBWForUpperBound << " available " << remainingBandwidth << endl;
    	}

        // Identification of the interval of ``yields without any bandwidth allocated'' including the optimal solution
        unsigned long indexMin = 0;
        unsigned long indexMax = yieldsWOBandwidth.size()-1;

//    	long double delta = deltaMin; // The optimal event date in the interval

    	// We check first if the largest yield-achieved-wo-bandwidth is achievable by all applications
    	if (totalRequiredBandwidth(nApps, get<0>(yieldsWOBandwidth[indexMax]), deltaMin, cTime, remCommVol, releaseDate, progress, maxBandwidth, constrainingApp)
    			- remainingBandwidth <= remainingBandwidth/INV_PRECISION)
    	{
    		if (ALLDEBUGNELEXMIN || PECULIAR_CASE_DEBUG) cout << "\t\t\tAll applications need bandwidth" << endl;

    		indexMin = indexMax;
    		for(unsigned long i = 0; i<nApps; i++)
    		{
    			applisNeedingBW[i] = true;
    			nbApplisNeedingBW++;
    		}
    		applisNeedingBW[constrainingApp] = false;
    	}
    	else
    	{
    		if (DEBUGNELEXMIN || PECULIAR_CASE_DEBUG) {
    			long double tmpTORBW = totalRequiredBandwidth(nApps, get<0>(yieldsWOBandwidth[indexMax]), deltaMin, cTime, remCommVol, releaseDate, progress, maxBandwidth, constrainingApp) ;
    			cout << "\t\t\t\t\t\tTotal required BW= " << tmpTORBW << "\tremainingBandwidth= " << remainingBandwidth << endl;
    		}
    		if (DEBUGNELEXMIN || PECULIAR_CASE_DEBUG) cout << "\t\t\t\t\t\tindexMin= " << indexMin << "\tindexMax= " << indexMax << endl;

    		while (indexMax-indexMin>1)
    		{
    			if (ALLDEBUGNELEXMIN) cout << "\t\t\t\t\t\t\tindexMin= " << indexMin << "\tindexMax= " << indexMax << endl;

    			unsigned long indexMiddle = (indexMin+indexMax)/2;
    			long double totalRequiredBWForMiddleIndex = totalRequiredBandwidth(nApps, get<0>(yieldsWOBandwidth[indexMiddle]), deltaMin, cTime, remCommVol, releaseDate, progress, maxBandwidth,constrainingApp);
    			if (DEBUGNELEXMIN || PECULIAR_CASE_DEBUG) cout << "\t\t\t\t\t\t\tBW excess: " << totalRequiredBWForMiddleIndex - remainingBandwidth << " precision: " << remainingBandwidth/INV_PRECISION << endl;
    			if (DEBUGNELEXMIN || PECULIAR_CASE_DEBUG) cout << "\t\t\t\t\t\t\tRequired: " << totalRequiredBWForMiddleIndex << " available: " << remainingBandwidth << endl;

//    			if (totalRequiredBWForMiddleIndex - remainingBandwidth > ERROR)
//    			{
//    				indexMax = indexMiddle;
//    			}
//    			else
//    			{
//    				indexMin = indexMiddle;
//    			}
    			if (totalRequiredBWForMiddleIndex - remainingBandwidth > remainingBandwidth/INV_PRECISION)
    			{
    				indexMax = indexMiddle;
    			}
    			else
    			{
    				indexMin = indexMiddle;
    			}
//    			if (totalRequiredBWForMiddleIndex - remainingBandwidth > -remainingBandwidth/INV_PRECISION)
//    			{
//    				indexMax = indexMiddle;
//    			}
//    			else
//    			{
//    				indexMin = indexMiddle;
//    			}
    		}

			if (ALLDEBUGNELEXMIN || PECULIAR_CASE_DEBUG) cout << "\t\t\t\t\t\t\tindexMin= " << indexMin << "\tindexMax= " << indexMax << endl;

    		long double refYield = get<0>(yieldsWOBandwidth[indexMin]);
			if (ALLDEBUGNELEXMIN || PECULIAR_CASE_DEBUG) cout << "\t\t\t\t\t\t\trefYield= " << refYield << endl;

    		for(unsigned long i = 0; i<nApps; i++)
    		{
    			if(i!=constrainingApp){
    				long double tmpYield = fYield(i, deltaMin, cTime, 0, releaseDate, progress, maxBandwidth);
    				if (ALLDEBUGNELEXMIN || PECULIAR_CASE_DEBUG) {
    					cout << "\t\tAppli " << i << " yield(0) = " << tmpYield << " diff: " << tmpYield -refYield << " vs. " << refYield/INV_PRECISION << endl;
    				}
    				if (tmpYield - refYield <= refYield/INV_PRECISION)
    				{
    					applisNeedingBW[i] = true;
    					nbApplisNeedingBW++;
    				}
    				else
    				{
    					applisNeedingBW[i] = false;
    				}
    			}
    		}
    	}

    	if (DEBUGNELEXMIN || PECULIAR_CASE_DEBUG){
    		cout << "\t\tApplis needing bandwidth: ";
    		for(unsigned long i = 0; i<nApps; i++){
    			cout << i << ":" << applisNeedingBW[i] << "  ";
    		}
    		cout << endl;
    	}

    	// We check whether the function y^opt_I is increasing
    	bool yieldIsIncreasing = isYieldIncreasing(cTime, constrainingApp,
    			nApps, remCommVol, releaseDate, progress, maxBandwidth,
				availableBandwidth, applisNeedingBW);
    	// If y^opt_I is not increasing, we conclude
    	if(!yieldIsIncreasing){
    		*bestDelta = deltaMin;
    		computeYields(*bestDelta, cTime, constrainingApp,
    				nApps, remCommVol, releaseDate, progress, maxBandwidth,
					availableBandwidth, applisNeedingBW,
					allocatedBandwidth, yields);
    	}
    	// y^opt_I is increasing
    	else{
    		// Data structure used to store the last date at which the solution is valid, and how the structure of the solution changes
    		bool intersectionFound = false;
    		long double lastValidDelta = deltaMax;
    		bool firstIntersectionIsWithNextApp = false;
    		bool firstIntersectionIsWithConstrainingApp = false;
    		bool firstIntersectionIsWithBandwidthBound = false;

    		// We compute the coefficient of the common yield function of the applications receiving bandwidth
    		long double commonYieldNumLin = availableBandwidth;
    		long double commonYieldNumCnt = -remCommVol[constrainingApp];
    		long double commonYieldDenomLin = 0;
    		long double commonYieldDenomCnt = 0;

    		for(unsigned long i = 0; i<nApps; i++)
    		{
    			if (applisNeedingBW[i])
    			{
    				commonYieldNumCnt += progress[i]*maxBandwidth[i];
    				commonYieldDenomLin += maxBandwidth[i];
    				commonYieldDenomCnt += (cTime-releaseDate[i])*maxBandwidth[i];
    			}
    		}


    		bool inSingularityCase = false;
			bool maxYieldIsLarger = false;
			long double deltaSecondIntersection = deltaMax;
    		bool solutionFound = false;
    		bool secondIntersectionFound = false;
    		int boundingApp = 0;
			long double  boundingAppYieldAtDeltaMax = 0;
			bool actualIntersectionFound = false;

    		// We compute the potential intersection with the yield function of an application receiving its maximum bandwidth
    		for(unsigned long i = 0; i<nApps; i++)
    		{
    			if((i!=constrainingApp)&&(applisNeedingBW[i])){
    	    		// We compute the coefficients of the maximum yield function
    				// First case: the maximum yield is achieved by receiving the maximum bandwidth
    				long double maximumYieldNumLin = 1;
    				long double maximumYieldNumCnt = progress[i];
    				long double maximumYieldDenomLin = 1;
    				long double maximumYieldDenomCnt = cTime-releaseDate[i];
    				// Second case: the maximum yield is achieved by the end of the communication
    				if((deltaMin+deltaMax)/2 > remCommVol[i]/maxBandwidth[i]){
    					maximumYieldNumLin = 0;
    					maximumYieldNumCnt = progress[i]+remCommVol[i]/maxBandwidth[i];
    				}

    				long double deltaIntersectWithMaxYield = deltaMax;
    				bool intersectionFoundWithMaxYield = firstIntersection(deltaMin, deltaMax,
    	    				commonYieldNumLin, commonYieldNumCnt,
    						commonYieldDenomLin, commonYieldDenomCnt,
							maximumYieldNumLin, maximumYieldNumCnt,
							maximumYieldDenomLin, maximumYieldDenomCnt,
    						&deltaIntersectWithMaxYield);

//    				if (PECULIAR_CASE_DEBUG){
//    					if (intersectionFoundWithMaxYield)
//    						cout << "\t\tintersectionFoundWithMaxYield FOUND at time " << deltaIntersectWithMaxYield << " for appli : " << i << endl;
//    					else
//    						cout << "\t\tintersectionFoundWithMaxYield not found for appli " << i << endl;
//    				}

    	    		if ((intersectionFoundWithMaxYield) &&
    	    				((!intersectionFound) || (deltaIntersectWithMaxYield < lastValidDelta))){
    	    			// We check whether the intersection found defines a singularity time
    	        		// If the intersection of the two yield functions is at the beginning of the interval (or almost)
    	        		// then the beginning of the interval is a singularity time where the common yield of the applications needing bandwidth
    	        		// is equal to the yield of the application defining the upper-bound on the minimum yield.
    	        		// In that case, *inside* the interval the common yield may define the minimum yield

    	        		if ((deltaIntersectWithMaxYield - deltaMin < SINGULARITY_PRECISION)&&(deltaIntersectWithMaxYield - deltaMin < deltaMin/INV_PRECISION)){
    	        			// We are in the singularity case: we record this and the application defining it
    	        			// If we had already identified that we were in the singularity we check whether the current application implies stringer constraints
    	        			long double tmpYieldAtDmax = fYieldMax(i, deltaMax, cTime, remCommVol, releaseDate, progress, maxBandwidth);

    	        			if((!inSingularityCase) || (tmpYieldAtDmax < boundingAppYieldAtDeltaMax)){
    	        				boundingAppYieldAtDeltaMax = tmpYieldAtDmax;
    	        				boundingApp = i;
    	        				inSingularityCase = true;

    	        				if(DEBUGNELEXMIN || PECULIAR_CASE_DEBUG){
    	        					cout << "\t\t\tIn singularity case: ";
    	        					if(inSingularityCase) {
    	        						cout << " YES " << endl;
    	        					}
    	        					else{
    	        						cout << " no " << endl;
    	        					}
    	        				}
    	        				// If the common yield function is non-increasing, we can conclude
    	        				bool yieldIsIncreasing = isYieldIncreasing(cTime, constrainingApp,
    	        						nApps, remCommVol, releaseDate, progress, maxBandwidth,
										availableBandwidth, applisNeedingBW);
    	        				if(!yieldIsIncreasing){
    	        					// This case should never happen... (as it has already checked for earlier)
    	        					if (PECULIAR_CASE_DEBUG) cout << "In the impossible case" << endl;
    	        					*bestDelta = deltaMin;
    	        					computeYields(*bestDelta, cTime, constrainingApp,
    	        							nApps, remCommVol, releaseDate, progress, maxBandwidth,
											availableBandwidth, applisNeedingBW,
											allocatedBandwidth, yields);
    	        					solutionFound = true;
    	        				}
    	        				else{
    	        					solutionFound = false;
    	        					// We check which function is the largest between the common-yield function and the maximum yield function
    	        					secondIntersectionFound = secondIntersection(deltaMin, deltaMax,
    	        							maximumYieldNumLin, maximumYieldNumCnt,
											maximumYieldDenomLin, maximumYieldDenomCnt,
											commonYieldNumLin, commonYieldNumCnt,
											commonYieldDenomLin, commonYieldDenomCnt,
											&deltaSecondIntersection, &maxYieldIsLarger);
    	        					if(DEBUGNELEXMIN || PECULIAR_CASE_DEBUG){
    	        						cout << "\t\t\tsecondIntersectionFound = " << secondIntersectionFound << endl;
    	        						cout << "\t\t\tat time : " << deltaSecondIntersection << endl;
    	        						cout << "\t\t\tmaxYieldIsLarger : " << maxYieldIsLarger << endl;
    	        					}
    	        					if(secondIntersectionFound){
    	        						intersectionFound = true;
    	        						if(deltaSecondIntersection < lastValidDelta){
    	        							lastValidDelta = deltaSecondIntersection;
    	        	    	    			firstIntersectionIsWithNextApp = false;
    	        	    	    			firstIntersectionIsWithConstrainingApp = false;
    	        	    	    			firstIntersectionIsWithBandwidthBound = true;
    	        						}
    	        					}
    	        				}
    	        			}
    	        		}
    	        		// We are not at a singularity. Hence, the intersection is an actual intersection time and we record it
    	        		else{
        	    			lastValidDelta = deltaIntersectWithMaxYield;
        	    			actualIntersectionFound = true;
        	    			firstIntersectionIsWithNextApp = false;
        	    			firstIntersectionIsWithConstrainingApp = false;
        	    			firstIntersectionIsWithBandwidthBound = true;
        	    			intersectionFound = true;
    	        		}
    	    		}

    	    		if(DEBUGNELEXMIN || PECULIAR_CASE_DEBUG){
    	    			if(intersectionFoundWithMaxYield){
    	    				cout << "\t\tIntersection between maximum yield function (for application " << i << ") and common yield found at time " << deltaIntersectWithMaxYield;
    	    				cout << " in [" << deltaMin << ", " << deltaMax << "]" << endl;
    	    			}
    	    			else{
    	    				cout << "\t\tIntersection between maximum yield function (for application " << i << ") and common yield NOT found." << endl;
    	    			}
    	    		}
    			}
    		}

    		if(actualIntersectionFound && (!inSingularityCase)){
    			intersectionFound = true;
    		}

    		// If we are not in the singularity case or if no solution was found, we continue to look for the best solution
    		if (!solutionFound){
	    		// We compute the coefficients of the maximum yield function
				long double maximumYieldNumLin = 1;
				long double maximumYieldNumCnt = progress[boundingApp];
				long double maximumYieldDenomLin = 1;
				long double maximumYieldDenomCnt = cTime-releaseDate[boundingApp];

    			// We identify the first application not receiving bandwidth
    			// Among the applications not receiving bandwidth it is one which achieves the minimum yield in the middle of the interval (and, thus, on any point in the interval)
    			bool nextAppFound = false;
    			int nextApp = -1;
    			long double yieldAtMiddle = 1.0;
    			long double deltaMiddle = (deltaMin+deltaMax)/2;
    			for(unsigned long i = 0; i<nApps; i++)
    			{
    				if((i!=constrainingApp)&&(!applisNeedingBW[i])){
    					if(nextAppFound){
    						long double tempYieldAtMiddle = fYield(i, deltaMiddle, cTime, 0, releaseDate, progress, maxBandwidth);
    						if (tempYieldAtMiddle < yieldAtMiddle){
    							yieldAtMiddle = tempYieldAtMiddle;
    							nextApp = i;
    						}
    					}
    					else{
    						yieldAtMiddle = fYield(i, deltaMiddle, cTime, 0, releaseDate, progress, maxBandwidth);
    						nextApp = i;
    						nextAppFound = true;
    					}
    				}
    			}

    			if ((ALLDEBUGNELEXMIN  || PECULIAR_CASE_DEBUG) && nextAppFound){
    				cout << "\t\t\tFirst application not receiving bandwidth " << nextApp << endl;
    			}

				bool intersectionFoundWithNextApp = false;
    			// We compute the coefficients of the yield function of this first application if it exists
    			if(nextAppFound){
    				long double nextAppYieldNumLin = 0;
    				long double nextAppYieldNumCnt = progress[nextApp];
    				long double nextAppYieldDenomLin = 1;
    				long double nextAppYieldDenomCnt = cTime-releaseDate[nextApp];

    				long double deltaIntersectionWithNextApp = deltaMax;

    				// We check whether this yield function intersects the ``common yield'' function
    				// except if we are in the singularity case and the maximum yield function is smaller
    				if((!inSingularityCase) || maxYieldIsLarger){
    					intersectionFoundWithNextApp = firstIntersection(deltaMin, deltaMax,
    							commonYieldNumLin, commonYieldNumCnt,
								commonYieldDenomLin, commonYieldDenomCnt,
								nextAppYieldNumLin, nextAppYieldNumCnt,
								nextAppYieldDenomLin, nextAppYieldDenomCnt,
								&deltaIntersectionWithNextApp);
    				}
    				else{
    					intersectionFoundWithNextApp = firstIntersection(deltaMin, deltaMax,
    							maximumYieldNumLin, maximumYieldNumCnt,
								maximumYieldDenomLin, maximumYieldDenomCnt,
								nextAppYieldNumLin, nextAppYieldNumCnt,
								nextAppYieldDenomLin, nextAppYieldDenomCnt,
								&deltaIntersectionWithNextApp);
    				}

    				if (intersectionFoundWithNextApp && ((!intersectionFound) ||(deltaIntersectionWithNextApp < lastValidDelta))){
    					firstIntersectionIsWithNextApp = true;
    					firstIntersectionIsWithConstrainingApp = false;
    					firstIntersectionIsWithBandwidthBound = false;
        				intersectionFound = true;
    					lastValidDelta = deltaIntersectionWithNextApp;
    				}

    				if(DEBUGNELEXMIN || PECULIAR_CASE_DEBUG){
    					if (intersectionFoundWithNextApp){
    						cout << "\t\tIntersection first application (application " << nextApp << ") and common yield found at time " << deltaIntersectionWithNextApp;
    						cout << " in [" << deltaMin << ", " << deltaMax << "]" << endl;

    						long double tmpCommonYield = (commonYieldNumLin*deltaIntersectionWithNextApp+commonYieldNumCnt)/(commonYieldDenomLin*deltaIntersectionWithNextApp+commonYieldDenomCnt);
    						cout <<"\t\tCommon yield at intersection: " << tmpCommonYield << endl;
    						long double tmpNextYield = fYield(nextApp, deltaIntersectionWithNextApp, cTime, 0, releaseDate, progress, maxBandwidth);
    						cout <<"\t\tYield next App at intersection: " << tmpNextYield << endl;
    						cout << "\t\tDifference: " << tmpCommonYield - tmpNextYield << endl;

    					}
    					else{
    						cout << "\t\tIntersection between first application (application " << nextApp << ") and common yield NOT found." << endl;
    					}
    				}
    			}

    			// We compute the coefficients of the yield function of the constraining application
    			long double constrainingAppYieldNumLin = 0;
    			long double constrainingAppYieldNumCnt = progress[constrainingApp]+remCommVol[constrainingApp]/maxBandwidth[constrainingApp];
    			long double constrainingAppYieldDenomLin = 1;
    			long double constrainingAppYieldDenomCnt = cTime-releaseDate[constrainingApp];

    			// We check whether this yield function intersects the ``common yield'' function
    			// except if we are in the singularity case and the maximum yield function is smaller
    			long double deltaIntersectWithConstrainingApp = deltaMax;
    			bool intersectionFoundWithConstrainingApp;
    			if((!inSingularityCase) || maxYieldIsLarger){
    				intersectionFoundWithConstrainingApp = firstIntersection(deltaMin, deltaMax,
    						commonYieldNumLin, commonYieldNumCnt,
							commonYieldDenomLin, commonYieldDenomCnt,
							constrainingAppYieldNumLin, constrainingAppYieldNumCnt,
							constrainingAppYieldDenomLin, constrainingAppYieldDenomCnt,
							&deltaIntersectWithConstrainingApp);
    			}
    			else{
    				intersectionFoundWithConstrainingApp = firstIntersection(deltaMin, deltaMax,
    						maximumYieldNumLin, maximumYieldNumCnt,
							maximumYieldDenomLin, maximumYieldDenomCnt,
							constrainingAppYieldNumLin, constrainingAppYieldNumCnt,
							constrainingAppYieldDenomLin, constrainingAppYieldDenomCnt,
							&deltaIntersectWithConstrainingApp);
    			}


    			if ((intersectionFoundWithConstrainingApp) &&
    					((!intersectionFound) || (deltaIntersectWithConstrainingApp < lastValidDelta+ERROR))){
    				firstIntersectionIsWithNextApp = false;
    				firstIntersectionIsWithConstrainingApp = true;
    				firstIntersectionIsWithBandwidthBound = false;
    				intersectionFound = true;

    				lastValidDelta = deltaIntersectWithConstrainingApp;
    			}

    			if(DEBUGNELEXMIN || PECULIAR_CASE_DEBUG){
    				cout << "\t\tIntersection between constraining application and common yield";
    				if (intersectionFoundWithConstrainingApp){
    					cout << " found at time " << deltaIntersectWithConstrainingApp << " in [" << deltaMin << ", " << deltaMax << "]" << endl;
    				}
    				else{
    					cout << " NOT found." << endl;
    				}
    			}

    			if (PECULIAR_CASE_DEBUG){
    				cout << "\tRIGHT before the potential conclusion" << endl;
    			}

    			if ((!intersectionFound) || (inSingularityCase && (!secondIntersectionFound) && (!intersectionFoundWithNextApp) && (!intersectionFoundWithConstrainingApp))){
        			if (PECULIAR_CASE_DEBUG){
        				cout << "\tIn the IF case" << endl;
        			}

    				*bestDelta = deltaMax;
    				computeYields(*bestDelta, cTime, constrainingApp,
    						nApps, remCommVol, releaseDate, progress, maxBandwidth,
							availableBandwidth, applisNeedingBW,
							allocatedBandwidth, yields);
    			}
    			else{
    				if (firstIntersectionIsWithConstrainingApp){
    					*bestDelta = lastValidDelta;
    					computeYields(*bestDelta, cTime, constrainingApp,
    							nApps, remCommVol, releaseDate, progress, maxBandwidth,
								availableBandwidth, applisNeedingBW,
								allocatedBandwidth, yields);
    				}
    				else{
    					// If the end of the validity interval is not distinguishable from the beginning of the interval
    					// the recursive call is likely to lead to an infinite loop...
    					if ((recurDepth > LIMITRECURSIVEDEPTH) && (lastValidDelta - deltaMin < SINGULARITY_PRECISION)){
    						bool HackActivation = true;
    						if(HackActivation && firstIntersectionIsWithBandwidthBound && ((remainingBandwidth < ERROR)||(remCommVol[nextApp]<0.001)) && (nbApplisNeedingBW==1) && applisNeedingBW[nextApp]){
    							long double newDelta = (remCommVol[constrainingApp]+remCommVol[nextApp])/availableBandwidth;
    							//    						*bestDelta = newDelta;
    							//    						computeYields(*bestDelta, cTime, constrainingApp,
    							//    								nApps, remCommVol, releaseDate, progress, maxBandwidth,
    							//									availableBandwidth, applisNeedingBW,
    							//									allocatedBandwidth, yields);
    							cout << " *** HACKFORPOTENTIAL BUG (1) *** " << endl;
    							MinYieldAtAppEventInInterval(newDelta, deltaMax, cTime,
    									constrainingApp,
										nApps, remCommVol, releaseDate, progress, maxBandwidth,
										availableBandwidth,
										allocatedBandwidth, yields, bestDelta, recurDepth+1);
    						}
    						else{
    							cout << HackActivation << " " << firstIntersectionIsWithBandwidthBound << " " << remainingBandwidth;
    							cout << " ERROR= " << ERROR << " " << remCommVol[nextApp] << " " << nbApplisNeedingBW << " " << applisNeedingBW[nextApp] << endl;

    							*bestDelta = deltaMin;
    							computeYields(*bestDelta, cTime, constrainingApp,
    									nApps, remCommVol, releaseDate, progress, maxBandwidth,
										availableBandwidth, applisNeedingBW,
										allocatedBandwidth, yields);
    							cerr << " *** POTENTIAL BUG (1) *** " << endl;
    							cout << "\t\t\t[" << deltaMin << ", " << deltaMax << "] : " << lastValidDelta << endl;
    							if(DEBUGNELEXMIN || PECULIAR_CASE_DEBUG){
    								cout << "\t\t\t\tThere should have been a recursive call (1)";
    								if (firstIntersectionIsWithBandwidthBound){
    									cout << " because first intersection is with bandwidth bound.";
    								}
    								if (firstIntersectionIsWithNextApp){
    									cout << " because first intersection is with next application";
    								}
    								cout << endl;
    							}
    							abort();
    							//    					MinYieldAtAppEventInInterval(deltaMin+ERROR, deltaMax, cTime,
    							//    							constrainingApp,
    							//								nApps, remCommVol, releaseDate, progress, maxBandwidth,
    							//								availableBandwidth,
    							//								allocatedBandwidth, yields, bestDelta);
    						}
    					}
    					else{
    						if(DEBUGNELEXMIN || PECULIAR_CASE_DEBUG){
    							cout << "\t\t\t\tRecursive call (1)";
    							if (firstIntersectionIsWithBandwidthBound){
    								cout << " because first intersection is with bandwidth bound.";
    							}
    							if (firstIntersectionIsWithNextApp){
    								cout << " because first intersection is with next application";
    							}
    							cout << endl;
    						}
    						MinYieldAtAppEventInInterval(lastValidDelta, deltaMax, cTime,
    								constrainingApp,
									nApps, remCommVol, releaseDate, progress, maxBandwidth,
									availableBandwidth,
									allocatedBandwidth, yields, bestDelta, recurDepth+1);
    					}
    				}
    			}
    		}
    	}
    }
    // The upper bound is achievable
    else{
    	if (DEBUGNELEXMIN || PECULIAR_CASE_DEBUG){
    		cout << "\t\tThe upper bound is achievable: required " << requiredBWForUpperBound << " available " << remainingBandwidth << endl;
    	}

    	// We identify the applications requiring bandwidth
		for(unsigned long i = 0; i<nApps; i++)
		{
			if(i!=constrainingApp){
				long double theYield = fYield(i, deltaMin, cTime, 0, releaseDate, progress, maxBandwidth);
				if (ALLDEBUGNELEXMIN || PECULIAR_CASE_DEBUG){
					cout << "\t\t\tApp " << i << " Y = " << theYield << " margin = " << theYield - upperBoundYield << " vs. " << upperBoundYield/INV_PRECISION << endl;
				}
				if (theYield - upperBoundYield < upperBoundYield/INV_PRECISION)
				{
					applisNeedingBW[i] = true;
				}
				else
				{
					applisNeedingBW[i] = false;
				}
			}
		}

    	if (DEBUGNELEXMIN || PECULIAR_CASE_DEBUG){
    		cout << "\t\tApplis needing bandwidth: ";
    		for(unsigned long i = 0; i<nApps; i++){
    			cout << i << ":" << applisNeedingBW[i] << "  ";
    		}
    		cout << endl;
    	}

		// We check whether the upper bound is achieved by the constraining application
		if (fYieldAtCompletion(constrainingApp, deltaMin, cTime, remCommVol, releaseDate, progress, maxBandwidth) < upperBoundYield+ERROR){
			*bestDelta = deltaMin;

			// The best solution is then defined at time deltaMin. As the whole bandwidth may not be distributed in a solution just achieving the minimum yield
			// we proceed by maximizing the minimum yield
			LexMinYield(*bestDelta, cTime,
					constrainingApp,
					nApps, remCommVol, releaseDate, progress, maxBandwidth,
					allocatedBandwidth, yields, availableBandwidth);
			if (DEBUGNELEXMIN || PECULIAR_CASE_DEBUG) cout << "\t\tUpper bound achieved by the constraining application." << endl;
		}
		else{
			// Identification of an application defining the upper-bound on the yield
			int boundingApp = -1;
			bool noAppFound = true;
			long double deltaMaxYield = -1;

			for(unsigned long i = 0; i<nApps; i++)
			{
				if ((i!=constrainingApp) && (fYieldMax(i, deltaMin, cTime, remCommVol, releaseDate, progress, maxBandwidth) < upperBoundYield+ERROR)) {
					long double tmpDeltaMaxYield = fYieldMax(i, deltaMax, cTime, remCommVol, releaseDate, progress, maxBandwidth);
					if ((noAppFound) || (tmpDeltaMaxYield < deltaMaxYield)){
						noAppFound = false;
						deltaMaxYield = tmpDeltaMaxYield;
						boundingApp = i;
					}
				}
			}

			// We must have found an application defining the upper bound
			assert(!noAppFound);
			assert(applisNeedingBW[boundingApp]);

			// This application must need bandwidth
			assert(applisNeedingBW[boundingApp]);
//			applisNeedingBW[boundingApp] = true;
//			// Note: we could have added an assert instead. The assignment is to be on the safe side (in case of rounding errors)

			if (DEBUGNELEXMIN || PECULIAR_CASE_DEBUG) cout << "\t\tUpper bound defined by application " << boundingApp << endl;

			// If the application defining the upper bound completes its communication at time deltaMin we can conclude
			if(requiredBandwidth(boundingApp, upperBoundYield, deltaMin, cTime, remCommVol, releaseDate, progress, maxBandwidth) > remCommVol[boundingApp]/deltaMin-ERROR){
				*bestDelta = deltaMin;

				// The best solution is then defined at time deltaMin. As the whole bandwidth may not be distributed in a solution just achieving the minimum yield
				// we proceed by maximizing the minimum yield
				if(PECULIAR_CASE_DEBUG)
				LexMinYield(*bestDelta, cTime,
						constrainingApp,
						nApps, remCommVol, releaseDate, progress, maxBandwidth,
						allocatedBandwidth, yields, availableBandwidth, true);
				else
					LexMinYield(*bestDelta, cTime,
							constrainingApp,
							nApps, remCommVol, releaseDate, progress, maxBandwidth,
							allocatedBandwidth, yields, availableBandwidth);
				if (DEBUGNELEXMIN || PECULIAR_CASE_DEBUG) cout << "\t\t\tThe application defining the upper bound completes its communication." << endl;
			}
			else{
				if (DEBUGNELEXMIN || PECULIAR_CASE_DEBUG) cout << "\t\t\tThe application defining the upper bound has a saturated bandwidth." << endl;

	    		// Data structure used to store the last date at which the solution is valid, and how the structure of the solution changes
	    		bool intersectionFound = false;
	    		long double lastValidDelta = deltaMax;
	    		bool firstIntersectionIsWithCommonYield = false;
	    		bool firstIntersectionIsWithConstrainingApp = false;
	    		bool firstIntersectionIsWithNextApp = false;

	    		// We compute the coefficients of the maximum yield function
				long double maximumYieldNumLin = 1;
				long double maximumYieldNumCnt = progress[boundingApp];
				long double maximumYieldDenomLin = 1;
				long double maximumYieldDenomCnt = cTime-releaseDate[boundingApp];

	    		// We compute the coefficient of the common yield function of the applications receiving bandwidth
	    		long double commonYieldNumLin = availableBandwidth;
	    		long double commonYieldNumCnt = -remCommVol[constrainingApp];
	    		long double commonYieldDenomLin = 0;
	    		long double commonYieldDenomCnt = 0;

	    		int nAppWithBW = 0;
	    		for(unsigned long i = 0; i<nApps; i++)
	    		{
	    			if (applisNeedingBW[i])
	    			{
	    				commonYieldNumCnt += progress[i]*maxBandwidth[i];
	    				commonYieldDenomLin += maxBandwidth[i];
	    				commonYieldDenomCnt += (cTime-releaseDate[i])*maxBandwidth[i];
	    				nAppWithBW++;
	    			}
	    		}

        		// We check whether this yield function intersects the ``maximum yield'' function
	    		// ONLY if there is more than 1 application requesting bandwidth: otherwise the two functions are identical !

	    		if (nAppWithBW>1){
	    			intersectionFound = firstIntersection(deltaMin, deltaMax,
	    					maximumYieldNumLin, maximumYieldNumCnt,
							maximumYieldDenomLin, maximumYieldDenomCnt,
							commonYieldNumLin, commonYieldNumCnt,
							commonYieldDenomLin, commonYieldDenomCnt,
							&lastValidDelta);
	    			if (intersectionFound){
	    				firstIntersectionIsWithCommonYield = true;
	    				firstIntersectionIsWithConstrainingApp = false;
	    				firstIntersectionIsWithNextApp = false;
	    			}
	    		}

        		if(DEBUGNELEXMIN || PECULIAR_CASE_DEBUG){
        			if(intersectionFound){
        				cout << "\t\tIntersection between common yield and maximum yield found at time " << lastValidDelta;
        				cout << " in [" << deltaMin << ", " << deltaMax << "] (offset of " << lastValidDelta-deltaMin  << ")" << endl;
        			}
        			else{
        				cout << "\t\tIntersection between common yield and maximum yield NOT found. " << endl;
        			}
        		}

        		// If the intersection of the two yield functions is at the beginning of the interval (or almost)
        		// then the beginning of the interval is a singularity time where the common yield of the applications needing bandwidth
        		// is equal to the yield of the application defining the upper-bound on the minimum yield.
        		// In that case, *inside* the interval the common yield may define the minimum yield

        		bool inSingularityCase = false;
        		bool solutionFound = false;
        		bool secondIntersectionFound = false;
				bool maxYieldIsLarger = false;
				long double deltaSecondIntersection = deltaMax;

				if (DEBUGNELEXMIN || PECULIAR_CASE_DEBUG){
				cout << "\t\t\tSingularity detection: " << lastValidDelta << " vs. " << deltaMin << " : " << lastValidDelta - deltaMin ;
				cout << " (with a precision limit of " << SINGULARITY_PRECISION << ") " << endl;
				}
        		if (intersectionFound && (lastValidDelta - deltaMin < SINGULARITY_PRECISION) && (lastValidDelta - deltaMin < deltaMin/INV_PRECISION)){
//            		if (intersectionFound && (lastValidDelta - deltaMin < deltaMin/INV_PRECISION)){
        			inSingularityCase = true;
        			if(DEBUGNELEXMIN || PECULIAR_CASE_DEBUG){
        				cout << "\t\t\tIn singularity case: ";
        				if(inSingularityCase) {
        					cout << " YES " << endl;
        				}
        				else{
        					cout << " no " << endl;
        				}
        			}
        			// If the common yield function is non-increasing, we can conclude
        			bool yieldIsIncreasing = isYieldIncreasing(cTime, constrainingApp,
        					nApps, remCommVol, releaseDate, progress, maxBandwidth,
							availableBandwidth, applisNeedingBW);
        			if(!yieldIsIncreasing){
        				*bestDelta = deltaMin;
        				computeYields(*bestDelta, cTime, constrainingApp,
        						nApps, remCommVol, releaseDate, progress, maxBandwidth,
								availableBandwidth, applisNeedingBW,
								allocatedBandwidth, yields);
        				solutionFound = true;
        			}
        			else{
        				// We check which function is the largest between the common-yield function and the maximum yield function
        				secondIntersectionFound = secondIntersection(deltaMin, deltaMax,
        						maximumYieldNumLin, maximumYieldNumCnt,
								maximumYieldDenomLin, maximumYieldDenomCnt,
								commonYieldNumLin, commonYieldNumCnt,
								commonYieldDenomLin, commonYieldDenomCnt,
								&deltaSecondIntersection, &maxYieldIsLarger);
        				if(DEBUGNELEXMIN || PECULIAR_CASE_DEBUG){
        				cout << "\t\t\tsecondIntersectionFound = " << secondIntersectionFound << endl;
        				cout << "\t\t\tat time : " << deltaSecondIntersection << endl;
        				cout << "\t\t\tmaxYieldIsLarger : " << maxYieldIsLarger << endl;
        				}
        				if(secondIntersectionFound){
        					intersectionFound = true;
        					lastValidDelta = deltaSecondIntersection;
        				}
        				else{
        					intersectionFound = false;
        					lastValidDelta = deltaMax;
    	    				firstIntersectionIsWithCommonYield = false;
    	    				firstIntersectionIsWithConstrainingApp = false;
    	    				firstIntersectionIsWithNextApp = false;
        				}
        			}
        		}

        		// If we are not in the singularity case or if no solution was found, we continue to look for the best solution
        		if (!solutionFound){
        			// Identification of an application not receiving any bandwidth and of minimum yield in that case
        			bool noNextApp = true;
        			int nextApp = -1;
        			long double minEmptyYield = 1;
        			long double deltaMiddle = (deltaMin+deltaMax)/2;

        			for(unsigned long i = 0; i<nApps; i++)
        			{
        				if ((!applisNeedingBW[i])&&(i!=constrainingApp)){
        					long double tmpMinEmptyYield = fYield(i, deltaMiddle, cTime, 0, releaseDate, progress, maxBandwidth);
        					if ((noNextApp) || (tmpMinEmptyYield < minEmptyYield)){
        						noNextApp = false;
        						nextApp = i;
        						minEmptyYield = tmpMinEmptyYield;
        					}
        				}
        			}

    				bool intersectionFoundWithMax = false;
        			// If such an application was found
        			if(!noNextApp){
        				long double nextAppYieldNumLin = 0;
        				long double nextAppYieldNumCnt = progress[nextApp];
        				long double nextAppYieldDenomLin = 1;
        				long double nextAppYieldDenomCnt = cTime-releaseDate[nextApp];
        				long double deltaIntersectionWithNextApp = deltaMax;

        				// We check whether this yield function intersects the ``maximum yield'' function,
        				// except if we are in the singularity case and the common yield function is smaller
        				if(inSingularityCase && maxYieldIsLarger){
            				intersectionFoundWithMax = firstIntersection(deltaMin, deltaMax,
    								commonYieldNumLin, commonYieldNumCnt,
    								commonYieldDenomLin, commonYieldDenomCnt,
    								nextAppYieldNumLin, nextAppYieldNumCnt,
    								nextAppYieldDenomLin, nextAppYieldDenomCnt,
    								&deltaIntersectionWithNextApp);

        				}
        				else{
        					intersectionFoundWithMax = firstIntersection(deltaMin, deltaMax,
        							maximumYieldNumLin, maximumYieldNumCnt,
									maximumYieldDenomLin, maximumYieldDenomCnt,
									nextAppYieldNumLin, nextAppYieldNumCnt,
									nextAppYieldDenomLin, nextAppYieldDenomCnt,
									&deltaIntersectionWithNextApp);
        				}


        				if ((intersectionFoundWithMax) &&
        						((!intersectionFound) || (deltaIntersectionWithNextApp < lastValidDelta))){
        					firstIntersectionIsWithCommonYield = false;
        					firstIntersectionIsWithNextApp = true;
        					firstIntersectionIsWithConstrainingApp = false;
        					intersectionFound = true;
        					lastValidDelta = deltaIntersectionWithNextApp;
        				}
        				if(DEBUGNELEXMIN || PECULIAR_CASE_DEBUG){
        					if(intersectionFoundWithMax){
        						cout << "\t\tIntersection between next application (application " << nextApp << ") and maximum yield found at time " << lastValidDelta;
        						cout << " in [" << deltaMin << ", " << deltaMax << "]" << endl;
        					}
        					else{
        						cout << "\t\tIntersection between next application and maximum yield NOT found." << endl;
        					}
        				}
        			}
        			else{
        				if(DEBUGNELEXMIN){
        					cout << "\t\tNo intersection because no next application." << endl;
        				}
        			}



        			// We compute the coefficients of the yield function of the constraining application
        			long double constrainingAppYieldNumLin = 0;
        			long double constrainingAppYieldNumCnt = progress[constrainingApp]+remCommVol[constrainingApp]/maxBandwidth[constrainingApp];
        			long double constrainingAppYieldDenomLin = 1;
        			long double constrainingAppYieldDenomCnt = cTime-releaseDate[constrainingApp];

    				// We check whether this yield function intersects the ``maximum yield'' function,
    				// except if we are in the singularity case and the common yield function is smaller
        			long double deltaIntersectWithConstrainingApp = deltaMax;
        			bool intersectionFoundWithConstrainingApp = false;
        			if(inSingularityCase && maxYieldIsLarger){
        				intersectionFoundWithConstrainingApp = firstIntersection(deltaMin, deltaMax,
        						commonYieldNumLin, commonYieldNumCnt,
								commonYieldDenomLin, commonYieldDenomCnt,
								constrainingAppYieldNumLin, constrainingAppYieldNumCnt,
								constrainingAppYieldDenomLin, constrainingAppYieldDenomCnt,
								&deltaIntersectWithConstrainingApp);

        			}
        			else{
        				intersectionFoundWithConstrainingApp = firstIntersection(deltaMin, deltaMax,
        						maximumYieldNumLin, maximumYieldNumCnt,
								maximumYieldDenomLin, maximumYieldDenomCnt,
								constrainingAppYieldNumLin, constrainingAppYieldNumCnt,
								constrainingAppYieldDenomLin, constrainingAppYieldDenomCnt,
								&deltaIntersectWithConstrainingApp);
        			}

        			if ((intersectionFoundWithConstrainingApp) &&
        					((!intersectionFound) || (deltaIntersectWithConstrainingApp < lastValidDelta+ERROR))){
        				firstIntersectionIsWithCommonYield = false;
        				firstIntersectionIsWithNextApp = false;
        				firstIntersectionIsWithConstrainingApp = true;
        				intersectionFound = true;

        				lastValidDelta = deltaIntersectWithConstrainingApp;
        			}

        			if(DEBUGNELEXMIN || PECULIAR_CASE_DEBUG){
        				if(intersectionFoundWithConstrainingApp){
        					cout << "\t\tIntersection between constraining application and maximum yield found at time " << deltaIntersectWithConstrainingApp;
        					cout << " in [" << deltaMin << ", " << deltaMax << "]" << endl;
        				}
        				else{
        					cout << "\t\tIntersection between constraining application and maximum yield NOT found." << endl;
        				}
        			}


        			if ((!intersectionFound) || (inSingularityCase && (!secondIntersectionFound) && (!intersectionFoundWithMax) && (!intersectionFoundWithConstrainingApp))){
        				*bestDelta = deltaMax;
        				//	    			computeYields(*bestDelta, cTime, constrainingApp,
        				//	    					nApps, remCommVol, releaseDate, progress, maxBandwidth,
        				//							availableBandwidth, applisNeedingBW,
        				//							allocatedBandwidth, yields);
        				// The best solution is then defined at time deltaMax. As the whole bandwidth may not be distributed in a solution just achieving the minimum yield
        				// we proceed by maximizing the minimum yield
        				LexMinYield(*bestDelta, cTime,
        						constrainingApp,
								nApps, remCommVol, releaseDate, progress, maxBandwidth,
								allocatedBandwidth, yields, availableBandwidth);
        			}
        			else{
        				if (firstIntersectionIsWithConstrainingApp){
        					*bestDelta = lastValidDelta;
        					//	    				computeYields(*bestDelta, cTime, constrainingApp,
        					//	    						nApps, remCommVol, releaseDate, progress, maxBandwidth,
        					//								availableBandwidth, applisNeedingBW,
        					//								allocatedBandwidth, yields);
        					// The best solution is then defined at time lastValidDelta. As the whole bandwidth may not be distributed in a solution just achieving the minimum yield
        					// we proceed by maximizing the minimum yield
        					LexMinYield(*bestDelta, cTime,
        							constrainingApp,
									nApps, remCommVol, releaseDate, progress, maxBandwidth,
									allocatedBandwidth, yields, availableBandwidth);	    			}
        				else{
        					// If the end of the validity interval is not distinguishable from the beginning of the interval
        					// the recursive call is likely to lead to an infinite loop...

        					//	    				if(HackActivation && firstIntersectionIsWithBandwidthBound && (remainingBandwidth < ERROR) && (nbApplisNeedingBW==1) && applisNeedingBW[nextApp]){
        					//	    				    						long double newDelta = (remCommVol[constrainingApp]+remCommVol[nextApp])/availableBandwidth;
        					//	    				//    						*bestDelta = newDelta;
        					//	    				//    						computeYields(*bestDelta, cTime, constrainingApp,
        					//	    				//    								nApps, remCommVol, releaseDate, progress, maxBandwidth,
        					//	    				//									availableBandwidth, applisNeedingBW,
        					//	    				//									allocatedBandwidth, yields);
        					//	    				    						cerr << " *** HACKFORPOTENTIAL BUG (1) *** " << endl;
        					//	    				    						MinYieldAtAppEventInInterval(newDelta, deltaMax, cTime,
        					//	    				    								constrainingApp,
        					//	    													nApps, remCommVol, releaseDate, progress, maxBandwidth,
        					//	    													availableBandwidth,
        					//	    													allocatedBandwidth, yields, bestDelta, recurDepth+1);
        					//	    				    					}


        					if ((recurDepth>LIMITRECURSIVEDEPTH) && (lastValidDelta - deltaMin < SINGULARITY_PRECISION)){
        						if (remCommVol[boundingApp] < 0.001){
        							if(maxBandwidth[constrainingApp]+maxBandwidth[boundingApp] >= availableBandwidth){
        								*bestDelta = (remCommVol[constrainingApp]+remCommVol[boundingApp])/availableBandwidth;
        								allocatedBandwidth[constrainingApp] = remCommVol[constrainingApp]/(*bestDelta);
        								allocatedBandwidth[boundingApp] = remCommVol[boundingApp]/(*bestDelta);
        							}
        							else{
        								double newDelta = deltaMin+remCommVol[boundingApp]/maxBandwidth[boundingApp];
        								*bestDelta = newDelta;
        								allocatedBandwidth[constrainingApp] = remCommVol[constrainingApp]/(*bestDelta);
        								allocatedBandwidth[boundingApp] = maxBandwidth[boundingApp];
        								remainingBandwidth = availableBandwidth-allocatedBandwidth[constrainingApp]-allocatedBandwidth[boundingApp];
        								double totalRemainingRequestedBW = 0;
        								for(unsigned long i = 0; i<nApps; i++){
        									if ((i!=constrainingApp)&&(i!=boundingApp)){
        										totalRemainingRequestedBW += maxBandwidth[i];
        									}
        								}
        								for(unsigned long i = 0; i<nApps; i++){
        									if ((i!=constrainingApp)&&(i!=boundingApp)){
        										allocatedBandwidth[i] = maxBandwidth[i]*remainingBandwidth/totalRemainingRequestedBW;
        									}
        								}
        							}
        							cout << "***PECULIAR LIMIT behavior " << endl;
        							LexMinYield(*bestDelta, cTime,
        									constrainingApp,
											nApps, remCommVol, releaseDate, progress, maxBandwidth,
											allocatedBandwidth, yields, availableBandwidth);
        						}
        						else{
        							if((recurDepth<=LIMITRECURSIVEDEPTH) && (recurDepth>LIMITRECURSIVEDEPTH/2)){
        								MinYieldAtAppEventInInterval(lastValidDelta+SINGULARITY_PRECISION, deltaMax, cTime,
        										constrainingApp,
												nApps, remCommVol, releaseDate, progress, maxBandwidth,
												availableBandwidth,
												allocatedBandwidth, yields, bestDelta, recurDepth+1);
        							}
        							else{


        								*bestDelta = deltaMin;
        								//	    					computeYields(*bestDelta, cTime, constrainingApp,
        								//	    							nApps, remCommVol, releaseDate, progress, maxBandwidth,
        								//									availableBandwidth, applisNeedingBW,
        								//									allocatedBandwidth, yields);
        								// The best solution is then defined at time deltaMin. As the whole bandwidth may not be distributed in a solution just achieving the minimum yield
        								// we proceed by maximizing the minimum yield
        								LexMinYield(*bestDelta, cTime,
        										constrainingApp,
												nApps, remCommVol, releaseDate, progress, maxBandwidth,
												allocatedBandwidth, yields, availableBandwidth);
        								cerr << " *** POTENTIAL BUG (2) *** " << endl;
        								if(DEBUGNELEXMIN || PECULIAR_CASE_DEBUG){
        									cout << "\t\t\t\tThere should have been a recursive call (2)... ";
        									if (firstIntersectionIsWithCommonYield){
        										cout << " because first intersection is with the total bandwidth bound.";
        									}
        									if (firstIntersectionIsWithNextApp){
        										cout << " because first intersection is with next application, namely " << nextApp;
        									}
        									cout << endl;
        								}

        								cout << "\t\t\t\tThe intersection was found at time " << lastValidDelta << " in [" << deltaMin << ", " << deltaMax << "]" << endl;


        								abort();
        							}
        						}
        					}
        					else{
        						if(DEBUGNELEXMIN || PECULIAR_CASE_DEBUG){
        							cout << "\t\t\t\tRecursive call (2)";
        							if (firstIntersectionIsWithCommonYield){
        								cout << " because first intersection is with the total bandwidth bound.";
        							}
        							if (firstIntersectionIsWithNextApp){
        								cout << " because first intersection is with next application, namely " << nextApp;
        							}
        							cout << endl;
        						}
        						MinYieldAtAppEventInInterval(lastValidDelta, deltaMax, cTime,
        								constrainingApp,
										nApps, remCommVol, releaseDate, progress, maxBandwidth,
										availableBandwidth,
										allocatedBandwidth, yields, bestDelta, recurDepth+1);
        					}
        				}
        			}
        		}
			}
		}
    }



//
//
//
//		if (ALLDEBUGNELEXMIN){
//			cout << "indexMin = " << indexMin << " indexMax = " << indexMax << " Remaining bandwidth: " << remainingBandwidth << endl;
//			for(unsigned long i = 0; i<nApps; i++)
//			{
//				cout << "\t Appli " << i << " needing BW : " << applisNeedingBW[i] << endl;
//				cout << "\t\t Yield at " << deltaMin << " WO BW " << fYield(i, deltaMin, cTime, 0, releaseDate, progress, maxBandwidth) << endl;
//			}
//
//			long double yieldConstrainingAppAtDeltaMin = fYield(constrainingApp, deltaMin, cTime, remCommVol[constrainingApp]/deltaMin, releaseDate, progress, maxBandwidth);
//
//			cout << "\t Yield of constraining app at " << deltaMin << " = " << yieldConstrainingAppAtDeltaMin << endl;
//		}
//
//
//		bool optimalFound = false;
//		int firstAppNotIncluded = get<1>(yieldsWOBandwidth[indexMax]);
//		int indexFirstAppNotIncluded = indexMax;
//		delta = deltaMin;
//		while(!optimalFound){
//			// We check whether the yield is increasing
//			long double yieldIsIncreasing = isYieldIncreasing(cTime, constrainingApp,
//					nApps, remCommVol, releaseDate, progress, maxBandwidth,
//					availableBandwidth, applisNeedingBW);
//			// If the yield is decreasing the optimal has been found and we are done
//			if(!yieldIsIncreasing){
//				optimalFound = true;
//				if (ALLDEBUGNELEXMIN) cout << "\t Optimal found : yield is decreasing" << endl;
//			}
//			// Otherwise the yield is increasing
//			else{
//				// We check whether at the end of the interval the yield of this solution is still smaller than or equal
//				// to the yield of the first application to which bandwidth has not been allocated
//				// Computation of the optimal yield
//				long double sumInNumerator = 0;
//				long double sumInDenominator = 0;
//				for(unsigned long i = 0; i<nApps; i++)
//				{
//					if ((applisNeedingBW[i])&&(i!=constrainingApp))
//					{
//						sumInNumerator += progress[i]*maxBandwidth[i];
//						sumInDenominator += (cTime+deltaMax-releaseDate[i])*maxBandwidth[i];
//					}
//				}
//				long double yieldAtDeltaMax = (deltaMax*availableBandwidth-remCommVol[constrainingApp]+sumInNumerator)/sumInDenominator;
//				// If the yield of the set of application is smaller at deltaMax we have found the optimum
//				if ((indexFirstAppNotIncluded==yieldsWOBandwidth.size())
//						|| (yieldAtDeltaMax - fYield(firstAppNotIncluded, deltaMax, cTime, 0, releaseDate, progress, maxBandwidth) <= yieldAtDeltaMax/INV_PRECISION)){
//					delta = deltaMax;
//					optimalFound = true;
//					if (ALLDEBUGNELEXMIN) cout << "\t Optimal found : delta = deltaMax" << endl;
//				}
//				// Otherwise the two functions intersects
//				else{
//					// We compute the intersection point
//					// We compute the coefficients of the solution defined by the applications receiving bandwidth (so far)
//					long double numConstant = -remCommVol[constrainingApp];
//					long double denomConstant = 0;
//					long double denomDeltaFactor = 0;
//					for(unsigned long i = 0; i<nApps; i++)
//					{
//						if(applisNeedingBW[i] && (i!=constrainingApp)){
//							numConstant += maxBandwidth[i]*progress[i];
//							denomConstant += maxBandwidth[i]*(cTime-releaseDate[i]);
//							denomDeltaFactor += maxBandwidth[i];
//						}
//					}
//					// We solve the second degree polynomial of the intersection
//					long double aSecondDegreePolynom = availableBandwidth;
//					long double bSecondDegreePolynom = availableBandwidth*(cTime-releaseDate[firstAppNotIncluded])+numConstant-denomDeltaFactor*progress[firstAppNotIncluded];
//					long double cSecondDegreePolynom = numConstant*(cTime-releaseDate[firstAppNotIncluded])-progress[firstAppNotIncluded]*denomConstant;
//					long double discriminantSQ = bSecondDegreePolynom*bSecondDegreePolynom-4*aSecondDegreePolynom*cSecondDegreePolynom;
//					assert(discriminantSQ >= -ERROR);
//					if (discriminantSQ <= ERROR){
//						discriminantSQ = 0;
//					}
//					long double discriminant = sqrt(discriminantSQ);
//					long double firstRoot = (-bSecondDegreePolynom-discriminant)/(2*aSecondDegreePolynom);
//					long double secondRoot = (-bSecondDegreePolynom+discriminant)/(2*aSecondDegreePolynom);
//					// We want to identify the first solution belonging in the interval [delta, deltaMax] (which exists by definition)
//					// We first sort the roots
//					if (secondRoot < firstRoot){
//						long double tmp = firstRoot;
//						firstRoot = secondRoot;
//						secondRoot = tmp;
//					}
//					long double newDelta = 0;
//
//					if (ALLDEBUGNELEXMIN){
//						cout << "deltaMin= " << deltaMin << "\tdelta= " << delta << "\tdeltaMax = " << deltaMax << endl;
//						cout << "firstRoot = " << firstRoot << "\tsecondRoot = " << secondRoot << endl;
//					}
//
//					assert(firstRoot -deltaMax <= deltaMax/INV_PRECISION); // Otherwise both roots are strictly greater than deltaMax and there is no solution; hence, a contradiction
//					assert(delta - secondRoot <= delta/INV_PRECISION); // Otherwise both roots are strictly smaller than delta and there is no solution; hence, a contradiction
//					if ((delta - firstRoot < delta/INV_PRECISION) && (firstRoot - deltaMax< deltaMax/INV_PRECISION) ){
////						cout << "\tIn the first case " << endl;
//						newDelta = firstRoot;
//					}
//					else{
//						if ((delta - secondRoot < delta/INV_PRECISION) && (secondRoot - deltaMax < deltaMax/INV_PRECISION)){
//							newDelta = secondRoot;
////							cout << "\tIn the second case " << endl;
//						}
//						else{
//							std::cerr << "***Impossible case***\nNo root in the interval [" << delta << ", " << deltaMax << "]";
//							std::cerr << "\tCandidate roots were " << firstRoot << " and " << secondRoot << std::endl;
//							abort();
//						}
//					}
//					if (newDelta - deltaMax >= - deltaMax/INV_PRECISION){
//						newDelta = deltaMax;
//						optimalFound = true;
//					}
//
//					delta = newDelta;
//
//					if(!optimalFound){
//						// We enroll all applications corresponding to the next step
//						long double limitYield = fYield(firstAppNotIncluded, newDelta, cTime, 0, releaseDate, progress, maxBandwidth);
////						cout << "limitYield = " << limitYield << endl;
//						do{
////							cout << "\tIncluding app " << firstAppNotIncluded << endl;
//							applisNeedingBW[firstAppNotIncluded] = true;
//							indexFirstAppNotIncluded++;
//							if (indexFirstAppNotIncluded < yieldsWOBandwidth.size()){
//								firstAppNotIncluded = get<1>(yieldsWOBandwidth[indexFirstAppNotIncluded]);
//							}
//						}
//						while( (indexFirstAppNotIncluded<yieldsWOBandwidth.size())
//								&&(fYield(firstAppNotIncluded, newDelta, cTime, 0, releaseDate, progress, maxBandwidth) < limitYield+ERROR));
////						if (indexFirstAppNotIncluded<yieldsWOBandwidth.size()){
////							cout << "fYield = " << fYield(firstAppNotIncluded, newDelta, cTime, 0, releaseDate, progress, maxBandwidth) << endl;
////						}
//					}
//				}
//			}
//		}
//
//	}
////
////	// If the yield is decreasing along the interval the best solution is found at deltaMin
////	if (trendIndicator>=0){
////		delta = deltaMin;
////	}
////	// Otherwise the yield is (at least initially) an increasing function of delta
////	else{
////	// If all applications must receive some bandwidth the optimal solution is known
////		if(indexMin==indexMax){
////			delta=deltaMax;
////		}
////		else{
////			// We must check whether
////
////			// Yield at deltaMax of applications receiving bandwidth
////			long double yieldAtDeltaMax = ;
////			if
////		}
////	}
//
//	// Sanity check
//	applisNeedingBW[constrainingApp] = false;
//	allocatedBandwidth[constrainingApp] = remCommVol[constrainingApp]/delta;
//	remainingBandwidth = availableBandwidth - allocatedBandwidth[constrainingApp];
//
//	// Computation of the optimal yield
//	long double sumInNumerator = 0;
//	long double sumInDenominator = 0;
//	for(unsigned long i = 0; i<nApps; i++)
//	{
//		if (applisNeedingBW[i])
//		{
//			sumInNumerator += progress[i]*maxBandwidth[i];
//			sumInDenominator += (cTime+delta-releaseDate[i])*maxBandwidth[i];
//		}
//	}
//	long double yield = (delta*remainingBandwidth+sumInNumerator)/sumInDenominator;
//
//	if (ALLDEBUGNELEXMIN)
//	{
//		cout << "\t\tYield = " << yield << "\t at time " << delta << endl;
//		cout << "\t\tYield of constraining app = " << fYield(constrainingApp, delta, cTime, allocatedBandwidth[constrainingApp], releaseDate, progress, maxBandwidth) << endl;
//	}
//
//	int nbAllocatedApps = 0;
//	int lastAllocatedApp = -1;
//	long double totalBWAllocated = 0;
//	for(unsigned long i = 0; i<nApps; i++)
//	{
//		if (ALLDEBUGNELEXMIN)
//		{
//			cout << "About application: " << i << " yield WO BW = " << fYield(i, delta, cTime, 0, releaseDate, progress, maxBandwidth) << endl;
//			cout << "\t yield with all the bandwidth = " << fYield(i, delta, cTime, remainingBandwidth, releaseDate, progress, maxBandwidth) << endl;
//		}
//		if ((remainingBandwidth > ERROR) && (applisNeedingBW[i]))
//		{
//			if (ALLDEBUGNELEXMIN)
//			cout << "\t need BW" << endl;
//			allocatedBandwidth[i] = requiredBandwidth(i, yield, delta, cTime, remCommVol, releaseDate, progress, maxBandwidth);
//			totalBWAllocated += allocatedBandwidth[i];
//			if (ALLDEBUGNELEXMIN)
//				cout << "\t allocated: " << allocatedBandwidth[i] << endl;
//			yields[i] = yield;
//			nbAllocatedApps++;
//			lastAllocatedApp = i;
//		}
//		else{
//			if (ALLDEBUGNELEXMIN)
//				cout << "\t does not need BW" << endl;
//			yields[i] = fYield(i, delta, cTime, allocatedBandwidth[i], releaseDate, progress, maxBandwidth);
//		}
//	}
//
//	if (ALLDEBUGNELEXMIN)
//	{
//		cout << "nbAllocatedApps = " << nbAllocatedApps << " totalBWAllocated=  " << totalBWAllocated << " for " << remainingBandwidth << endl;
//	}
//
//	if(remainingBandwidth>availableBandwidth/INV_PRECISION){
//		assert(nbAllocatedApps>0);
//	if (nbAllocatedApps==1){
//		allocatedBandwidth[lastAllocatedApp] = remainingBandwidth;
//	}
//	else{
//		if (totalBWAllocated==0){
//			for(unsigned long i = 0; i<nApps; i++){
//				if (applisNeedingBW[i])
//				{
//					allocatedBandwidth[i] = remainingBandwidth/nbAllocatedApps;
//				}
//			}
//		}
//		else{
//			if (totalBWAllocated > remainingBandwidth){
//				long double correctionRatio = remainingBandwidth/totalBWAllocated;
//				for(unsigned long i = 0; i<nApps; i++){
//					if (applisNeedingBW[i])
//					{
//						allocatedBandwidth[i] = allocatedBandwidth[i]*correctionRatio;
////						cout << " Allocating to " << i << " BW " << allocatedBandwidth[i] << endl;
//						yields[i] = fYield(i, delta, cTime, allocatedBandwidth[i], releaseDate, progress, maxBandwidth);
//					}
//				}
//			}
//			else{
//				if ((remainingBandwidth > ERROR) && (remainingBandwidth-totalBWAllocated > remainingBandwidth/INV_PRECISION)){
//					long double correctionRatio = remainingBandwidth/totalBWAllocated;
//					for(unsigned long i = 0; i<nApps; i++){
//						if (applisNeedingBW[i])
//						{
//							allocatedBandwidth[i] = allocatedBandwidth[i]*correctionRatio;
////							cout << " Allocating to " << i << " BW " << allocatedBandwidth[i] << endl;
//							yields[i] = fYield(i, delta, cTime, allocatedBandwidth[i], releaseDate, progress, maxBandwidth);
//						}
//					}
//				}
////			if (remainingBandwidth < ERROR){
////				cout << " ALMOST no remaining bandwidth" << endl;
//			}
//		}
//	}
//	}
//
//	if (ALLDEBUGNELEXMIN)
//		cout << "Everything is OK..." << endl;


    if (PECULIAR_CASE_DEBUG){
    	cout << "**END of MinYieldAtAppEventInInterval" << endl;
    	for(unsigned long i = 0; i<nApps; i++){
    		cout << "\tApplication " << i << " BW " << allocatedBandwidth[i] << " (" << maxBandwidth[i] << ")" << endl;
    	}
    }
}

void MinYieldAtAppEventInInterval(double deltaMin, double deltaMax, double cTime,
		int constrainingApp,
		unsigned long nApps, double * remCommVol, double * releaseDate, double * progress, double * maxBandwidth,
		double availableBandwidth,
		double * allocatedBandwidth, double * yields, double * bestDelta){

	MinYieldAtAppEventInInterval(deltaMin, deltaMax, cTime,
			constrainingApp,
			nApps, remCommVol, releaseDate, progress, maxBandwidth,
			availableBandwidth,
			allocatedBandwidth, yields, bestDelta, 0);

	// We sort the array tempYields
	sort(yields, yields + nApps);
}

// Computes the best mean yield that can be achieved between the dates cTime+firstBound and cTime+secondBound if application ``constraining''
// is the one defining the event
void MinYieldAtAppEventBUGGED(double firstBound, double secondBound, double cTime,
		int constrainingApp,
		unsigned long nApps, double * remCommVol, double * releaseDate, double * progress, double * maxBandwidth,
		double availableBandwidth,
		double * allocatedBandwidth, double * yields){

	if (ALLDEBUGNELEXMIN) cout << "\t\t\t\tMinYieldAtAppEvent computation in [ " << firstBound << ", " << secondBound << "] avail bw= " << availableBandwidth  << endl;

	if(ALLDEBUGNELEXMIN){
		for(unsigned long i = 0; i<nApps; i++){
			cout <<"\t\t\t\t\tAppli " << i << " max bw: " << maxBandwidth[i] << "\tRemaining vol: " << remCommVol[i] << endl;
		}
	}



	// We initialize the data structure which will hold the data allocation
	for(int i=0; i<nApps; i++){
		allocatedBandwidth[i] = 0;
	}
	double delta = (firstBound+secondBound)/2;
	allocatedBandwidth[constrainingApp] = remCommVol[constrainingApp]/delta;
	cout << "Initial allocation to constraining application (" << constrainingApp << "): " << allocatedBandwidth[constrainingApp] << endl;

	double remainingBandwidth = availableBandwidth - allocatedBandwidth[constrainingApp];

	// Build and sort the set of yields achieved when no bandwidth is allocated
	vector<long double> yieldsWOBandwidth;
	for(unsigned long i = 0; i<nApps; i++)
	{
		if (i != constrainingApp) yieldsWOBandwidth.push_back(fYield(i, delta, cTime, 0, releaseDate, progress, maxBandwidth));
	}
    sort(yieldsWOBandwidth.begin(), yieldsWOBandwidth.end());

	// Identification of applications which need bandwidth to reach the target yield
	bool applisNeedingBW[nApps];
    // Identification of the interval of ``yields without any bandwidth allocated'' including the optimal solution
    unsigned long indexMin = 0;
    unsigned long indexMax = yieldsWOBandwidth.size()-1;

    cout << "before if" << endl;

	if (totalRequiredBandwidth(nApps, yieldsWOBandwidth[indexMax], delta, cTime, remCommVol, releaseDate, progress, maxBandwidth, constrainingApp) <= remainingBandwidth+ERROR)
	{
	    cout << "in then " << endl;


		indexMin = indexMax;
		for(unsigned long i = 0; i<nApps; i++)
		{
			applisNeedingBW[i] = true;
		}
		applisNeedingBW[constrainingApp] = false;
	}
	else
	{
	    cout << "in else " << endl;

		cout << "\tindexMin= " << indexMin << "\tindexMax= " << indexMax << endl;

		while (indexMax-indexMin>1)
		{

			cout << "indexMin= " << indexMin << "\tindexMax= " << indexMax << endl;

			unsigned long indexMiddle = (indexMin+indexMax)/2;
			if (totalRequiredBandwidth(nApps, yieldsWOBandwidth[indexMiddle], delta, cTime, remCommVol, releaseDate, progress, maxBandwidth,constrainingApp) > remainingBandwidth-ERROR)
			{
				indexMax = indexMiddle;
			}
			else
			{
				indexMin = indexMiddle;
			}
		}
		for(unsigned long i = 0; i<nApps; i++)
		{
			if(i!=constrainingApp){
				if (fYield(i, delta, cTime, 0, releaseDate, progress, maxBandwidth) < yieldsWOBandwidth[indexMin]+ERROR)
				{
					applisNeedingBW[i] = true;
				}
				else
				{
					applisNeedingBW[i] = false;
				}
			}
		}
	}

	// Sanity check
	applisNeedingBW[constrainingApp] = false;

	// To identify the maximum yield in the interval we study the yield variation
	long double sumMaxBandwidths = 0;
	long double sumPseudoCommVol = 0;
	long double sumNegativePart = 0;
	for(int i=0; i<nApps; i++){
		if (applisNeedingBW[i]){
			sumMaxBandwidths += maxBandwidth[i];
			sumPseudoCommVol += maxBandwidth[i]*progress[i];
			sumNegativePart += maxBandwidth[i]*(availableBandwidth*(cTime-releaseDate[i])+remCommVol[constrainingApp]);
		}
	}

	long double bestDelta = 0;
	long double numerator = sumMaxBandwidths*sumPseudoCommVol-sumNegativePart;
	if (numerator <0 ){
		bestDelta = secondBound;
	}
	else{
		bestDelta = firstBound;
	}

	allocatedBandwidth[constrainingApp] = remCommVol[constrainingApp]/bestDelta;
	remainingBandwidth = availableBandwidth-allocatedBandwidth[constrainingApp];
	if (remainingBandwidth < ERROR){
		remainingBandwidth = 0;
	}

	cout << "*&*&*&* allocatedBandwidth[" << constrainingApp << "] = " << allocatedBandwidth[constrainingApp] << " for time: " << bestDelta << " in [" <<firstBound << ", " << secondBound << "]"  << endl;
	cout << "\tnumerator was " << numerator << endl;
	cout << sumMaxBandwidths << " " << sumPseudoCommVol << " " << sumNegativePart << endl;

	// Computation of the optimal yield
	long double sumInNumerator = 0;
	long double sumInDenominator = 0;
	for(unsigned long i = 0; i<nApps; i++)
	{
		if (applisNeedingBW[i])
		{
			sumInNumerator += progress[i]*maxBandwidth[i];
			sumInDenominator += (cTime+bestDelta-releaseDate[i])*maxBandwidth[i];
		}
	}
	long double yield = (delta*remainingBandwidth+sumInNumerator)/sumInDenominator;

cout << "\t\tYield = " << yield << endl;

	for(unsigned long i = 0; i<nApps; i++)
	{
		cout << "About application: " << i << endl;
		if (applisNeedingBW[i])
		{
			cout << "\t need BW" << endl;
			allocatedBandwidth[i] = requiredBandwidth(i, yield, delta, cTime, remCommVol, releaseDate, progress, maxBandwidth);
			if(allocatedBandwidth[i] > maxBandwidth[i]){
				allocatedBandwidth[i] = maxBandwidth[i];
			}
			cout << "\t allocated: " << allocatedBandwidth[i] << endl;
			yields[i] = yield;
		}
		else{
			cout << "\t does not need BW" << endl;
			yields[i] = fYield(i, delta, cTime, allocatedBandwidth[i], releaseDate, progress, maxBandwidth);
		}
	}
	cout << "Everything is OK..." << endl;
}


// Computes the best mean yield that can be achieved between the dates cTime+firstBound and cTime+secondBound
// knowing that not all applications can have a yield greater than or equal to the constraining application
void MinYieldAtEvent(double firstBound, double secondBound, double cTime,
		unsigned long nApps, double * remCommVol, double * releaseDate, double * progress, double * maxBandwidth,
		bool * processedApplis,
		double availableBandwidth,
		double * allocatedBandwidth, double * yields){

	// Data structure to hold temporary results
	double tmpAllocatedBandwidth[nApps];
	double tmpYields[nApps];
	// Initialization of results to be produced
	for(int i=0; i<nApps; i++){
		allocatedBandwidth[i] = 0;
		yields[i] = 0;
	}

	// We loop over all the (non already processed) applications as possible constraining application
	// if their maximum admissible bandwidth would enable them to complete in the considered time window
	for(int i=0; i<nApps; i++){
		if((!processedApplis[i])&&(firstBound*maxBandwidth[i]>=remCommVol[i]-ERROR)){
			MinYieldAtAppEvent(firstBound, cTime,
					i,
					nApps, remCommVol, releaseDate, progress, maxBandwidth,
					availableBandwidth,
					tmpAllocatedBandwidth, tmpYields);

			keepBestSolution(nApps, tmpYields, tmpAllocatedBandwidth, yields, allocatedBandwidth);
			MinYieldAtAppEvent(secondBound, cTime,
					i,
					nApps, remCommVol, releaseDate, progress, maxBandwidth,
					availableBandwidth,
					tmpAllocatedBandwidth, tmpYields);

			keepBestSolution(nApps, tmpYields, tmpAllocatedBandwidth, yields, allocatedBandwidth);
		}
	}
}

//void MinYieldAtEventInInterval(double firstBound, double secondBound, double cTime,
//		unsigned long nApps, double * remCommVol, double * releaseDate, double * progress, double * maxBandwidth,
//		bool * processedApplis,
//		double availableBandwidth,
//		double * allocatedBandwidth, double * yields, double * bestDelta){
//
//	// Data structure to hold temporary results
//	double tmpAllocatedBandwidth[nApps];
//	double tmpYields[nApps];
//	// Initialization of results to be produced
//	for(int i=0; i<nApps; i++){
//		allocatedBandwidth[i] = 0;
//		yields[i] = 0;
//	}
//
//	// We loop over all the (non already processed) applications as possible constraining application
//	// if their maximum admissible bandwidth would enable them to complete in the considered time window
//	for(int i=0; i<nApps; i++){
//		if((!processedApplis[i])&&(firstBound*maxBandwidth[i]>=remCommVol[i]-ERROR)){
//			MinYieldAtAppEventInInterval(firstBound, secondBound, cTime,
//					i,
//					nApps, remCommVol, releaseDate, progress, maxBandwidth,
//					availableBandwidth,
//					tmpAllocatedBandwidth, tmpYields, bestDelta);
//
//			keepBestSolution(nApps, tmpYields, tmpAllocatedBandwidth, yields, allocatedBandwidth);
//		}
//	}
//}




bool LexMinYield(double delta, double cTime,
		unsigned long nApps, double * remCommVol, double * releaseDate, double * progress, double * maxBandwidth,
		double * allocatedBandwidths, double * yields, double &totBand){

	if (ALLDEBUGNELEXMIN)
		cout << endl << "\t\t***Entering LexMinYield***" << endl;

	assert(delta >= 0);

	// Amount of bandwidth that remains to be distributed
	double remainingBandwidth = totBand;
	// Data structure to hold the QUALITY of the solution
	vector<double> achievedYields;
	// Data structure to record for which applications the bandwidth (and thus the yield) has not already been fixed
	bool freeApplis[nApps];
	for(unsigned long i=0; i<nApps; i++)
	{
		freeApplis[i] = true;
	}
	int nbFreeApps = nApps;
	bool constrainingApps[nApps];

	// While there are free applications and bandwidth to be distributed
	while((nbFreeApps>0) && (remainingBandwidth > ERROR))
	{
		int formerNBFreeApps = nbFreeApps;

		if (ALLDEBUGNELEXMIN)
			cout << "\t\t   remainingBandwidth= " << remainingBandwidth << " for " << nbFreeApps << " free applications" << endl;

		// If there is a simple free application the solution is obvious
		if (nbFreeApps==1){
			for(unsigned long i=0; i<nApps; i++)
			{
				if (freeApplis[i])
				{
					if (remainingBandwidth<maxBandwidth[i]){
						allocatedBandwidths[i] = remainingBandwidth;
					}
					else{
						allocatedBandwidths[i] = maxBandwidth[i];
					}
					freeApplis[i] = false;
					nbFreeApps--;
					// We record its yield
					achievedYields.push_back(fYield(i, delta, cTime, allocatedBandwidths[i], releaseDate, progress, maxBandwidth));
				}
			}
		}
		else{

			// Compute the best achievable minimum yield for the free applications
			bool bwIsExhausted;
			long double yield = MinYield(delta, cTime,
					nApps,
					freeApplis,
					remCommVol,  releaseDate,  progress,  maxBandwidth,
					remainingBandwidth, totBand,
					&bwIsExhausted, constrainingApps);

			if (ALLDEBUGNELEXMIN) cout << "\t\t\tTarget yield= " << yield << endl;

			// Compute the needed bandwidth
			double requiredBW = totalRequiredBandwidth(nApps,
					freeApplis,
					yield, delta, cTime,
					remCommVol, releaseDate, progress, maxBandwidth);

			//		// We should always use some of the bandwidth
			//		assert(requiredBW>0);

			// Sanity check
			if (ALLDEBUGNELEXMIN) cout << "\t\t\trequiredBW= " << requiredBW << "\tremainingBandwidth= " << remainingBandwidth << "\tBandwidth left: " << remainingBandwidth-requiredBW<< endl;
			//		cout <<"\t\t\tError margin: " << nApps*ERROR << endl;
			//		assert(requiredBW <= remainingBandwidth+nApps*ERROR);
			// If the whole remaining bandwidth is used we are done
			//		if (requiredBW >= remainingBandwidth - nbFreeApps*ERROR)
			if (bwIsExhausted||(requiredBW >= remainingBandwidth - nbFreeApps*ERROR)||(remainingBandwidth -requiredBW < remainingBandwidth/10000))
			{
				if (ALLDEBUGNELEXMIN){
					cout << "\t\t\tWhole bandwidth is used" << endl;
					if (bwIsExhausted & (!(requiredBW >= remainingBandwidth - nbFreeApps*ERROR))){
						cout << "\t\tPOTENTIAL PROBLEM: " << requiredBW << " ?>=? " << remainingBandwidth - nbFreeApps*ERROR << endl;
					}
				}
				long double scalingFactor = remainingBandwidth/requiredBW;
				for(unsigned long i=0; i<nApps; i++)
				{
					if (freeApplis[i])
					{
						// We set the bandwidth allocated to application i
						if (yield > fYield(i, delta, cTime, 0, releaseDate, progress, maxBandwidth))
						{
							allocatedBandwidths[i] = (maxBandwidth[i]/delta)*(yield*(cTime+delta-releaseDate[i])-progress[i])*scalingFactor;
							// To avoid rounding problems
							if (allocatedBandwidths[i] > maxBandwidth[i]){
								allocatedBandwidths[i] = maxBandwidth[i];
							}
							if (allocatedBandwidths[i] > remainingBandwidth){
								allocatedBandwidths[i] = remainingBandwidth;
							}
						}
						else
						{
							allocatedBandwidths[i] = 0;
						}
						if (ALLDEBUGNELEXMIN) cout << "\t\t\t\tBW allocated to application "<<i<<": " << allocatedBandwidths[i] << endl;
						remainingBandwidth -= allocatedBandwidths[i];
						freeApplis[i] = false;
						nbFreeApps--;
						// We record its yield
						achievedYields.push_back(fYield(i, delta, cTime, allocatedBandwidths[i], releaseDate, progress, maxBandwidth));
					}
				}
				//			remainingBandwidth = 0;
			}
			// Otherwise we identify the applications whose yield is equal to the optimal minimum yield
			else
			{
				// If requiredBW is equal to 0, this means that the remaining bandwidth is so small that allocating it
				// does not change the minimum yield. The processing of this case is a hack to cope with rounding problems
				if(requiredBW==0){
					// We allocate the remaining bandwidth equally among the free apps
					int nbAllocatedApps = 0;
					for(unsigned long i=0; i<nApps; i++)
					{
						if (freeApplis[i])
						{
							allocatedBandwidths[i] = remainingBandwidth/nbFreeApps;
							freeApplis[i] = false;
							nbAllocatedApps++;
							// We record its yield
							achievedYields.push_back(fYield(i, delta, cTime, allocatedBandwidths[i], releaseDate, progress, maxBandwidth));
						}
					}
					assert(nbAllocatedApps==nbFreeApps);
					remainingBandwidth = 0;
					nbFreeApps=0;
				}
				else{
					if (ALLDEBUGNELEXMIN) cout << "\t\t\tSome bandwidth will be left" << endl;
					for(unsigned long i=0; i<nApps; i++)
					{
						if (freeApplis[i])
						{
							// Compute the bandwidth needed by application i to achieve the target yield
							double bwForAppli = 0;
							if (yield > fYield(i, delta, cTime, 0, releaseDate, progress, maxBandwidth))
							{
								bwForAppli = (maxBandwidth[i]/delta)*(yield*(cTime+delta-releaseDate[i])-progress[i]);
								if(bwForAppli > maxBandwidth[i]){
									bwForAppli = maxBandwidth[i];
								}
								if(bwForAppli > remainingBandwidth){
									bwForAppli = remainingBandwidth;
								}
							}
							else
							{
								bwForAppli = 0;
							}
							// Sanity check
							assert(bwForAppli < maxBandwidth[i]+ERROR);

							if (ALLDEBUGNELEXMIN)
							{
								cout << "\t\t\t\tApplication " << i << " BW: " << bwForAppli << " (<= max: " << maxBandwidth[i] << ") ";
								cout << "\tYield: " << fYield(i, delta, cTime, bwForAppli, releaseDate, progress, maxBandwidth);
								cout << "\tRemaining comm volume: " << remCommVol[i] - delta*bwForAppli << endl;
							}

							//					cout << "**MYDEBUG**:\t" << bwForAppli << "?>=?" << maxBandwidth[i]-ERROR << " \tOR\t " << bwForAppli*delta << "?>=?" <<  remCommVol[i]-ERROR << endl;


							// Check whether application i is constraining
							//					if ((bwForAppli >= maxBandwidth[i]-ERROR) || (bwForAppli*delta >= remCommVol[i]-ERROR) || (bwForAppli == 0))
							//						cout << "\t\t\t\t" << bwForAppli << " ?>=? " << maxBandwidth[i]-ERROR << "   OR   " << bwForAppli << "*" << delta << " = " << bwForAppli*delta << " ?>=? " << remCommVol[i]-ERROR << endl;
							if ((bwForAppli >= maxBandwidth[i]-ERROR) || (bwForAppli*delta >= remCommVol[i]-ERROR) || constrainingApps[i])
							{

								if (ALLDEBUGNELEXMIN) cout << "\t\t\t\t\tAppli " << i << " is tight" << endl;
								// We set the bandwidth allocated to application i
								allocatedBandwidths[i] = bwForAppli;
								if (ALLDEBUGNELEXMIN) cout << "\t\t\t\t\tAllocated bandwidth: " << allocatedBandwidths[i] << endl;
								remainingBandwidth -= allocatedBandwidths[i];
								freeApplis[i] = false;
								nbFreeApps--;
								// We record its yield
								achievedYields.push_back(fYield(i, delta, cTime, allocatedBandwidths[i], releaseDate, progress, maxBandwidth));
							}
						}
					}
				}
			}
		}

		// If there remains some free applications (because the bandwidth was exhausted and they did not need any)
		// we complete the solution
		if(nbFreeApps>0){
			for(unsigned long i=0; i<nApps; i++){
				if (freeApplis[i]){
					allocatedBandwidths[i] = 0;
					achievedYields.push_back(fYield(i, delta, cTime, allocatedBandwidths[i], releaseDate, progress, maxBandwidth));
				}
			}
		}


		//		cout << formerNBFreeApps << " > " << nbFreeApps << " OR " << remainingBandwidth << " < " << nbFreeApps*ERROR << endl;
		if(CATCHERROR){
			if ((formerNBFreeApps<=nbFreeApps)&(remainingBandwidth>nbFreeApps*ERROR)){
				return false;
			}
		}
		else{
			assert((formerNBFreeApps>nbFreeApps) || (remainingBandwidth<nbFreeApps*ERROR));
		}

	}

	// We sort the achieved yields
	sort(achievedYields.begin(), achievedYields.end());
	// We store the lexicographically sorted vector of achieved yields
	for(unsigned long i=0; i<nApps; i++)
	{
		yields[i] = achievedYields[i];
	}
	return true;
}


void computeGreedyCom(double cTime,vector<unsigned long> appInCom, double remainingBW, std::vector<Application> &apps)
{
    const unsigned long m=appInCom.size();
    pair<double,long> yields[m];
    double bandwidths[m];

    unsigned long i,j;
    unsigned long ind;

    // Store all bandwidth values of apps communicating and sort those values by increasing order
    for(i=0; i<m; i++)
    {
        ind=appInCom[i];
        yields[i]=make_pair(apps[ind].getRemVol()/apps[ind].getMaxBandwidth(),i);
        apps[ind].setBW(-1);
        bandwidths[i]=apps[ind].getMaxBandwidth();
    }

    size_t len = sizeof(yields) / sizeof(yields[0]);
    sort(yields, yields + len);
    
    
    j=0;
    while (remainingBW>0 and j<m){
        i=get<1>(yields[j]);
        ind=appInCom[i];
        //cout << apps[ind].name << " AAA " << get<0>(yields[j])<< "\n";
        if (remainingBW>=bandwidths[i]){
            apps[ind].setBW(bandwidths[i]);
            remainingBW-=bandwidths[i];
        }
        else
        {
            apps[ind].setBW(remainingBW);
            remainingBW=0;
        }
        j+=1;
    }
    while(j<m){
        i=get<1>(yields[j]);
        ind=appInCom[i];
        apps[ind].setBW(0);
        j+=1;
    }
}

void computeFixedWindow(double cTime, vector<unsigned long> appInCom, std::vector<Application> &apps, double &totBand, double delta, bool &failed){
    
    bool succeedLex;
	if (DEBUGNELEXMIN) cout << "***computeFixedWindow***" << endl;

    const unsigned long m=appInCom.size();  // Number of applications communicating
    double releaseDate[m];             // List of release dates, whose index corresponds to the one of appInCom
    double progress[m];           // List of values representing the time required by the application to progress as much as it would do if it was running solo.
    double maxBandwidth[m];            // List of maxBandwidth of all applications
    double remCommVol[m];               // List of the volume of communications remaining

    // The total bandwidth available is a global variable called "totBand"
    // The end of the window is a global variable called "endWin"
    // All application data are stored in a global vector called "apps"


    // Get all input values
    for (unsigned long i=0; i<m; i++) // For more details, check/test the function printApplication
    {
    	unsigned long ind=appInCom[i];
        releaseDate[i]=apps[ind].release;
        maxBandwidth[i]=apps[ind].getMaxBandwidth();
        progress[i]=apps[ind].getTProgress();
        remCommVol[i]=apps[ind].getRemVol();
    }

    // Initialization to be on the safe side (the bandwidth allocated to non communicating applications will therefore be zero)
	
    // Compute the sum of the requested bandwidths
    double totalRequestedBandwidth = 0;
    for(unsigned long i = 0; i<m; i++){
    	totalRequestedBandwidth += maxBandwidth[i];
    }
    // If the total requested bandwidth does not exceed the total available bandwidth the solution is trivial
    if (totalRequestedBandwidth <= totBand + ERROR){
        //Allocate its maximum bandwidth to any communicating application
        for (unsigned long i=0; i<m; i++)
        {
        	unsigned long ind=appInCom[i];
            apps[ind].setBW(apps[ind].getMaxBandwidth());
        }
        return;
    }

    // Data structures to hold the current solution and the current best solution
    // Initialized to dummy default
    
    double tempBandwidths[m];
    double tempYields[m];

    // We check whether the end of the processing window can be a solution
	if (DEBUGNELEXMIN) cout << "\tChecking potential solution at end of processing window";
    succeedLex=LexMinYield(delta, cTime, m, remCommVol, releaseDate, progress, maxBandwidth, tempBandwidths, tempYields, totBand);
    
    if (not succeedLex){
        computeGreedyCom(cTime,appInCom, totBand, apps);
        failed =true;
        return;
    }
    
    double totBandwi=0;
    for (unsigned long i=0; i<m; i++)
    {
        unsigned long ind=appInCom[i];
        apps[ind].setBW(tempBandwidths[i]);
        totBandwi+=tempBandwidths[i];
    }
    totBandwi=totBand-totBandwi;
    unsigned long i=0;
    while(totBandwi>ERROR and i<m){
        unsigned long ind=appInCom[i];
        double xx=apps[ind].getMaxBandwidth()-tempBandwidths[i];
        if(xx<totBandwi){
            apps[ind].setBW(apps[ind].getMaxBandwidth());
            totBandwi-=xx;
        }
        else{
            apps[ind].setBW(tempBandwidths[i]+xx);
            totBandwi=0;
        }
        i+=1;
    }
}



void computeNextEvLexMin(double cTime, vector<unsigned long> appInCom, std::vector<Application> &apps, double &totBand, double &endWin, bool &failed)  // cTime is the current time
{
	// Variable to enable the debug of a peculiar case
	bool PECULIAR_CASE_DEBUG = false;
//	if ((cTime > 315) && (cTime<350)){
//		PECULIAR_CASE_DEBUG = true;
//	}
//	if (cTime>350){
//		abort();
//	}

	// appInCom contains the list of indexes of applications communicating

    bool succeedLex;
	if (DEBUGNELEXMIN || PRINTCOMPUTENEXTEVENT)
		cout << "\n***computeNextEvLexMin*** at time " << cTime << endl;

    const unsigned long m=appInCom.size();  // Number of applications communicating

    // We check that we are not in a degenerate case
    if (m==0){
    	return;
    }

    double releaseDate[m];             // List of release dates, whose index corresponds to the one of appInCom
    double progress[m];           // List of values representing the time required by the application to progress as much as it would do if it was running solo.
    double maxBandwidth[m];            // List of maxBandwidth of all applications
    double remCommVol[m];               // List of the volume of communications remaining

    // The total bandwidth available is a global variable called "totBand"
    // The end of the window is a global variable called "endWin"
    // All application data are stored in a global vector called "apps"




    // Get all input values
    for (unsigned long i=0; i<m; i++) // For more details, check/test the function printApplication
    {
    	unsigned long ind=appInCom[i];
        releaseDate[i]=apps[ind].release;
        maxBandwidth[i]=apps[ind].getMaxBandwidth();
        progress[i]=apps[ind].getTProgress();
        remCommVol[i]=apps[ind].getRemVol();
    }

    for (unsigned long i=0; i<m; i++){
    	if (ALLDEBUGNELEXMIN || CHECKOPTIMISATION) cout << "\tRemaining Comm Volume for Application " << i << ": " << remCommVol[i] << " y= " << progress[i]/(cTime-releaseDate[i]) << " rel= " << releaseDate[i] << endl;
    	assert(remCommVol[i]>0);
    	assert(maxBandwidth[i]>0);
    }

    // Initialization to be on the safe side (the bandwidth allocated to non communicating applications will therefore be zero)
	for (unsigned long i=0; i<apps.size(); i++)
    {
        apps[i].setBW(0);
    }

    // Data structures to hold the current solution and the current best solution
    // Initialized to dummy default
    double bestTempBandwidths[m];
    double bestTempYields[m];
    for (unsigned long i=0; i<m; i++)
    {
    	bestTempYields[i] = -1;
    	bestTempBandwidths[i] = 0;
    }

    double tempBandwidths[m];
    double tempYields[m];



    // Compute the sum of the requested bandwidths
    double totalRequestedBandwidth = 0;
    for(unsigned long i = 0; i<m; i++){
    	totalRequestedBandwidth += maxBandwidth[i];
    }
    // If the total requested bandwidth does not exceed the total available bandwidth the solution is trivial
    if (totalRequestedBandwidth <= totBand + ERROR){
        //Allocate its maximum bandwidth to any communicating application
        for (unsigned long i=0; i<m; i++)
        {
        	unsigned long ind=appInCom[i];
            apps[ind].setBW(apps[ind].getMaxBandwidth());
        }
    	if (DEBUGNELEXMIN) cout << "\tEach application is allocated its maximum bandwidth" << endl;
        return;
    }


	// We eliminate degenerate cases: if a communication lasts less 0.01 second we use a crude solution
	if(HACKSHORTCOMM){
		if (ALLDEBUGNELEXMIN) cout << "In HACKSHORTCOMM " << endl;
		// We compute the duration of the shortest communication
		double shortestComm = remCommVol[0]/maxBandwidth[0];
		unsigned long shortestApp = 0;
	    for(unsigned long i = 1; i<m; i++){
	    	double commLength = remCommVol[i]/maxBandwidth[i];
	    	if (shortestComm > commLength){
	    		shortestComm = commLength;
	    		shortestApp = i;
	    	}
	    }
	    // If it is too short we use a simplified bandwidth allocation
	    if(shortestComm < 0.001){
	    	if (ALLDEBUGNELEXMIN) cout << "\tIn HACKSHORTCOMM : there is a short communication." << endl;
	    	bestTempBandwidths[shortestApp] = maxBandwidth[shortestApp];
	    	double remainingBandwidth = totBand - bestTempBandwidths[shortestApp];
	    	// We compute the total bandwidth reuested by the other applications
	    	double totRemainingRequestedBW = totalRequestedBandwidth-bestTempBandwidths[shortestApp];
	    	assert(totalRequestedBandwidth > remainingBandwidth);
		    for(unsigned long i = 1; i<m; i++){
		    	if(i!=shortestApp){
		    		bestTempBandwidths[i] = remainingBandwidth*maxBandwidth[i]/totRemainingRequestedBW;
		    		if(bestTempBandwidths[i] > maxBandwidth[i]){
		    			bestTempBandwidths[i] = maxBandwidth[i];
		    		}
		    	}
		    }
	    	// Finally we store the solution
	        for (unsigned long i=0; i<m; i++)
	        {
	        	unsigned long ind=appInCom[i];
	        	apps[ind].setBW(bestTempBandwidths[i]);
	        }
	        return;
	    }
	    else{
	    	if (DEBUGNELEXMIN) cout << "\tIn HACKSHORTCOMM : NO short communication. Shortest communication of length : " << shortestComm << endl;
	    }
	}


    // We check whether the end of the processing window can be a solution
	if (DEBUGNELEXMIN) cout << "\tChecking potential solution at end of processing window";
    succeedLex=LexMinYield(endWin, cTime, m, remCommVol, releaseDate, progress, maxBandwidth, tempBandwidths, tempYields, totBand);
    if (not succeedLex){
        computeGreedyCom(cTime,appInCom, totBand, apps);
        failed=true;
        return;
    }
    
    double tempBandwidthUsage = 0;
    for (unsigned long i=0; i<m; i++)
    {
    	tempBandwidthUsage += tempBandwidths[i];
    }

//    if ((tempYields[0] > 0) && (tempBandwidthUsage > totBand-ERROR)){
    if (tempYields[0] > 0){
    	if (DEBUGNELEXMIN) cout << "\t\tEnd of processing window is a potential optimal solution with bandwidth usage : " << tempBandwidthUsage << endl;
        for (unsigned long i=0; i<m; i++)
        {
        	bestTempBandwidths[i] = tempBandwidths[i];
        	bestTempYields[i] = tempYields[i];
        }
    }
    else
    {
    	if (DEBUGNELEXMIN) cout << "\t\tEnd of processing window is NOT a potential optimal solution (yield:  " << tempYields[0] << ", bandwidth usage: " << tempBandwidthUsage << ")" << endl ;
    }


    // We build the set of interesting times
    std::vector<double> interestingTimes = {endWin};
    // We consider the first time each communication can end
    bool noPossibleCommEndingTimeInWindow = true;
    double overallEarliestCommEndingTime = endWin;
    for(unsigned long i = 0; i<m; i++){
    	double earliestCommEndingTime = remCommVol[i]/maxBandwidth[i];
    	if (earliestCommEndingTime < endWin)
    	{
    		interestingTimes.push_back(earliestCommEndingTime);
    		if(noPossibleCommEndingTimeInWindow || (earliestCommEndingTime < overallEarliestCommEndingTime)){
    			overallEarliestCommEndingTime = earliestCommEndingTime;
    		}
    		noPossibleCommEndingTimeInWindow = false;
    	}

    }

    // If no communication can end during the time window the next event can only be at the window end
    // Which is the solution stored in bestTempBandwidths
    if (noPossibleCommEndingTimeInWindow)
    {
        for (unsigned long i=0; i<m; i++)
            {
                unsigned long ind=appInCom[i];
                apps[ind].setBW(bestTempBandwidths[i]);
            }
        return;
    }

    // We consider the dates at which two communications can end simultaneously
    for(unsigned long i = 0; i<m-1; i++)
    {
    	double denominatorFirstComm = progress[i] + remCommVol[i]/maxBandwidth[i];

    	for(unsigned long j = i+1; j<m; j++)
	    {
    		double denominatorSecondComm = progress[j] + remCommVol[j]/maxBandwidth[j];
    		// If the denominator will not be equal to zero
    		if ((denominatorFirstComm-denominatorSecondComm >= ERROR) || (denominatorFirstComm-denominatorSecondComm <= -ERROR))
    		{
    			double candiDate = (denominatorSecondComm*(cTime-releaseDate[i])-denominatorFirstComm*(cTime-releaseDate[j]))
    					/(denominatorFirstComm-denominatorSecondComm);
    			if ((endWin >= candiDate) && (candiDate >= remCommVol[i]/maxBandwidth[i]) && (candiDate >= remCommVol[j]/maxBandwidth[j]))
    			{
    	    		interestingTimes.push_back(candiDate);
    			}
    		}
	    }
    }

    // We consider the dates at which a communication can end while another is receiving its maximum bandwidth and both have the same yield
    for(unsigned long i = 0; i<m; i++)
    {
    	for(unsigned long j = 0; j<m; j++)
    	{
    		if (i != j){
    			double polynomB = progress[i]-progress[j]+cTime-releaseDate[j]-remCommVol[j]/maxBandwidth[j];
    			double polynomC = progress[i]*(cTime-releaseDate[j])-(progress[j]+remCommVol[j]/maxBandwidth[j])*(cTime-releaseDate[i]);
    			double discriminant = polynomB*polynomB-4*polynomC;

    			// If the second degree polynomial admits real solutions
    			if (discriminant >= 0){
    				// We record the first root if its belong to the window
    				double firstRoot = (-polynomB-sqrt(discriminant))/2;
    				if((firstRoot >= 0) && (firstRoot <= endWin) && (firstRoot >= remCommVol[j]/maxBandwidth[j]))
    				{
        	    		interestingTimes.push_back(firstRoot);
    				}
    				// We record the second root if its belong to the window and its different from the first
    				double secondRoot = (-polynomB+sqrt(discriminant))/2;
    				if((discriminant > 0) && (secondRoot >= 0) && (secondRoot <= endWin) && (secondRoot >= remCommVol[j]/maxBandwidth[j]))
    				{
        	    		interestingTimes.push_back(secondRoot);
    				}
    			}
    		}
    	}
    }

    // We consider the dates at which two communications receive their maximum bandwidth and have the same yield
    for(unsigned long i = 0; i<m; i++)
    {
    	for(unsigned long j = 0; j<m; j++)
    	{
    		if (i != j){
    			double denominator = progress[i]-progress[j]+releaseDate[i]-releaseDate[j];
    			if ((denominator > ERROR) || (denominator < -ERROR)){
    				double candiDate = (progress[j]*(cTime-releaseDate[i])-progress[i]*(cTime-releaseDate[j]))/denominator;
    				if ((candiDate >= overallEarliestCommEndingTime) && (candiDate <= endWin)){
    					interestingTimes.push_back(candiDate);
    				}
    			}
    		}
    	}
    }

    // We consider the dates at which a communication can end while achieving the same yield than an application that does not receive any bandwidth
    for(unsigned long i = 0; i<m; i++)
    {
    	for(unsigned long j = 0; j<m; j++)
    	{
    		if (i != j){
    			double denominator = progress[i]-progress[j]-remCommVol[j]/maxBandwidth[j];
    			if ((denominator > ERROR) || (denominator < -ERROR)){
    				double candiDate = ((progress[j]+remCommVol[j]/maxBandwidth[j])*(cTime-releaseDate[i])-progress[i]*(cTime-releaseDate[j]))/denominator;
    				if ((candiDate >= remCommVol[j]/maxBandwidth[j]) && (endWin >= candiDate))
    				{
        	    		interestingTimes.push_back(candiDate);
    				}
    			}
    		}
    	}
    }

    // We consider the dates at which two applications that do not receive any bandwidth have the same yield
    for(unsigned long i = 0; i<m; i++)
    {
    	for(unsigned long j = 0; j<m; j++)
    	{
    		if (i != j){
    			double denominator = progress[j]-progress[i];
    			if ((denominator > ERROR) || (denominator < -ERROR)){
    				double candiDate = (progress[i]*(cTime-releaseDate[j])-progress[j]*(cTime-releaseDate[i]))/denominator;
    				if ((candiDate >= overallEarliestCommEndingTime) && (endWin >= candiDate))
    				{
        	    		interestingTimes.push_back(candiDate);
    				}
    			}
    		}
    	}
    }

    // We sort the set of interesting times
    sort(interestingTimes.begin(), interestingTimes.end());

    // Sanity check debug
    if (DEBUGNELEXMIN)
    {
    	cout << "End of window: " << endWin << "\t Intervals: ";
    	for(unsigned long k=0; k<interestingTimes.size(); k++)
    	{
    		cout << interestingTimes[k] << "\t";
    	}
    	cout << endl;
    }

    // Data structure to record which applications have already been processed
    bool processedApplis[m];
    for(unsigned long i=0; i<m; i++){
    	processedApplis[i] = false;
    }
    int nbUnprocessedApplis = m;

    bool OverallMaximumNotFound = true;

    // We loop other the intervals defined by interestingTimes
//    for(unsigned long k=0; (k<interestingTimes.size()-1)&&(nbUnprocessedApplis>0); k++)
        for(unsigned long k=0; (k<interestingTimes.size()-1)&&(nbUnprocessedApplis>0)&&(OverallMaximumNotFound); k++)
    {

    	if (PECULIAR_CASE_DEBUG){
    		cout << "*** for time " << interestingTimes[k] << " BW = " << bestTempBandwidths[57] << endl;
    	}

	    if((DEBUGNELEXMIN)||  PECULIAR_CASE_DEBUG ){
	    	cout << "\tBeginning of step " << k << endl;
	        for (unsigned long i=0; i<m; i++){
	        	cout << i << ": " << bestTempBandwidths[i] << " (" << maxBandwidth[i] << ") ";
	        }
	        cout << endl;
	    }


    	// Sanity check to avoid side effects: we ignore intervals of negligible lengths
    	if(interestingTimes[k+1]-interestingTimes[k]<ERROR)
    	{
    		continue;
    	}

    	// The middle of the interval
    	double middle = (interestingTimes[k]+interestingTimes[k+1])/2;

    	if (DEBUGNELEXMIN) cout << endl << "   Search interval: [" << interestingTimes[k] << ", " << interestingTimes[k+1] << "]" << endl;
    	if (DEBUGNELEXMIN) cout << "\t\tMiddle of search interval: " << middle << endl;

    	// Identification of an application that defines the minimum achievable yield on this interval
    	int constrainingApp = -1;
    	double targetMinYield = -1;
    	bool initializedMinYield = false;
        for(unsigned long i = 0; i<m; i++)
        {
        	// Any constraining application should be able to complete its communication by time middle (or by time interestingTimes[k])
        	if (middle >= remCommVol[i]/maxBandwidth[i]){
        		double currentYield = (progress[i]+remCommVol[i]/maxBandwidth[i])/(cTime+middle-releaseDate[i]);
        		if ((!initializedMinYield) || (currentYield < targetMinYield)){
        			constrainingApp = i;
        			targetMinYield = currentYield;
        			initializedMinYield = true;
        		}
        	}
        }
        assert(initializedMinYield);

        if (DEBUGNELEXMIN) cout << "\t\tTarget min yield: " << targetMinYield << endl;
        if (DEBUGNELEXMIN) cout << "\t\tMost constraining application: " << constrainingApp << endl;


        // If the constraining applications has not already been processed we try it
        if(!processedApplis[constrainingApp]){
        	double proposedDelta;
        	MinYieldAtAppEventInInterval(interestingTimes[k], interestingTimes[k+1], cTime,
        			constrainingApp,
					m, remCommVol, releaseDate, progress, maxBandwidth,
					totBand,
					tempBandwidths, tempYields, &proposedDelta);
        	// If the solution found is achieved with the constrainingApp achieving the minimum yield, it is marked as processed
        	if (tempYields[0] >= fYield(constrainingApp, proposedDelta, cTime, remCommVol[constrainingApp]/proposedDelta, releaseDate, progress, maxBandwidth) - ERROR){
        		processedApplis[constrainingApp] = true;
        		nbUnprocessedApplis--;
        		OverallMaximumNotFound = false;
        		if (CHECKOPTIMISATION) cout << "\t*** Found last possible optimal(1): application " << constrainingApp << " at time " << proposedDelta << endl;
        	}
        	keepBestSolution(m, tempYields, tempBandwidths, bestTempYields, bestTempBandwidths);
        }

	    if((DEBUGNELEXMIN)||  PECULIAR_CASE_DEBUG ){
	    	cout << "\tAfter most constraining application " << endl;
	        for (unsigned long i=0; i<m; i++){
	        	cout << i << ": " << bestTempBandwidths[i] << " (" << maxBandwidth[i] << ") ";
	        }
	        cout << endl;
	    }

        // We then try any other application that can define an event in (the middle of) the target interval
        for(unsigned long i = 0; i<m; i++){
        	if ((!processedApplis[i]) && (i!=constrainingApp) && (middle >= remCommVol[i]/maxBandwidth[i])){
        		double proposedDelta;
        		MinYieldAtAppEventInInterval(interestingTimes[k], interestingTimes[k+1], cTime,
        				i,
						m, remCommVol, releaseDate, progress, maxBandwidth,
						totBand,
						tempBandwidths, tempYields, &proposedDelta);
        		// If the solution found is achieved with application i achieving the minimum yield, it is marked as processed
        		if (tempYields[0] >= fYield(i, proposedDelta, cTime, remCommVol[i]/proposedDelta, releaseDate, progress, maxBandwidth) - ERROR){
        			processedApplis[i] = true;
            		nbUnprocessedApplis--;
            		OverallMaximumNotFound=false;
            		if (CHECKOPTIMISATION) cout << "\t*** Found last possible optimal(2): application " << i << " at time " << proposedDelta << endl;
        		}

        	    if((DEBUGNELEXMIN)||  PECULIAR_CASE_DEBUG ){
        	    	cout << "\tCandidate solution: for application " << i << endl;
        	        for (unsigned long i=0; i<m; i++){
        	        	cout << i << ": " << tempBandwidths[i] << " (" << maxBandwidth[i] << ") ";
        	        }
        	        cout << endl;


        	    	cout << "\tBefore keepBestSolution: for application " << i << endl;
        	        for (unsigned long i=0; i<m; i++){
        	        	cout << i << ": " << bestTempBandwidths[i] << " (" << maxBandwidth[i] << ") ";
        	        }
        	        cout << endl;
        	    }

        		keepBestSolution(m, tempYields, tempBandwidths, bestTempYields, bestTempBandwidths);

        	    if((DEBUGNELEXMIN)||  PECULIAR_CASE_DEBUG ){
        	    	cout << "\tAfter keepBestSolution: for application " << i << endl;
        	        for (unsigned long i=0; i<m; i++){
        	        	cout << i << ": " << bestTempBandwidths[i] << " (" << maxBandwidth[i] << ") ";
        	        }
        	        cout << endl;
        	    }

        	}
        }

//        // If the constraining applications has already been processed we skip it
//        if(processedApplis[constrainingApp]){
//        	continue;
//        }
//
//        // Identification of an application that defines the maximum achievable yield on this interval
//        double upperBoundYield = -1;
//        for(unsigned long i = 0; i<m; i++)
//        {
//        	double currentYield = fYieldMax(i, middle, cTime, remCommVol, releaseDate, progress, maxBandwidth);
//        	if ((i==0) || (currentYield < upperBoundYield)){
//        		upperBoundYield = currentYield;
//        	}
//        }
//
//        if (DEBUGNELEXMIN) cout << "\t\tUpper bound yield: " << upperBoundYield << endl;
//
//        // If the upper bound is smaller than the target yield there is a solution where some applications have a yield smaller than the constraining application
//        if (upperBoundYield < targetMinYield-ERROR)
//        {
////        	cout << "CALL 1 to MinYieldAtEvent" << endl;
//
////        	MinYieldAtEvent(interestingTimes[k], interestingTimes[k+1], cTime,
////        			m, remCommVol, releaseDate, progress, maxBandwidth,
////					processedApplis,
////					totBand,
////					tempBandwidths, tempYields);
//
//        	MinYieldAtEventInInterval(interestingTimes[k], interestingTimes[k+1], cTime,
//        			m, remCommVol, releaseDate, progress, maxBandwidth,
//					processedApplis,
//					totBand,
//					tempBandwidths, tempYields);
//
//            keepBestSolution(m, tempYields, tempBandwidths, bestTempYields, bestTempBandwidths);
//        	continue;
//        }
//
//        // We check whether interestingTimes[k] is a possible solution
//        double attemptedYield = fYieldAtCompletion(constrainingApp, interestingTimes[k], cTime, remCommVol, releaseDate, progress, maxBandwidth);
//        // We compute the bandwidth required
//        double totalBW = totalRequiredBandwidth(m, attemptedYield, interestingTimes[k], cTime,
//        		remCommVol, releaseDate, progress, maxBandwidth);
//
//        if (DEBUGNELEXMIN) cout << "\t\tBW at interestingTimes[" << k << "]: " << totalBW << endl;
//
//        // Check whether the solution is feasible
//        if (totalBW <= totBand + ERROR)
//        {
//        	// If the solution is feasible we compute it and store it if its the best solution found so far
//        	processedApplis[constrainingApp] = true;
//        	nbUnprocessedApplis--;
//        	succeedLex=LexMinYield(interestingTimes[k], cTime, m, remCommVol, releaseDate, progress, maxBandwidth, tempBandwidths, tempYields, totBand);
//
//            if (not succeedLex){
//                computeGreedyCom(cTime,appInCom, totBand, apps);
//                failed =true;
//                return;
//            }
//            keepBestSolution(m, tempYields, tempBandwidths, bestTempYields, bestTempBandwidths);
//        }
//        // interestingTimes[k] is not a feasible solution
//        // Look in the interval (interestingTimes[k], interestingTimes[k+1])
//        else
//        {
//        	// Set of applications to which bandwidth must be allocated
//        	bool needBandwidth[m];
//        	for(unsigned long i=0; i<m; i++)
//        	{
//        		if (fYield(i, middle, cTime, 0, releaseDate, progress, maxBandwidth) <= targetMinYield)
//        		{
//            		needBandwidth[i] = true;
//        		}
//        		else
//        		{
//            		needBandwidth[i] = false;
//        		}
//        	}
//
//        	if (DEBUGNELEXMIN)
//        	{
//        		cout << "\t\tApplications requiring bandwidth: ";
//        		for(unsigned long i=0; i<m; i++)
//        		{
//        			cout << needBandwidth[i] << "\t";
//        		}
//        		cout << endl;
//        	}
//
//        	// Parameters to solve the degree two polynomial
//        	double lambda = - totBand*(cTime-releaseDate[constrainingApp]);
//        	for(unsigned long i = 0; i<m; i++)
//        	{
//        		if(needBandwidth[i])
//        		{
//        		lambda += (progress[constrainingApp]+remCommVol[constrainingApp]/maxBandwidth[constrainingApp]-progress[i])*maxBandwidth[i];
//        		}
//        	}
//        	double mu = 0;
//        	for(unsigned long i = 0; i<m; i++)
//        	{
//        		if(needBandwidth[i])
//        		{
//        			mu += maxBandwidth[i]*((progress[constrainingApp]+remCommVol[constrainingApp]/maxBandwidth[constrainingApp])*(cTime-releaseDate[i])
//        					- progress[i]*(cTime-releaseDate[constrainingApp]));
//        		}
//        	}
//        	// If the second degree equation has real solutions we proceed
//        	double discriminant = lambda*lambda+4*totBand*mu;
//
//
//        	if (DEBUGNELEXMIN)
//        	{
//        		double yieldAtK = fYieldAtCompletion(constrainingApp, interestingTimes[k], cTime, remCommVol, releaseDate, progress, maxBandwidth);
//        		// We compute the bandwidth required
//        		double BWAtK = totalRequiredBandwidth(m, yieldAtK, interestingTimes[k], cTime,
//        				remCommVol, releaseDate, progress, maxBandwidth);
//        		double yieldAtKPlus = fYieldAtCompletion(constrainingApp, interestingTimes[k+1], cTime, remCommVol, releaseDate, progress, maxBandwidth);
//        		// We compute the bandwidth required
//        		double BWAtKPlus = totalRequiredBandwidth(m, yieldAtKPlus, interestingTimes[k+1], cTime,
//        				remCommVol, releaseDate, progress, maxBandwidth);
//
//        		cout << "\t\tSanity check, bandwidths at the interval bounds: " << BWAtK << " (yield: " << yieldAtK << ") and " << BWAtKPlus << " (yield: " << yieldAtKPlus << ")" << endl;
//        	}
//
//        	if (discriminant >= 0){
//        		// We build and sort the roots
//        		double rootOne = (lambda+sqrt(discriminant))/(2*totBand);
//        		double rootTwo = (lambda-sqrt(discriminant))/(2*totBand);
//        		double smallestRoot, largestRoot;
//        		if (rootOne < rootTwo)
//        		{
//        			smallestRoot = rootOne;
//        			largestRoot = rootTwo;
//        		}
//        		else
//        		{
//        			smallestRoot = rootTwo;
//        			largestRoot = rootOne;
//        		}
//
//        		if (DEBUGNELEXMIN) cout << "\t\tThe roots are " << smallestRoot << " and " << largestRoot << endl;
//
//        		// Check whether any root can define a (better) solution
//        		if ((interestingTimes[k] <= smallestRoot) && (smallestRoot<= interestingTimes[k+1]))
//        		{
//                    // We  look for a solution in the interval [interestingTimes[k], smallestRoot]
////                	MinYieldAtEvent(interestingTimes[k], smallestRoot, cTime,
////                			m, remCommVol, releaseDate, progress, maxBandwidth,
////							processedApplis,
////        					totBand,
////        					tempBandwidths, tempYields);
//
//                	MinYieldAtEventInInterval(interestingTimes[k], smallestRoot, cTime,
//                			m, remCommVol, releaseDate, progress, maxBandwidth,
//        					processedApplis,
//        					totBand,
//        					tempBandwidths, tempYields);
//                	keepBestSolution(m, tempYields, tempBandwidths, bestTempYields, bestTempBandwidths);
//
//
//                	processedApplis[constrainingApp] = true;
//                	nbUnprocessedApplis--;
//                	succeedLex=LexMinYield(smallestRoot, cTime, m, remCommVol, releaseDate, progress, maxBandwidth, tempBandwidths, tempYields, totBand);
//
//                    if (not succeedLex){
//                        computeGreedyCom(cTime,appInCom, totBand, apps);
//                        failed =true;
//                        return;
//                    }
//                    keepBestSolution(m, tempYields, tempBandwidths, bestTempYields, bestTempBandwidths);
//
////                	cout << "CALL 2 to MinYieldAtEvent" << endl;
//
//
//
//        		}
//        		else
//        		{
//            		if ((interestingTimes[k] <= largestRoot) && (largestRoot <= interestingTimes[k+1]))
//            		{
//                        // We look for a solution in the interval [interestingTimes[k], largestRoot]
////                    	MinYieldAtEvent(interestingTimes[k], largestRoot, cTime,
////                    			m, remCommVol, releaseDate, progress, maxBandwidth,
////								processedApplis,
////            					totBand,
////            					tempBandwidths, tempYields);
//
//                    	MinYieldAtEventInInterval(interestingTimes[k], largestRoot, cTime,
//                    			m, remCommVol, releaseDate, progress, maxBandwidth,
//            					processedApplis,
//            					totBand,
//            					tempBandwidths, tempYields);
//
//                    	keepBestSolution(m, tempYields, tempBandwidths, bestTempYields, bestTempBandwidths);
//
////            			cout << "Before" << endl;
//
//                    	processedApplis[constrainingApp] = true;
//                    	nbUnprocessedApplis--;
//                    	succeedLex=LexMinYield(largestRoot, cTime, m, remCommVol, releaseDate, progress, maxBandwidth, tempBandwidths, tempYields, totBand);
//
////            			cout << "After : " << succeedLex << endl;
//
//                        if (not succeedLex){
//                            computeGreedyCom(cTime,appInCom, totBand, apps);
//                            failed =true;
//                            return;
//                        }
//                        keepBestSolution(m, tempYields, tempBandwidths, bestTempYields, bestTempBandwidths);
//
////                    	cout << "CALL 3 to MinYieldAtEvent" << endl;
//
//
//
////                        cout << "Evrything was OK for CALL 3" << endl;
//            		}
//        		}
//        	}
//        	else{
////            	cout << "CALL 4 to MinYieldAtEvent" << endl;
//
//                // We look for a solution in the interval [interestingTimes[k], interestingTimes[k+1]]
////            	MinYieldAtEvent(interestingTimes[k], interestingTimes[k+1], cTime,
////            			m, remCommVol, releaseDate, progress, maxBandwidth,
////						processedApplis,
////    					totBand,
////    					tempBandwidths, tempYields);
//            	MinYieldAtEventInInterval(interestingTimes[k], interestingTimes[k+1], cTime,
//            			m, remCommVol, releaseDate, progress, maxBandwidth,
//    					processedApplis,
//    					totBand,
//    					tempBandwidths, tempYields);
//
//                keepBestSolution(m, tempYields, tempBandwidths, bestTempYields, bestTempBandwidths);
//        	}
//        }
    }

    // Finally we store the solution
    for (unsigned long i=0; i<m; i++)
    {
    	unsigned long ind=appInCom[i];
    	apps[ind].setBW(bestTempBandwidths[i]);
    }

    if((DEBUGNELEXMIN)||  PECULIAR_CASE_DEBUG || CHECKOPTIMISATION){
    	long double deltaOpt = endWin;
        for (unsigned long i=0; i<m; i++){
        	if((bestTempBandwidths[i]>ERROR)&&(remCommVol[i]/bestTempBandwidths[i]<deltaOpt)){
        		deltaOpt = remCommVol[i]/bestTempBandwidths[i];
        	}
        }

        cout << "bestTempYields[0]: " << bestTempYields[0] << " for time " << cTime+deltaOpt << endl;

    	cout << "\tThe solution is: ";
        for (unsigned long i=0; i<m; i++){
        	cout << i << ": " << bestTempBandwidths[i] << " (" << maxBandwidth[i] << ") ";
        	cout << fYield(i, deltaOpt, cTime, bestTempBandwidths[i], releaseDate, progress, maxBandwidth) << "  ";
        }
        cout << endl;
    }

}

void computeEqual(vector<unsigned long> appInCom, std::vector<Application> &apps, double &totBand)
{
    const unsigned long m=appInCom.size();
    double totalBand=0;
    unsigned long i;
    unsigned long ind;
    double ratio=0;

    // Store all bandwidth values of apps communicating and sort those values by increasing order
    for(i=0; i<m; i++)
    {
        ind=appInCom[i];
        totalBand+=apps[ind].getMaxBandwidth();
    }
    ratio=totBand/totalBand;
    if (ratio>1){
        ratio=1;
    }
    
    //Allocate this bandwidth to all Applications
    for (i=0; i<m; i++)
    {
        ind=appInCom[i];
        apps[ind].setBW(apps[ind].getMaxBandwidth()*ratio);
    }
}

    
void computeFIFOCom(vector<unsigned long> appInCom, std::vector<Application> &apps, double &totBand)
{
    const unsigned long m=appInCom.size();
    pair<double,long> releases[m];
    double bandwidths[m];

    unsigned long i,j;
    unsigned long ind;
    double remainingBW=totBand;

    // Store all bandwidth values of apps communicating and sort those values by increasing order
    for(i=0; i<m; i++)
    {
        ind=appInCom[i];
        releases[i]=make_pair(apps[ind].getReleasePhase(),i);
        apps[ind].setBW(-1);
        bandwidths[i]=apps[ind].getMaxBandwidth();
    }

    size_t len = sizeof(releases) / sizeof(releases[0]);
    sort(releases, releases + len);
    
    
    j=0;
    while (remainingBW>0 and j<m){
        i=get<1>(releases[j]);
        ind=appInCom[i];
        if (remainingBW>=bandwidths[i]){
            apps[ind].setBW(bandwidths[i]);
            remainingBW-=bandwidths[i];
        }
        else
        {
            apps[ind].setBW(remainingBW);
            remainingBW=0;
        }
        j+=1;
    }
    while(j<m){
        i=get<1>(releases[j]);
        ind=appInCom[i];
        apps[ind].setBW(0);
        j+=1;
    }
}

void computeGreedyYield(double cTime,vector<unsigned long> appInCom, double remainingBW, std::vector<Application> &apps)
{
    const unsigned long m=appInCom.size();
    pair<double,long> yields[m];
    double bandwidths[m];

    unsigned long i,j;
    unsigned long ind;

    // Store all bandwidth values of apps communicating and sort those values by increasing order
    for(i=0; i<m; i++)
    {
        ind=appInCom[i];
        yields[i]=make_pair(apps[ind].getYield(cTime),i);
        apps[ind].setBW(-1);
        bandwidths[i]=apps[ind].getMaxBandwidth();
    }

    size_t len = sizeof(yields) / sizeof(yields[0]);
    sort(yields, yields + len);
    
    
    j=0;
    while (remainingBW>0 and j<m){
        i=get<1>(yields[j]);
        ind=appInCom[i];
        if (remainingBW>=bandwidths[i]){
            apps[ind].setBW(bandwidths[i]);
            remainingBW-=bandwidths[i];
        }
        else
        {
            apps[ind].setBW(remainingBW);
            remainingBW=0;
        }
        j+=1;
    }
    while(j<m){
        i=get<1>(yields[j]);
        ind=appInCom[i];
        apps[ind].setBW(0);
        j+=1;
    }
}



void computeLookAheadGreedyYield(double cTime,vector<unsigned long> appInCom, std::vector<Application> &apps, double &totBand, double &endWin){
    
    vector<unsigned long> appInCom2;
    const unsigned long m=appInCom.size();
    double minYields[m];
    double PBallocated;
    double yield;
    double minYield;
    double maxMinYield;
    double nextEvent,endPhase,newProgress;
    
    unsigned long i,j,ind,choice;
    unsigned long n=apps.size();
    
    // Store all bandwidth values of apps communicating and sort those values by increasing order
    if(m==0){
        return;
    }
    if(m==1){
        ind=appInCom[0];
        if(totBand>apps[ind].getMaxBandwidth()){
            apps[ind].setBW(apps[ind].getMaxBandwidth());
        }
        else{
            apps[ind].setBW(totBand);
        }
    }
    else{
        for(i=0; i<m; i++)
        {
            ind=appInCom[i];
            if(totBand>apps[ind].getMaxBandwidth()){
                PBallocated=apps[ind].getMaxBandwidth();
            }
            else{
                PBallocated=totBand;
            }
            
            apps[ind].setBW(PBallocated);
            
            appInCom2={};
            for (j=0;j<m;j++){
                if (j!=i){
                    appInCom2.push_back(appInCom[j]);
                }
            }
            
            computeGreedyYield(cTime,appInCom2,totBand-PBallocated, apps);
            
            nextEvent=endWin;
            
            for (j=0; j<n; j++){
                if (apps[j].getType())
                {
                    endPhase=apps[j].getRemVol()/apps[j].getBW();
                    if(endPhase+cTime<nextEvent)
                    {
                        nextEvent=cTime+endPhase;
                    }
                }
            }
            
            minYield=1;
            for(j=0;j<n;j++){
                if (apps[j].getType()){
                    newProgress=(nextEvent-cTime)*apps[j].getBW()/apps[j].getMaxBandwidth();
                    yield=(apps[j].getTProgress()+newProgress)/(nextEvent-apps[j].release);
                    if (yield<minYield){
                        minYield=yield;
                    }
                }
            }
            
            minYields[i]=minYield;
        }
        
        maxMinYield=0;
        choice=0;
        for(i=0;i<m;i++){
            if(minYields[i]>maxMinYield){
                choice=i;
                maxMinYield=minYields[i];
            }
        }
        
        
        ind=appInCom[choice];
        if(totBand>apps[ind].getMaxBandwidth()){
            PBallocated=apps[ind].getMaxBandwidth();
        }
        else{
            PBallocated=totBand;
        }
        
        apps[ind].setBW(PBallocated);
        
        appInCom2={};
        for (j=0;j<m;j++){
            if (j!=choice){
                appInCom2.push_back(appInCom[j]);
            }
        }
        
        computeGreedyYield(cTime,appInCom2,totBand-PBallocated, apps);
    }
}

void computeSet10(vector<unsigned long> appInCom, std::vector<Application> &apps, double &totBand)  // cTime is the current time
{
    // appInCom contains the list of indexes of applications communicating
    //cout<<"new\n";

    const unsigned long m=appInCom.size();  // Number of applications communicating
    double maxBandwidth[m];            // List of maxBandwidth of all applications
    double remCommVol[m];               // List of the volume of communications remaining
    double tempBand[m];
    long priorities[m];
    double remainingBW;
    double remBand;
    long Phase;
    long count;
    unsigned long doneprev=0;
    double Av=0;
    double bandj;
    bool conti=false;
    vector<pair<long, double> > BandSet={};
    bool okSet[m];
    
    bool found = false;
    unsigned long done=0;
    pair<double,long> releases[m];
    unsigned long ind;
    // The total bandwidth available is a global variable called "totBand"
    // The end of the window is a global variable called "endWin"
    // All application data are stored in a global vector called "apps"
    if(m==0){
        return;
    }
    unsigned long i;
    // Get all input values
    for (i=0; i<m; i++) // For more details, check/test the function printApplication
    {
    	ind=appInCom[i];
        maxBandwidth[i]=apps[ind].getMaxBandwidth();
        remCommVol[i]=apps[ind].getRemVol();
        Phase=apps[ind].getPhaseId();
        

        releases[i]=make_pair(apps[ind].getReleasePhase(),i);

        count=apps[ind].getCountIter();
        if(not apps[ind].knownWiter and apps[ind].getNew()){
            //cout << apps[ind].name << " lol " << Av << " " << apps[ind].getNew() << "\n";
            apps[ind].setNew(false);
            if (count<100000000){
                if (count==0){
                    //Av=apps[ind].Volumes[Phase]/apps[ind].getMaxBandwidth();
                    Av=0.01;
                }
                else{
                    Av=(apps[ind].getWiter()*(count-1)+apps[ind].Volumes[Phase-1]+apps[ind].Volumes[Phase]/apps[ind].getMaxBandwidth())/count;
                }
            }
            else{
                Av=0.1*(apps[ind].Volumes[Phase-1]+apps[ind].Volumes[Phase])+0.9*apps[ind].getWiter();
            }
            apps[ind].setCountIter(count+1);
            apps[ind].setWiter(Av);
        }
        else{
            Av=apps[ind].getWiter();
        }
        //cout << apps[ind].name << " lol " << Av << " " << apps[ind].getNew() << "\n";
        priorities[i]=round(log10(Av));
       // cout << apps[ind].name << " " << Av << "\n";
        found=false;
        for (unsigned long ii = 0; ii < BandSet.size(); ii++) {
            if (BandSet[ii].first == priorities[i]) {
                BandSet[ii].second +=maxBandwidth[i];
                found = true;
                break;
            }
        }

        if (!found) {
            BandSet.push_back(make_pair(priorities[i],maxBandwidth[i]));
        }
    }
    
    
    size_t len = sizeof(releases) / sizeof(releases[0]);
    sort(releases, releases + len);
    
    double totP=0;
    for (unsigned long j=0;j<BandSet.size();j++){
        okSet[j]=false;
        totP+=pow(10,-BandSet[j].first);
    }
    remBand=totBand;
    
    for (unsigned long i=0; i<m; i++)
    {
        tempBand[i]=0;
    }
            
    while(done<BandSet.size()){
        doneprev=done;
        /*cout << "new loop !!! " << done << "\n";
        for (unsigned long i=0; i<m; i++)
        {
            cout <<tempBand[i] << "\n";
        }*/
        conti=false;
        for (unsigned long j=0;j<BandSet.size();j++){ 
            if(conti){
                continue;
            }
            if(not okSet[j]){
                bandj=pow(10,-BandSet[j].first)/totP*remBand;
                if(bandj>BandSet[j].second){
                    for (i=0;i<m;i++){
                        if (priorities[i]==BandSet[j].first){
                            ind=appInCom[i];
                            tempBand[i]=maxBandwidth[i];
                            remBand-=maxBandwidth[i];
                        }
                    }
                    totP-=pow(10,-BandSet[j].first);
                    okSet[j]=true;
                    done+=1;
                    conti=true;
                    continue;
                }
                else{
                    
                    unsigned long jj=0;
                    remainingBW=bandj;
                    while (remainingBW>0 and jj<m){
                            i=get<1>(releases[jj]);
                            if(priorities[i]==BandSet[j].first){
                                if (remainingBW>=maxBandwidth[i]){
                                    tempBand[i]=maxBandwidth[i];
                                    remainingBW-=maxBandwidth[i];
                                }
                                else
                                {
                                    tempBand[i]=remainingBW;
                                    remainingBW=0;
                                }
                            }
                            jj+=1;
                    }
                    while(jj<m){
                        i=get<1>(releases[jj]);
                        if(priorities[i]==BandSet[j].first){
                            tempBand[i]=0;
                        }
                        jj+=1;
                    }
                }
            }
        }
        if (doneprev==done){
            break;
        }
        else{
            doneprev=done;
        }
    }
    
	for (unsigned long i=0; i<m; i++)
    {
        apps[appInCom[i]].setBW(tempBand[i]);
    }

}

//************** OTHER CORE FUNCTIONS****************

double compute(double cTime, string heuristic, vector<unsigned long> appInCom, string outName, std::vector<Application> &apps, double &totBand, double &endWin, double delta, bool &failed)    // Returns the next moment it needs to take a decision
{
    // The return is only useful for the algorithm that takes
    // new decisions regularly 
    if (heuristic=="fairShare")
    {
        computeEqual(appInCom, apps, totBand);
        return endWin;
    }
    else if(heuristic=="nextEvLexMin")
    {
        computeNextEvLexMin(cTime, appInCom, apps, totBand, endWin, failed);
        return endWin;
    }
    else if(heuristic=="FCFS"){
        computeFIFOCom(appInCom, apps, totBand);
        return endWin;
    }
    else if(heuristic=="fixedWindow")
    {
        computeGreedyYield(cTime, appInCom, totBand, apps);
        if ((cTime-((floor(cTime/delta)+1)*delta))/cTime<ERROR){
            return (cTime+delta);
        }
        else{
            return ((floor(cTime/delta)+1)*delta);
        }
    }
    else if(heuristic=="greedyYield"){
        computeGreedyYield(cTime,appInCom, totBand, apps);
        return endWin;
    }
    else if(heuristic=="greedyCom"){
        computeGreedyCom(cTime,appInCom, totBand, apps);
        return endWin;
    }
    else if(heuristic=="lookAheadGreedyYield"){
        computeLookAheadGreedyYield(cTime,appInCom, apps, totBand, endWin);
        return endWin;
    }
    else if(heuristic=="Set10" or heuristic=="Set10Learn"){
        computeSet10(appInCom, apps, totBand);
        return endWin;
    }
    else
    {
        std::stringstream sb;
        sb << "FATAL ERROR: Heuristic '" << heuristic << "' Not Recognized";
        throw std::runtime_error(sb.str());
    }
}
//************** SIMULATOR ****************

template <int debuga>
double simulate(string heuristic2, double totBand, double &endWin, std::vector<Application> &apps, string outName, bool changeHeuristic,double delta, bool &failed, string outdetails)
{
    double cTime = 0; //Current time
    vector<unsigned long> appInCom; //List of applications communicating
    const unsigned long n=apps.size(); //// Number of applications  now we may not add or terminate applications, this can change easily though
    // I used n because I let the m for your function, i.e. m=number of applications communicating
    
    unsigned long i; // Utilitary
    //For step B
    double nextDecision;
    double nextEvent;
    double endPhase;
    //For step C
    double newDVolume;
    double newRVolume;
    double newRVolume2;
    double newProgress;
    // Compute objective function
    double minyield=1;
    double curyield=0;
    double tempForTest=0;
    double totalBW=0;
    double totalPB=0;
    double initialProgress=0;
    double finalProgress=0;
    double timeChange=0;
    double start=0;
    double count=0;
    
    double counttest = 0;  //for DEBUG
    
    string heuristic=heuristic2;
    
    bool changedHeuristic=true;
    
    
    if(changeHeuristic){
        changedHeuristic=false;
        timeChange=endWin/2;
        heuristic="fairShare";
    }
    else{
        for (i=0;i<n;i++){
            initialProgress+=apps[i].getTProgress()*apps[i].getMaxBandwidth();
        }
    }
    
    while (cTime<endWin-ERROR and counttest<10)
    {
        //counttest++;
        if (cTime>timeChange-ERROR and not changedHeuristic){
            changedHeuristic=true;
            for (i=0;i<n;i++){
                initialProgress+=apps[i].getTProgress()*apps[i].getMaxBandwidth();
            }
            heuristic=heuristic2;
            start=cTime;
        }
        
        
        if constexpr (debuga>0)
        {
            cout << "\n\n*****NEW LOOP ! ("<<cTime<<")********\n";
        }
        if constexpr (debuga>1)
        {
            printIni<debuga>(cTime,apps);
        }

        // STEP A : Apply heuristic (Set the allocated bandwidth, initial value : 0)
        appInCom= {};
        for(i=0; i<n; i++)
        {
            if(apps[i].getType())
            {
                appInCom.push_back(i);
            }
        }
        nextDecision=compute(cTime,heuristic,appInCom,outName, apps, totBand, endWin,delta,failed);

        if constexpr (debuga>0)
        {
            printStepA(appInCom, apps);
        }

        // STEP B : Compute next event date (Update cTime)
        nextEvent=endWin;
        tempForTest=endWin-cTime;
        if(nextDecision<nextEvent)
        {
            nextEvent=nextDecision;
            tempForTest=nextDecision-cTime;
        }
        totalBW=0;
        totalPB=0;
        for (i=0; i<n; i++)
        {
            if (apps[i].getBW() <-ERROR or apps[i].getBW() > apps[i].getMaxBandwidth()+ERROR ){
                std::stringstream sb;
                sb << "FATAL ERROR : Illegal bandwidth allocated : " << apps[i].getBW() << " for maximum of " << apps[i].getMaxBandwidth();
                throw std::runtime_error(sb.str());
            }
            if (apps[i].getType())
            {
                /*if(apps[i].name=="ALS24"){
                    cout << apps[i].getBW() << " " << apps[i].getWiter() << " " << "\n";
                }*/
                endPhase=apps[i].getRemVol()/apps[i].getBW();
                totalPB+=apps[i].getBW();
                totalBW+=apps[i].getMaxBandwidth();
                if constexpr (debuga>1){
                    cout << "end "<< apps[i].name << " " << endPhase << "\n";
                }
                if(endPhase<tempForTest){
                    tempForTest=endPhase;
                }
            }
            else
            {
                endPhase=apps[i].getRemVol();
            }
            if(endPhase+cTime<nextEvent)
            {
                nextEvent=cTime+endPhase;
            }
        }
        
        if(totalPB/totBand<0.9998 and totalPB/totalBW<0.9998){
            std::stringstream sb;
            cout << "WARNING : You could allocate more bandwidth, which should never happen I think\n";
            cout << totalPB << " " << totBand-ERROR << " " << (totalPB<totBand-ERROR) << " " << (totalPB<totalBW-ERROR);
            //exit(0);aa
        }
        if constexpr (debuga>1){
            cout << "Total Bandwidth: " << totalPB << "\n";
            cout<< "Next known event in "<< tempForTest <<"\n";
            for(i=0;i<n;i++){
                if (apps[i].getType()){
                    newProgress=(tempForTest)*apps[i].getBW()/apps[i].getMaxBandwidth();
                    cout<<i << " " << newProgress << " " << (apps[i].getTProgress()+newProgress)/(tempForTest+cTime-apps[i].release)<<"\n";
                }
            }
        }
        
        if(cTime<timeChange-ERROR2 and nextEvent > timeChange+ERROR2){
            nextEvent=timeChange;
        }
        if constexpr (debuga>0)
        {
            cout<<"\n### STEP B : Next Event is at date : " << nextEvent << " ###\n\n";
        }

        // STEP C : Update all datas to prepare for next loop
        for(i=0; i<n; i++)
        {
            if(apps[i].getType())
            {
                newDVolume=(nextEvent-cTime)*apps[i].getBW();
                newProgress=(nextEvent-cTime)*apps[i].getBW()/apps[i].getMaxBandwidth();
            }
            else
            {
                newDVolume=nextEvent-cTime;
                newProgress=newDVolume;
            }
            apps[i].setTProgress(apps[i].getTProgress()+newProgress);
            newRVolume=apps[i].getRemVol()-newDVolume;
            if (apps[i].getType()){
                newRVolume2=newRVolume/apps[i].getMaxBandwidth();
            }
            else{
                newRVolume2=newRVolume;
            }
            if(newRVolume2<-100*ERROR2)
            {
                throw std::logic_error("FATAL ERROR : Next Event is After end of a phase");
            }
            else if(newRVolume2<100*ERROR2)
            {
                if constexpr (debuga>0)
                {
                    cout << "Application " << apps[i].name << " completed Phase "<<apps[i].getPhaseId() <<" of volume "<< apps[i].getCurAct()<<".\n";
                }
                if(not apps[i].getType() and cTime>timeChange){
                    if (cTime-apps[i].getCurAct()>timeChange){
                        apps[i].setUtilGopi(apps[i].getUtilGopi()+apps[i].getCurAct());
                    }
                    else{
                        apps[i].setUtilGopi(apps[i].getUtilGopi()+cTime-timeChange);
                    }
                }
                apps[i].setPhaseId(apps[i].getPhaseId()+1);
                apps[i].setRemVol(apps[i].getCurAct());
                apps[i].setReleasePhase(nextEvent);
                apps[i].setNew(true);
                apps[i].setType(not apps[i].getType());
            }
            else
            {
                apps[i].setRemVol(newRVolume);
            }
        }
        if(cTime == nextEvent){
            count+=1;
            if (cTime+ERROR2> endWin){
                break;
            }
            else if (count>10){
                throw std::logic_error("FATAL ERROR: cTime is not progressing");
            }
        }
        else{
            cTime=nextEvent;
            count=0;
        }
    }

    if constexpr (debuga>0){
        cout << "\n\n ### END !!! ###\n";
    }
    //GET MINIMUM YIELD
    if(cTime>endWin+100*ERROR or cTime < endWin-100*ERROR)
    {
        cout << endWin << " " << cTime<< " " << cTime-endWin << " " << ERROR <<  "\n";
        throw std::logic_error("FATAL ERROR: We went outside the window");
    }
    
    double bandi=0;
    double finalProgressGopi=0;
    
    for (i=0; i<n; i++)
    {
        if(not apps[i].getType()){
            if(cTime-apps[i].getCurAct()+apps[i].getRemVol()<start){
                apps[i].setUtilGopi(cTime-start);
            }
            else{
                apps[i].setUtilGopi(apps[i].getUtilGopi()+apps[i].getCurAct()-apps[i].getRemVol());
            }
        }
        finalProgressGopi+=apps[i].getUtilGopi();
        bandi+=apps[i].getMaxBandwidth();
        
        curyield=apps[i].getYield(cTime);
        if(curyield<minyield)
        {
            minyield=curyield;
        }
    }
    
    
    for (i=0;i<n;i++){
        finalProgress+=apps[i].getTProgress()*apps[i].getMaxBandwidth();
    }
    
    // Store it
    ofstream outFile;
    outFile.open(outName, ios::app);
    outFile << heuristic<<" " << minyield << " " << (finalProgress-initialProgress)/(cTime-start)/bandi << " " << finalProgressGopi/(cTime-start)/n << "\n";
    outFile.close();
    
    
    if (outdetails!=""){
        
        double nB=60.0;
        if(outdetails.find("Trilab")!= std::string::npos){
            nB=15.0;
        }
        pair<double,double> releases[n];
        pair<double,double> releases2[n];
        
        for(i=0;i<n;i++){
            releases[i]=make_pair(apps[i].realWiter,apps[i].getYield(cTime));
            releases2[i]=make_pair(apps[i].realRatio,apps[i].getYield(cTime));
        }
        size_t len = sizeof(releases) / sizeof(releases[0]);
        sort(releases, releases + len);
        len = sizeof(releases2) / sizeof(releases2[0]);
        sort(releases2, releases2 + len);
        
        ofstream outFile;
        outFile.open(outdetails+"-wIter.txt", ios::app);
        outFile << heuristic << " ";
        for(i=0;i<(unsigned long)(nB);i++){
            outFile << get<1>(releases[int(i/nB*n)])<<" ";
        }
        outFile << "\n";

        outFile.close();
        
        outFile.open(outdetails+"-Ratio.txt", ios::app);
        outFile << heuristic << " ";
        for(i=0;i<(unsigned long)(nB);i++){
            outFile << get<1>(releases2[int(i/nB*n)])<<" ";
        }
        outFile << "\n";

        outFile.close();
    }
    
    return minyield;
}

int main(int argc, char** argv)
{
    const int debuga=0; // Set to 3 : Full details,
    int parallelism;
    bool failed;
    // Set to 2 : Print timestep, state of all communicating applications at each decision point (+ the decisions),
    // Set to 1 : Only print the decisions and their timesteps
    // Set to 0 : Print nothing
    string aimedLoad;
    double scaleFactor;
    ifstream input(argv[1],ios::in);
    std::vector<std::tuple<string, double, double, string, string,string>> inputs;

    {
        string appname;
        string heuristic;
        string typeWin;
        string outName;
        double endWin;
        double totBand;
        while(input >> appname) {
            input >> totBand;
            input >> endWin;
            input >> heuristic;
            input >> typeWin;
            input >> outName;
            inputs.emplace_back(std::make_tuple(appname, totBand, endWin, heuristic, typeWin,outName));
        }
    }

#ifdef _OPENMP
    if constexpr (debuga > 0) 
        parallelism = 1;
    else {
        char *v = getenv("OMP_NUM_THREADS");
        if( NULL == v ) {
            std::cerr << "Set OMP_NUM_THREADS to an appropriate value to get parallelism. Defaulting to 1 thread." << std::endl;
            parallelism = 7;
        } else {
            parallelism = atoi(v);
            if(parallelism <= 0) {
                std::cerr << "Invalid value '" << v << "' for OMP_NUM_THREADS, defaulting to 1 thread." << std::endl;
                parallelism = 7;
            } else {
                std::cerr << "Running with " << parallelism << " threads. Set OMP_NUM_THREADS to an appropriate value to change this" << std::endl;
            }
        }
    }
#else
    parallelism = 6;
#endif

    #pragma omp parallel for num_threads(parallelism)
    for(unsigned int inp =0; inp<inputs.size();inp++)
    {
        auto run=inputs[inp];
        string appname;
        string heuristic;
        string typeWin;
        bool changeHeuristic;
        string outName;
        double endWin;
        double totBand;
        
        vector<Application> apps;
        appname = std::get<0>(run);
        totBand = std::get<1>(run);
        endWin  = std::get<2>(run);
        heuristic = std::get<3>(run);
        typeWin = std::get<4>(run);
        outName   = std::get<5>(run);
        
        if(typeWin=="v1"){
            changeHeuristic=false;
        }

        else{
            changeHeuristic=true;
        }
        
        if constexpr (debuga>0){
            cout<<"\n\n\n\n\n";
            cout << "\nNEW " << appname << " " << totBand << " " << endWin << " " <<heuristic << " " << changeHeuristic << "\n";
        }

        ifstream appFile;
        appFile.open(appname);
        
        addApps2(appFile, apps,endWin,totBand,heuristic); // just computes endWin
        appFile.close();
        appFile.open(appname);
        pair<double,double> data=addApps(appFile, apps,endWin,totBand,heuristic);
        appFile.close();
        double delta=get<0>(data);
        double rescaleBef= get<1>(data);
        string version;
        #pragma omp critical
        {
            cout << "Communication Load before AUTOSCALE : " << rescaleBef << std::endl;
        }
        
        totBand=1;
        
        unsigned first;
        unsigned last;
        string outdetails="";
        //cout << endWin<<"\n";
        if (AUTOSCALE){
            first = outName.find("v");
            last = outName.find(".txt");
            version = outName.substr (first+1,last-first-1);

            first = outName.find("w");
            last = outName.find("-");   
            aimedLoad = outName.substr (first+1,last-first-1);
            scaleFactor=stod(aimedLoad)/rescaleBef;
            autoScale(apps,scaleFactor);
            totBand/=scaleFactor;
            
            
            
            if(scaleFactor<1){
                endWin*=scaleFactor;
            }
            appFile.close();
            if (outName.find("-nH20-b0.5-s0.5") != std::string::npos){
                outdetails="results/detailedGopi/w"+aimedLoad+"v"+version;
            }
            if(outName.find("Nersc")!= std::string::npos or outName.find("nersc")!= std::string::npos){
                outdetails="results/detailedNersc/w"+aimedLoad+"v"+version;
            }
            if(outName.find("Trilab")!= std::string::npos or outName.find("trilab")!= std::string::npos){
                outdetails="results/detailedTrilab/w"+aimedLoad+"v"+version;
            }
                //cout << aimedLoad << " " << scaleFactor << " " << totBand << "\n";
        }
        
        try {
            failed=false;
            double minyield = simulate<debuga>(heuristic, totBand, endWin, apps, outName,changeHeuristic,delta, failed, outdetails);
            #pragma omp critical 
            {
                cout << appname << " " << totBand << " " << endWin << " " << heuristic << " " << typeWin << " -> MinYield of " << minyield << std::endl;
            }
            if(failed){
                ofstream outFile;
                outFile.open("error.txt", ios::app);
                outFile << " ERROR " << appname << " " << totBand << " " << endWin << " " << heuristic << " " << typeWin << "\n";
                outFile.close();
            }
        } catch (const std::exception &e) {
            #pragma omp critical
            {
                std::cerr << "Failed to simulate " << appname << " with heuristic " << heuristic << " and typeWin " << typeWin << " : " << e.what() << std::endl;
            }
        }
    }
    return 0;
}
