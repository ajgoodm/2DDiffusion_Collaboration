// August 2016
// Program for simulating 2D Diffusion and annihilation

#include <iostream>
#include <vector>		// For dynamically allocated vectors
#include <string>		// For file naming
#include <fstream>		// For saving results to file for analysis in Matlab
#include <cmath>		// For Sine and Cosine
#include <cstdlib>

using namespace std;


int annihilate(double x[][3], int Nexcitons, double radius, vector < vector<double> > &hist, int tindex, double dhist, int NBINS);
												/* looks at all of the excitons stored in array x[][3]
												and annihilates an exciton from pairs of excitons
												separated by less than 2*radius. Returns the number of
												surviving excitons. Passes the 2D vector array hist by reference
												and if you have taken n*100 time steps stores the instantaneous
												distribution of pairwise distances between excitons that fall in (R.,1 um]
												histogrammed in hist. dhist is 1um/NBINS */

int removeZeros(double x[][3], int Nstart);
												/* takes the 2D array x[][3] which records the exciton's (x,y) coordinates
												in x=x[i][0] and y=x[i][1] and whether or not the exciton exists bool=x[i][2].
												It then sorts x by the third column so all of the remaining excitons are the
												first rows of x so later calculations can ignore the dead excitons at the
												last rows of x. Returns the number of living excitons. */

bool isOverlapped(double x[][3], double radius, int nthExciton);
												/* takes 2D array and checks the exciton indexed by nthExciton
												to see if it overlaps with the other excitons already placed.
												true means the excitons	overlaps a previous exciton.
												false means it does not. */

int main()
{												// declare variables for user input.
	const int NXSTART=30;						// exciton density per square micron
	double Diffusivity;							// diffusivity in um^2/s
	double Radius=0;							// radius in um
	double timens;								// time in ns

	string file1;								// declare filenames to write output to.
	string file2;
	string file3;

	cout << "What diffusivity would you like to simulate [um^2/s]? ";
	cin >> Diffusivity;
	cout << "What radius would you like to simulate [um]? ";
	cin >> Radius;
	cout << "How much time  would you like to simulate [s]? ";
	cin >> timens;

	/*------------------------------------ constructing file names -------------*/
	string Nfstr="N_";
	string Hfstr="H_";
	string Dsp="D1e";
	string Dstr=to_string(log10(Diffusivity)).substr(0,1);
	string Rsp="_R";
	string Rstr=to_string(Radius).substr(2,3);
	string Nsp="_Nstart";
	string Nstr=to_string(NXSTART);
	string Ifstr="Info_";
	string txt=".txt";

	file1=Nfstr+Dsp+Dstr+Rsp+Rstr+Nsp+Nstr+txt;
	file2=Hfstr+Dsp+Dstr+Rsp+Rstr+Nsp+Nstr+txt;
	file3=Ifstr+Dsp+Dstr+Rsp+Rstr+Nsp+Nstr+txt;
	/*------------------------------------ file names constructed -------------*/


	srand(time(NULL));
 	ofstream outputFile;
	ofstream histFile;
	ofstream infofile;

	outputFile.open(file1);
	histFile.open(file2);
	infofile.open(file3);

	double D=Diffusivity;				// Diffusivity in um^2/s
	double stepsize=0.001;				// step size in um (1 nm). Chosen to be much smaller
										// than exciton radius (roughly 5 nm).
	double krad=2e8;					// radiative decay rate in s^-1.
	double dt=pow(stepsize,2)/D;		// time increment is determined by the diffusivity
										// and the hop size.
	double dt1=timens/5000;				// maximum dt for meshing small diffusivity cases
	dt=min(dt,dt1);						// recompute dt so we calculate at least 50,000 time
										// points in our time period
	stepsize=sqrt(dt*D);				// recompute stepsize for meshing


	/*----------------- variable declarations ------------*/
	int Nexciton;
	const long int tsteps=floor(timens/dt);	// number of time steps
	const int NTRIALS=100;					// number of simulation trials
	const double RADIUS=Radius;				// annihilation radius in um;
	const int BOXSIDE=2;					// simulation box edge size in um
	double halfboxside=BOXSIDE/2.0;
	const double pi = 3.1415926535897;
	const int NBINS=1000;					// number of bins to mesh 1um with for
											// histogram
	double dhist=halfboxside/double(NBINS); // spacing between histogram bins.
	int attempt;
	double random;
	long int ii;							// generic loop index 1
	long int jj;							// generic loop index 2

	double x[NXSTART*BOXSIDE*BOXSIDE][3];	// x[0][i] is the x position of exciton i,
											// x[1][i] is the y position of exciton i,
											// x[2][i] is {0,1} saying whether exciton
											// i has decayed
	double histbins[NBINS]={};				// histbins is vector of histogram bin positions



	vector< vector<double> > hist;
	hist.resize(NBINS, vector<double>(static_cast<int>(floor(tsteps/25))));
	vector<double> Nsave(tsteps);

	/* hist is the histogram of pairwise distance distributions at variuous times
	2D array. Declared as vector for memory management. Nsave is the population
	decay vector saved over the multiple trials. */

	//Info file is a file that stores the the time increment and histogram bin positions
	//for data analysis purposes.
	infofile << "dt = " << dt << endl;
	for (ii=0; ii<NBINS; ii++)
	{
		histbins[ii]=ii*dhist;
		infofile << ii*dhist << ", ";
	}
	infofile << endl;
	/*------------------------------*/


	for (int trial=0; trial<NTRIALS; trial++)
	{
		cout << trial << endl;
		Nexciton=NXSTART*BOXSIDE*BOXSIDE;				// number of excitons at start
		for (ii=0; ii<Nexciton; ii++)						// initialize exciton positions in
		{																				// BOXSIDE x BOXSIDE square. density is 1/Nexciton
			attempt=0;
			do
			{
				x[ii][2]=1;
				random=BOXSIDE*double(rand())/double(RAND_MAX);
				x[ii][0]=random;
				random=BOXSIDE*double(rand())/double(RAND_MAX);
				x[ii][1]=random;
				attempt++;
			}while (isOverlapped(x,RADIUS,ii)&&(attempt<200));

		}


		for (ii=0; ii<tsteps; ii++)											// ii is looping over time index
		{
			for (jj=0; jj<Nexciton; jj++)									// updating exciton positions after hop
			{																// jj is looping over exciton index x[jj][]
				random=double(2*pi*double(rand())/double(RAND_MAX));
				x[jj][0]=x[jj][0]+(cos(random)*stepsize);
				x[jj][1]=x[jj][1]+(sin(random)*stepsize);

				if (x[jj][0]<0)
				{x[jj][0]=x[jj][0]+BOXSIDE;}
				if (x[jj][0]>BOXSIDE)
				{x[jj][0]=x[jj][0]-BOXSIDE;}
				if (x[jj][1]<0)
				{x[jj][1]=x[jj][1]+BOXSIDE;}
				if (x[jj][1]>BOXSIDE)
				{x[jj][1]=x[jj][1]-BOXSIDE;}
				
				// Radiative decay condition
				random=double(rand())/double(RAND_MAX);
				if (random>exp(-dt*krad/2))
					x[jj][2]=0;

			}

			Nexciton=annihilate(x,Nexciton,RADIUS,hist,ii,dhist,NBINS);		// tracks exciton population after every time step
			Nsave[ii]=Nsave[ii]+Nexciton;									// tracks exciton population aggregated over all trials
		}
	}

	/*------------------------
	Simulation over. Entering file output.
	------------------------*/

	for(ii=0; ii<tsteps; ii++)
	{
		outputFile << Nsave[ii]/static_cast<double>(BOXSIDE*BOXSIDE*NTRIALS) << endl;
	}


	for (ii=0; ii<NBINS; ii++)
	{
		for (jj=1; jj<floor(tsteps/25); jj++)
		{
			histFile << ", " << hist[ii][jj];
		}
		histFile << endl;
	}


	outputFile.close();
	histFile.close();
	infofile.close();

	return 0;
}



int annihilate(double x[][3], int Nexcitons, double radius, vector < vector<double> > &hist, int tindex, double dhist, int NBINS)
{
	int ii;
	int jj;
	double dist;


	for (ii=0; ii<Nexcitons; ii++)
	{
		for (jj=ii+1; jj<Nexcitons; jj++)
		{
			dist=sqrt(pow((x[ii][0]-x[jj][0]),2.0)+pow((x[ii][1]-x[jj][1]),2.0)); 	//Calculate distance between exciton ii and jj


			if ((tindex%25==0)&&(int(floor(dist/dhist))<NBINS)&&(dist>2*radius))	//If we've taken n*100 tstepts append hist.
			{																	  	//increment value in histogram bin by 1.
				hist[int(floor(dist/dhist))][tindex/25]=hist[int(floor(dist/dhist))][tindex/25]+1;
			}

			if (dist<=2*radius)														//If excitons are within radius of each other
			{																		//kill an exciton.
				x[jj][2]=0;
			}
		}
	}

	return removeZeros(x, Nexcitons);												//sorts x so living excitons are at the front
}																					//checks and returns the number of living excitons.

int removeZeros(double x[][3], int Nstart)
{
	int lastIndex=0;
	double temp;

	for (int counter=0; counter<Nstart; counter++)									//all rows beyond Nstart are dead excitons (x[ii][2]==0 for
	{																				//all ii>=Nstart.

		if (x[counter][2]==1.0)														//if you encounter a living exciton in row counter
		{																			//swap it with the lowest row containing a dead exciton
																					//which is stored as lastIndex. Increment last index by 1
			temp=x[counter][0];														//sorts all living excitons to to the top of x.
			x[counter][0]=lastIndex;
			x[lastIndex][0]=temp;
			temp=x[counter][1];
			x[counter][1]=lastIndex;
			x[lastIndex][1]=temp;
			temp=x[counter][2];
			x[counter][2]=lastIndex;
			x[lastIndex][2]=temp;
			lastIndex=lastIndex+1;
		}
	}

	return lastIndex++;																//the last index+1 (because array is zero-indexed) is
}																					//the number of living excitons.

bool isOverlapped(double x[][3], double radius, int nthExciton)
{
	for (int ii=0; ii<nthExciton-1; ii++)
	{
		if (sqrt(pow((x[ii][0]-x[nthExciton][0]),2)+pow((x[ii][1]-x[nthExciton][1]),2))<(2*radius))
		{return true;}
	}
	return false;
}
