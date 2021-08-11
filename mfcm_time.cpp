#include <bits/stdc++.h>

using namespace std;
using namespace std::chrono;

double** myImg;
double* oldMeans;
double* newMeans;
double*** membershipMatix;
double*** spatialMatix;
ifstream fin("sample.txt");
int w,h,p,q,K,C;

double square(double x){
	return x*x;
}

bool checkConv(double eps){				// checks if the convergence is reached 
	bool temp = true;
	for(int i = 0; i < K; i++){
		if(abs(newMeans[i] - oldMeans[i]) > eps)
			temp = false;
		oldMeans[i] = newMeans[i];
		newMeans[i] = 0;
	}
	return temp;
}

void fcm(int x_start, int x_end, int y_start, int y_end){			// calculates initial membership Matrix without considering spatial information
	for(int i = x_start; i <= x_end; i++){							
		for(int j = y_start; j <= y_end; j++){
			for(int k = 0; k < C; k++){
				double temp = 0;
				for(int l = 0; l < C; l++)
					temp += square((myImg[i][j] - oldMeans[k])/(myImg[i][j] - oldMeans[l] + 1e-5));
				membershipMatix[i][j][k] = 1.0/(temp + 1e-5);
			}
		}
	}
}

void sfcm(int x_start, int x_end, int y_start, int y_end){			//calculates membership Matrix after considering the spatial information using 
	for(int i = x_start; i <= x_end; i++){						    //using previously calculated membership matrix
		for(int j = y_start; j <= y_end; j++){
			for(int k = 0; k < C; k++){
				double temp = 0;
				for(int b1 = max(i-2,0); b1 < min(i+3, h); b1++){
					for(int b2 = max(j-2, 0); b2 < min(j+3,w); b2++)
						temp += membershipMatix[b1][b2][k];
				}
				spatialMatix[i][j][k] = temp;
			}
		}
	}
	for(int i = x_start; i <= x_end; i++){
		for(int j = y_start; j <= y_end; j++){
			double sum = 0;
			for(int k = 0; k < C; k++)
				sum += pow(membershipMatix[i][j][k],p) * pow(spatialMatix[i][j][k],q);
			for(int k = 0; k < K; k++)
				membershipMatix[i][j][k] = (pow(membershipMatix[i][j][k],p) * pow(spatialMatix[i][j][k],q))/(sum + 1e-5);
		}
	}
}

void update(int idx){										// calculates new mean for a given cluster
	double num = 0;
	double den = 0;
	for(int i = 0; i < h; i++){
		for(int j = 0; j < w; j++){
			num += square(membershipMatix[i][j][idx]) * myImg[i][j];
			den += square(membershipMatix[i][j][idx]);
		}
	}
	newMeans[idx] = num/(den+1e-5);
}

int main(){
	double epsilon;
    srand(time(NULL));
    cout << "Width and Height: ";
    cin >> w >> h;
    cout << "No. of clusters: ";
    cin >> C;
    cout << "No. of threads: ";
    cin >> K;
	cout << "Threshold for convergence: ";
	cin >> epsilon;
	cout << "p and q: ";
	cin >> p >> q;
    oldMeans = new double[C];							// to store means calculated in the previous iteration
    newMeans = new double[C];							// to store means calculated in the current iteration
    myImg = new double*[h];
    membershipMatix = new double**[h];
	spatialMatix = new double**[h];
    for(int i = 0; i < h; i++){									// initializing myImg, membershipMatrix, spatialMatrix
        myImg[i] = new double[w];
		string params;
		getline(fin, params);
		istringstream ss(params);
		membershipMatix[i] = new double*[w];
		spatialMatix[i] = new double*[w];
		for(int j = 0; j < w; j++){
			ss >> myImg[i][j];
			membershipMatix[i][j] = new double[C];
			spatialMatix[i][j] = new double[C];
			for(int k = 0; k < C; k++){
				membershipMatix[i][j][k] = 0;
				spatialMatix[i][j][k] = 0;
			}
		}
    }
	int randrows[h] = {0};								// used for selecting C random and distinct data points from myImg matrix
	int randcols[w] = {0};
	iota(randrows,randrows+h,0);
	iota(randcols,randcols+w,0);
	random_shuffle(randrows,randrows+h);
	random_shuffle(randcols,randcols+w);
	cout << "Initial means:" << endl;
    for(int i = 0; i < C; i++){									//randomly select C distinct datapoints as initial means 
        oldMeans[i] = myImg[randrows[i]][randcols[i]];
		cout << oldMeans[i] << endl;
        newMeans[i] = 0;
    }
	int iter = 0;
	auto tstart = high_resolution_clock::now();					// timer starts
	while(!checkConv(epsilon)){
		thread tids[K];
		int wstep = w/K;
		int hstep = h/K;
		int wleft = w%K;
		int hleft = h%K;
		int xstart = 0;
		int ystart = 0;
		int xend = 0;
		int yend = 0;
		for(int k = 0; k < K; k++){								// each thread calculates the membership matrix(w/o spatial information) for its own submatrix
			xend = xstart + hstep;
			yend = ystart + wstep;
			if(hleft <= 0){
				xend--;
			}
			if(wleft <= 0){
				yend--;
			}
			hleft--;
			wleft--;
			tids[k] = thread(fcm, xstart, xend, ystart, yend);
			xstart = xend+1;
			ystart = yend+1;
		}
		for(auto& t : tids)										
			t.join();
		xstart = 0;
		ystart = 0;
		xend = 0;
		yend = 0;
		for(int k = 0; k < K; k++){								// each thread calculates the membership matrix (spatial information taken into account) using
			xend = xstart + hstep;								// previously calculated membership matrix
			yend = ystart + wstep;
			if(hleft <= 0){
				xend--;
			}
			if(wleft <= 0){
				yend--;
			}
			hleft--;
			wleft--;
			tids[k] = thread(sfcm, xstart, xend, ystart, yend);
			xstart = xend+1;
			ystart = yend+1;
		}
		for(auto& t : tids)
			t.join();
        thread tid2[C];										// allocate a thread for each of the cluster
		for(int k = 0; k < C; k++)							// find the final means for all the clusters
			tid2[k] = thread(update, k);
		for(auto& t : tid2)
			t.join();
		iter++;
	}
	auto tend = high_resolution_clock::now();						// timer ends
	auto duration = duration_cast<microseconds>(tend-tstart);		// time taken
	cout<<"time : "<<duration.count()<<endl;	
	cout << iter << endl;
	cout << "New means:" << endl;
    for(int i = 0; i < C; i++)
		cout << oldMeans[i] << endl;
	fin.close();
	return 0;
}