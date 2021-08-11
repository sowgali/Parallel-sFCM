#include <bits/stdc++.h>

using namespace std;
using namespace std::chrono;

double** myImg;					//input image
double* oldMeans;				//means for each cluster in the previous iteration
double* newMeans;				//means for each cluster after current iteration
atomic<bool>*** updating;		//indicates the status of a data point whether it is being updated by some thread in the current iteration
atomic<bool>*** updated;		//indicates the status of a data point whether it is already updated by some thread or not in the current iteration
double*** membershipMatix;		//stores the membership values of each datapoint to each of the cluster
double*** spatialMatix;			
ifstream fin("sample.txt");
int w,h,p,q,K,C;
bool truth = true;
bool fallacy = false;
double square(double x){
	return x*x;
}

void sfcm(int x_start, int x_end, int y_start, int y_end){			//calculates membership matrix for specified region by also considering the spatial information
	for(int i = x_start; i <= x_end; i++){							//find the membership matrix without considering the spatial information
		for(int j = y_start; j <= y_end; j++){
			for(int k = 0; k < C; k++){
				if(updating[i][j][k].compare_exchange_strong(fallacy, truth)){
					if(!updated[i][j][k]){
						double temp = 0;
						for(int l = 0; l < C; l++)
							temp += square((myImg[i][j] - oldMeans[k])/(myImg[i][j] - oldMeans[l] + 1e-5));
						membershipMatix[i][j][k] = 1.0/(temp + 1e-5);
						updated[i][j][k].store(true, memory_order_seq_cst);
						updating[i][j][k].store(false, memory_order_seq_cst);
					}
					else{
						updating[i][j][k].store(false, memory_order_seq_cst);
					}
				}
			}
		}
	}
	for(int i = x_start; i <= x_end; i++){							//calculates the spatial matrix
		for(int j = y_start; j <= y_end; j++){
			for(int k = 0; k < C; k++){
				double temp = 0;
				for(int b1 = max(i-2,0); b1 < min(i+3, h); b1++){
					for(int b2 = max(j-2, 0); b2 < min(j+3,w); b2++){
						if(updated[b1][b2][k].load(memory_order_seq_cst))
							temp += membershipMatix[b1][b2][k];
						else if(updating[b1][b2][k].compare_exchange_strong(fallacy, truth)){
							double temp = 0;
							for(int l = 0; l < C; l++)
								temp += square((myImg[b1][b2] - oldMeans[k])/(myImg[b1][b2] - oldMeans[l] + 1e-5));
							membershipMatix[b1][b2][k] = 1.0/(temp + 1e-5);
							updated[b1][b2][k].store(true, memory_order_seq_cst);
							updating[b1][b2][k].store(false, memory_order_seq_cst);
						}
						else{
							while(!updated[b1][b2][k]){cout << "Hello" << endl;}
							temp += membershipMatix[b1][b2][k];
						}
					}
				}
				spatialMatix[i][j][k] = temp;
			}
		}
	}
	for(int i = x_start; i <= x_end; i++){						// calculates the membership matrix with spatial information using previously calculated 
		for(int j = y_start; j <= y_end; j++){					// membership matrix and spatial matrix
			double sum = 0;
			for(int k = 0; k < C; k++)
				sum += pow(membershipMatix[i][j][k],p) * pow(spatialMatix[i][j][k],q);
			for(int k = 0; k < C; k++)
				membershipMatix[i][j][k] = (pow(membershipMatix[i][j][k],p) * pow(spatialMatix[i][j][k],q))/(sum + 1e-5);
		}
	}
}

void update(int idx){											//calculates final means for all clusters 
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

bool checkConv(double eps){								//checks if the algorithm converged
	bool temp = true;
	for(int i = 0; i < C; i++){
		if(abs(newMeans[i] - oldMeans[i]) > eps)
			temp = false;
		oldMeans[i] = newMeans[i];
		newMeans[i] = 0;
	}
	return temp;
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
    oldMeans = new double[C];
    newMeans = new double[C];
    myImg = new double*[h];
    membershipMatix = new double**[h];
	spatialMatix = new double**[h];
	updating = new atomic<bool>**[h];
	updated = new atomic<bool>**[h];
    for(int i = 0; i < h; i++){
        myImg[i] = new double[w];
		string params;
		getline(fin, params);
		istringstream ss(params);
		membershipMatix[i] = new double*[w];
		spatialMatix[i] = new double*[w];
		updating[i] = new atomic<bool>*[w];
		updated[i] = new atomic<bool>*[w];
		for(int j = 0; j < w; j++){
			ss >> myImg[i][j];
			membershipMatix[i][j] = new double[C];
			spatialMatix[i][j] = new double[C];
			updating[i][j] = new atomic<bool>[C];
			updated[i][j] = new atomic<bool>[C];
			for(int k = 0; k < C; k++){
				membershipMatix[i][j][k] = 0;
				spatialMatix[i][j][k] = 0;
				updated[i][j][k].store(false, memory_order_seq_cst);
				updating[i][j][k].store(false, memory_order_seq_cst);
			}
		}
    }
	int randrows[h] = {0};
	int randcols[w] = {0};
	iota(randrows,randrows+h,0);
	iota(randcols,randcols+w,0);
	random_shuffle(randrows,randrows+h);
	random_shuffle(randcols,randcols+w);
	cout << "Initial means:" << endl;
    for(int i = 0; i < C; i++){
        oldMeans[i] = myImg[randrows[i]][randcols[i]];
		cout << oldMeans[i] << endl;
        newMeans[i] = 0;
    }
	auto tstart = high_resolution_clock::now();
	int iter = 0;
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
		for(int k = 0; k < K; k++){
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
			tids[k] = thread(sfcm, xstart, xend, ystart, yend);
			xstart = xend+1;
			ystart = yend+1;
		}
		for(auto& t : tids)
			t.join();
		for(int i=0;i<h;++i){												//make the updated and updating matrices false for next iteration
			for(int j=0;j<w;++j){
				for(int k=0;k<C;++k){
					updated[i][j][k].store(false,memory_order_seq_cst);
					updating[i][j][k].store(false,memory_order_seq_cst);
				}
			}
		}
        thread tid2[C];
		for(int k = 0; k < C; k++)
			tid2[k] = thread(update, k);
		for(auto& t : tid2)
			t.join();
		//iter++;
		cout << iter++ << endl;
	}
	auto tend = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(tend-tstart);
	cout<<"time: "<<duration.count()<<endl;
	cout<<"iters: "<<iter<<endl;
	cout << "New means:" << endl;
    for(int i = 0; i < C; i++)
		cout << oldMeans[i] << endl;
	fin.close();
	return 0;
}