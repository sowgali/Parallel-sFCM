#include <bits/stdc++.h>

using namespace std;

double** myImg;
double* oldMeans;
double* newMeans;
atomic<bool>*** updating;
atomic<bool>*** updated;
double*** membershipMatix;
double*** spatialMatix;
ifstream fin("sample.txt");
int w,h,p,q,K;
bool truth = true;
bool fallacy = false;
double square(double x){
	return x*x;
}

void sfcm(int x_start, int x_end, int y_start, int y_end){
	for(int i = x_start; i <= x_end; i++){
		for(int j = y_start; j <= y_end; j++){
			for(int k = 0; k < K; k++){
				if(updating[i][j][k].compare_exchange_strong(fallacy, truth)){
					if(!updated[i][j][k]){
						double temp = 0;
						for(int l = 0; l < K; l++)
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
	for(int i = x_start; i <= x_end; i++){
		for(int j = y_start; j <= y_end; j++){
			for(int k = 0; k < K; k++){
				double temp = 0;
				for(int b1 = max(i-2,0); b1 < min(i+3, h); b1++){
					for(int b2 = max(j-2, 0); b2 < min(j+3,w); b2++){
						if(updated[b1][b2][k].load(memory_order_seq_cst))
							temp += membershipMatix[b1][b2][k];
						else if(updating[b1][b2][k].compare_exchange_strong(fallacy, truth)){
							double temp = 0;
							for(int l = 0; l < K; l++)
								temp += square((myImg[b1][b2] - oldMeans[k])/(myImg[b1][b2] - oldMeans[l] + 1e-5));
							membershipMatix[b1][b2][k] = 1.0/(temp + 1e-5);
							updating[b1][b2][k].store(false, memory_order_seq_cst);
							updated[b1][b2][k].store(true, memory_order_seq_cst);
						}
						else{
							while(!updated[b1][b2][k]);
							temp += membershipMatix[b1][b2][k];
						}
					}
				}
				spatialMatix[i][j][k] = temp;
			}
		}
	}
	for(int i = x_start; i <= x_end; i++){
		for(int j = y_start; j <= y_end; j++){
			double sum = 0;
			for(int k = 0; k < K; k++)
				sum += pow(membershipMatix[i][j][k],p) * pow(spatialMatix[i][j][k],q);
			for(int k = 0; k < K; k++)
				membershipMatix[i][j][k] = (pow(membershipMatix[i][j][k],p) * pow(spatialMatix[i][j][k],q))/(sum + 1e-5);
		}
	}
}

void update(int idx){
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

bool checkConv(double eps){
	bool temp = true;
	for(int i = 0; i < K; i++){
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
    cin >> K;
	cout << "Threshold for convergence: ";
	cin >> epsilon;
	cout << "p and q: ";
	cin >> p >> q;
    oldMeans = new double[K];
    newMeans = new double[K];
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
			membershipMatix[i][j] = new double[K];
			spatialMatix[i][j] = new double[K];
			updating[i][j] = new atomic<bool>[K];
			updated[i][j] = new atomic<bool>[K];
			for(int k = 0; k < K; k++){
				membershipMatix[i][j][k] = 0;
				spatialMatix[i][j][k] = 0;
				updated[i][j][k].store(false, memory_order_seq_cst);
				updating[i][j][k].store(false, memory_order_seq_cst);
			}
		}
    }
	cout << "Initial means:" << endl;
    for(int i = 0; i < K; i++){
        oldMeans[i] = (((double)rand())/RAND_MAX);
		cout << oldMeans[i] << endl;
        newMeans[i] = 0;
    }
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
		for(int k = 0; k < K; k++)
			tids[k] = thread(update, k);
		for(auto& t : tids)
			t.join();
	}
	cout << "New means:" << endl;
    for(int i = 0; i < K; i++)
		cout << oldMeans[i] << endl;
	fin.close();
	return 0;
}