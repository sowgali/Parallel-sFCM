#include <bits/stdc++.h>

using namespace std;

double square(double x){
	return x*x;
}

int main(){
	double myImg[104][104];
	ifstream fin("input.txt");
	for(int i = 0; i < 104; i++){
		if(i < 2 || i > 101){
			for(int j = 0; j < 104; j++){
				myImg[i][j] = 0;
			}
		}
		else{
			string params;
			getline(fin, params);
			istringstream ss(params);
			myImg[i][0] = 0;
			myImg[i][1] = 0;
			for(int j = 2; j < 102; j++)
				ss >> myImg[i][j];
			myImg[i][102] = 0;
			myImg[i][103] = 0;
		}
	}
	cout << "No. of clusters:" << endl;
	int K;
	cin >> K;
	double membership_matrix[100][100][K];
	double spatialH[100][100][K];
	srand(time(NULL));
	double old_means[K];
	double new_means[K];
	cout << "Initial Means:" << endl;
	for(int i = 0; i < K; i++){
		old_means[i] = (((double)rand())/RAND_MAX);
		new_means[K] = 0;
		cout << old_means[i] << endl;
	}
	bool flag = true;
	int iter = 0;
	while(flag){
		for(int i = 0; i < 100; i++){
			for(int j = 0; j < 100; j++){
				for(int k = 0; k < K; k++){
					double temp = 0;
					for(int l = 0; l < K; l++)
						temp += square((myImg[i+2][j+2] - old_means[k])/(myImg[i+2][j+2] - old_means[l] + 1e-5));
					membership_matrix[i][j][k] = 1.0/(temp + 1e-5);
				}
			}
		}
		/*for(int i = 2; i < 102; i++){
			for(int j = 2; j < 102; j++){
				for(int k = 0; k < K; k++){
					double temp = 0;
					for(int l = i - 2; l <= i + 2; l++){
						for(int m = j - 2; m <= j + 2; m++){
							//temp += membership_matrix[l][m]
						}
					}
				}
			}
		}*/
		for(int k = 0; k < K; k++){
			double num = 0;
			double den = 0;
			for(int i = 0; i < 100; i++){
				for(int j = 0; j < 100; j++){
					num += square(membership_matrix[i][j][k]) * myImg[i+2][j+2];
					den += square(membership_matrix[i][j][k]);
				}
			}
			new_means[k] = num/(den+1e-5);
		}
		bool temp = true;
		for(int i = 0; i < K; i++){
			if(abs(new_means[i] - old_means[i]) <= 0.001)
				temp = temp && true;
			else
				temp = temp && false;
			old_means[i] = new_means[i];
			new_means[i] = 0;
		}
		flag = !temp;
		iter++;
	}
	cout << "The means are:" << endl;
	for(int i = 0; i < K; i++){
		cout << old_means[i] << endl;
	}
	ofstream fout("KMeans.txt");
	for(int i = 0; i < 100; i++){
		for(int j = 0; j < 100; j++){
			int cluster = 0;
			for(int k = 0; k < K; k++)
				if(membership_matrix[i][j][cluster] < membership_matrix[i][j][k])
					cluster = k;
			fout << cluster+1 << " ";
		}
		fout << endl;
	}
	return 0;
}