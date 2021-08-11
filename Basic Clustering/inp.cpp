#include <bits/stdc++.h>

using namespace std;

int main(){
	ofstream fout("input.txt");
	srand(time(NULL));
	for(int i = 0; i < 100; i++){
		for(int j = 0; j < 100; j++){
			fout << (((double)rand())/RAND_MAX) << " ";
		}
		fout << endl;
	}
	return 0;
}