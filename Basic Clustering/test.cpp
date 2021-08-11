#include <bits/stdc++.h>

using namespace std;

int main(){
    ifstream fin("test.txt");
    string params;
    getline(fin, params);
    istringstream ss(params);
    while(ss){
        string temp;
        ss >> temp;
        cout << temp << endl;
    }
    return 0;
}