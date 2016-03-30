#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;  
    
int main(void){
	vector<double> price;    
    ifstream myfile;
    int file = 1;
    string a = "../data/data";
    string b = ".csv";
    string name = a+to_string(file)+b;
    myfile.open(name);
    int count = 0;
    string line;
    while(myfile.good()){
        getline(myfile,line);
        int pos = line.find(",");
        if(pos>=3){
            price.push_back(stod(line.substr(pos+1)));
        }
        
    }
    cout<<price.size()<<endl;
    ofstream out;
    a = "../data/final";
    b = ".csv";
    name = a+to_string(file)+b;
    out.open(name);
    
    out<<price.size()<<endl;
    out.close();
    return 0;
}