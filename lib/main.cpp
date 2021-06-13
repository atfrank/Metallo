//Author Lichirui Zhang

#include "read_csv.h"
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

int main(){
    
    /*
    const char *filename1 = "./data/systems.txt";
    const char *filename2 = "./data/benchmark.txt";
    
    csv_reader f1(filename1);
    csv_reader f2(filename2);
    vector < vector<string> > vec1;
    vector < vector<string> > vec2;
    
    f1.csv_reader::read_csv(vec1, filename1);
    f2.csv_reader::read_csv(vec2, filename2);
    //cout << "aaa" << endl;
    
    
    
    
    for(int i=0; i < vec1.size(); i++){
        string pdbname = vec1[i][0];
        for(int j=0; j < vec2.size(); j++){
            if(vec2[j][0] == pdbname){
                ofstream file;
                file.open("./data/Mg_info/"+pdbname+".txt", ios_base::app);
                file << pdbname << " " << vec2[j][1] << " " << vec2[j][2] << " \n";
            }
        }
    }
    */
    

    const char *filename = "benchmark.txt";
    
    csv_reader f(filename);
    
    vector < vector<string> > vec;
    
    f.csv_reader::read_csv(vec, filename);
    
    for(int i=0; i < vec.size(); i++){
        cout << vec[i][2] << endl;
    }
    
    
    /*
    const char *f = "benchmark.txt";

    ifstream fp(f);
    string str;
    while(getline(fp, str)){
        cout << str << endl;
        str.clear();
    }
    */
    
    

    return 0;
}

