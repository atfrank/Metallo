//Author Lichirui Zhang
//A C++ program to read data from a txt/csv file into a 2-D vector
#include "read_csv.hpp"
#include<iostream>
#include<fstream>
#include<vector>
#include<sstream>
#include<string>

using namespace std;

csv_reader::csv_reader(const char* f){
    fp = fopen(f, "r");
    if (fp == NULL) {
        perror(f);
        exit(1);
    }
}

//destructor
csv_reader::csv_reader(){
    fclose(fp);
}

void csv_reader::read_csv(vector < vector<string> > &vec, const char *f, char delimiter){
    ifstream fp(f);
    string str;
    while(getline(fp, str)){
    	string word;
        vector <string> temp;
        
        for(auto c : str){
            if(isdigit(c) || isalpha(c)){
                word += c;
            }
            if(c == delimiter || c == ' ' || c == ','){ //set space and comma as default delimiters
                temp.push_back(word);
                word.clear();
            }
        }
        temp.push_back(word);
        word.clear();
        vec.push_back(temp);
        str.clear();
    }
}
