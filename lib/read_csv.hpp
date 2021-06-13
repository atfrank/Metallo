//Author Lichirui Zhang
//Author Lichirui Zhang
//A C++ program to read data from a txt/csv file into a 2-D vector
#include<iostream>
#include<fstream>
#include<vector>
#include<sstream>

using namespace std;

class csv_reader {
public:
    csv_reader(const char*);
    csv_reader();
    
    void read_csv(vector < vector<string> > &vec, const char *f, char delimiter = ' ');
    
private:
    csv_reader(const csv_reader&);
    csv_reader& operator=(const csv_reader&);
    FILE* fp;
};