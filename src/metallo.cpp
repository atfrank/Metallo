/*
 
 Copyright University of Michigan.
 This file is part of the Larmor software suite and is made available under license.
 University of Michigan (UM) TECHtransfer: phone: 734-763-0614 email: techtransfer@umich.edu.
 
 Author: Lichirui Zhang
 pending: improve efficiency of identifying remote grids
 
 */

#include "Molecule.hpp"
#include "Residue.hpp"
#include "Atom.hpp"
#include "Misc.hpp"
#include "Metallo.hpp"
#include "Analyze.hpp"
#include "Trajectory.hpp"
#include "AtomicFeaturizer.hpp"
#include "Chain.hpp"
#include "read_csv.hpp"

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <math.h>
#include <algorithm>
#include <stdlib.h>
#include <iomanip>

using namespace std;

void usage(){
    cerr << "====================================================" <<endl;
    cerr << "C++ Featurizer" << std::endl;
    cerr << "====================================================" <<endl;
    cerr << "Usage:   ./bin/test -spacing spacing -fc fc -row atomtype -sel atomtype -e elicf -treshold treshold -predictor predictor -dmyf -mgf -addedf <PDBfile>" << endl;
    cerr << "Options: [-spacing] Grid spacing (number) DEFAULT: 2.8Å" << endl;
    cerr << "         [-fc] featurization distance cutoff (number) DEFAULT: 8 Å" << endl;
    cerr << "         [-row] the atom type(s) for which fingerprints to be calculated (string) DEFAULT: None" << endl;
    cerr << "         [-sel] the atom type(s) used to describe the environment (string) DEFAULT: None"<< endl;
    cerr << "         [-predictor] use predictor 1 or 2. DEFAULT: 1" << endl;
    cerr << "         [-threshold] Place a Mg2+ for which the prediction is above the threshold. DEFAULT: 0.5" << endl;
    cerr << "         [-dmyf] Keep temporary file for adding dummy atoms if stated. DEFAULT: False" << endl;
    cerr << "         [-mgf] Keep temporary file for adding predicted Mg2+ if stated. DEFAULT: False" << endl;
    cerr << "         [-addedf] Write out the pdb file of the system after adding predicted Mg2+ if stated. DEFAULT: False" << endl;
    
    cerr << endl;
    exit(0);
}

int main(int argc, char **argv){
    
    vector<string> pdbs;
    string currArg;
    double spacing=2.8;
    double fc=8; // featurization distance cutoff
    string select_atm;
    string outfile = "./results/output.txt";
    double eta;
    int numEta;
    int etaStartPow;
    unsigned int predictor;
    
    double etaBase;
    vector<string> sel_atmname;
    vector<string> row_atmname;
    eta = 1.;
    numEta = 8;
    etaStartPow = -1;
    etaBase = 2.0;
    pdbs.clear();
    predictor = 1;
    double threshold = 0.50;
    bool keep_add_dummy = false;
    bool keep_add_mg = false;
    bool keep_addedmg_pdb = false;
    
    for (int i=1; i<argc; i++){
        currArg=argv[i];
        if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0){
            usage();
        }
        else if (currArg.compare("-fc") == 0){
            currArg = argv[++i];
            std::stringstream(currArg)>>fc;
        }

        else if (currArg.compare("-threshold") == 0 || currArg.compare("-t") == 0){
            currArg = argv[++i];
            std::stringstream(currArg)>>threshold;
        }
        else if (currArg.compare("-spacing") == 0 || currArg.compare("-s") == 0){
            currArg = argv[++i];
            std::stringstream(currArg) >> spacing;
        }
        else if (currArg.compare("-eta") == 0)
        {
            currArg=argv[++i];
            std::stringstream(currArg) >> eta;
        }
        else if (currArg.compare("-numEta") == 0)
        {
            currArg=argv[++i];
            std::stringstream(currArg) >> numEta;
        }
        else if (currArg.compare("-etaBase") == 0)
        {
            currArg=argv[++i];
            std::stringstream(currArg) >> etaBase;
        }
        else if (currArg.compare("-etaStartPow") == 0)
        {
            currArg=argv[++i];
            std::stringstream(currArg) >> etaStartPow;
        }
        else if (currArg.compare("-sel") == 0){
            currArg=argv[++i];
            select_atm = currArg;
            size_t pos = 0;
            string atm;
            string delim = " ";
            while ((pos = select_atm.find(delim)) != string::npos){
                atm = select_atm.substr(0, pos);
                sel_atmname.push_back(atm);
                select_atm.erase(0, pos + delim.length());
            }
            sel_atmname.push_back(select_atm);
        }
        else if (currArg.compare("-row") == 0){
            currArg=argv[++i];
            select_atm = currArg;
            size_t pos = 0;
            string atm;
            string delim = " ";
            while ((pos = select_atm.find(delim)) != string::npos){
                atm = select_atm.substr(0, pos);
                row_atmname.push_back(atm);
                select_atm.erase(0, pos + delim.length());
            }
            row_atmname.push_back(select_atm);
        }
        else if (currArg.compare("-output") == 0)
        {
            currArg=argv[++i];
            outfile = currArg;
        }
        else if (currArg.compare("-predictor") == 0)
        {
            currArg=argv[++i];
            predictor = atoi(currArg.c_str());
        }
        else if (currArg.compare("-dmyf") == 0)
        {
            keep_add_dummy = true;
        }
        else if (currArg.compare("-mgf") == 0)
        {
            keep_add_mg = true;
        }
        else if (currArg.compare("-addedf") == 0)
        {
            keep_addedmg_pdb = true;
        }
        else{
            pdbs.push_back(currArg);
        }
    }
    
    
    //Pop out warning if no pdb is provided
    if (pdbs.size() == 0 ){
        cerr<<endl<<"Error: Please provide an input PDB file"<<endl<<endl;
        usage();
    }
    
    
    //Read pdb
    for(unsigned int f=0; f<pdbs.size(); f++){
        Molecule *rna=NULL;
        rna = rna->readPDB(pdbs.at(f));
        cout<<endl<<"reading pdb complete."<<endl;
        
        string pdbfilename = pdbs.at(f);
        string pdbname = pdbfilename.substr(pdbfilename.length() - 8, 4);
        cout << pdbname << endl;
        
        
        //get Mg information
        //string Mg_info_filename1 = "./data/Mg_info/" + pdbname + ".txt";
        //const char *Mg_info_filename = &Mg_info_filename1[0];
        //csv_reader Mg_info_f(Mg_info_filename);
        
        //vector < vector<string> > Mg_info;
        //Mg_info_f.csv_reader::read_csv(Mg_info, Mg_info_filename);
        
        stringstream rnass;
    
        
        //for(unsigned int i=0; i<Mg_info.size(); ++i){
            //if(!count(NumericalChainID.begin(), NumericalChainID.end(), Mg_info[i][1])){
                //cout << Mg_info[i][2] << endl;
        //        rnass <<  ".:" << Mg_info[i][2] << ".MG_"; //A:1-5.CA_B:10-15.CA
            //}
        //}
        
        rnass<<".:.MG+C1'+C2+C2'+C3'+C4+C4'+C5+C5'+C6+C8+N1+N2+N3+N4+N6+N7+N9+O2+O2'+O3'+O4+O4'+O5'+O6+OP1+OP2+P+O"<<"&A+U+G+C+MG+HOH.:.";
        string rnasel = rnass.str();
        cout << rnasel <<endl;
        rna->select(rnasel);//segmentation fault when read 3v23
        rna = rna->copy(true);
        
        unsigned int natom=rna->getAtmVecSize();

        
        
        outfile = "./results/"+pdbname+".txt";

        //ofstream test;
        //test.open("./test_new.pdb");
        //test<<rna->writePDB(true,false,false);
        //test.close();
        
        cout<<"select rna atoms done! "<<natom<<" rna atoms in total"<<endl;
        
        double xmax = -9999;
        double ymax = -9999;
        double zmax = -9999;
        double xmin = 9999;
        double ymin = 9999;
        double zmin = 9999;
        
        //get the boundary of the system
        for (unsigned int i = 0; i<rna->getAtmVecSize(); i++){
            double xcoor = rna->getAtom(i)->getX();
            double ycoor = rna->getAtom(i)->getY();
            double zcoor = rna->getAtom(i)->getZ();
            
            if (xcoor>xmax){
                xmax = xcoor;
            }
            if (ycoor>ymax){
                ymax = ycoor;
            }
            if (zcoor>zmax){
                zmax = zcoor;
            }
            if (xcoor<xmin){
                xmin = xcoor;
            }
            if (ycoor<ymin){
                ymin = ycoor;
            }
            if (zcoor<zmin){
                zmin = zcoor;
            }
        }
        
        double xlength = fabs(xmax - xmin); // absolute value
        double ylength = fabs(ymax - ymin);
        double zlength = fabs(zmax - zmin);
        
        int nofxnode = ceil(xlength/spacing)+1;
        int nofynode = ceil(ylength/spacing)+1;
        int nofznode = ceil(zlength/spacing)+1;
        
        //cout << nofxnode * nofynode * nofznode << endl;
        
        vector <int> mgx;
        vector <int> mgy;
        vector <int> mgz;
        
        map<int, vector <float>> dist_dic;
        
        
        for (unsigned int n=0; n<natom; n++){
            double atmcoorx = rna->getAtom(n)->getX();
            double atmcoory = rna->getAtom(n)->getY();
            double atmcoorz = rna->getAtom(n)->getZ();
            
            int xtag = floor(fabs(atmcoorx - xmin)/spacing);
            int ytag = floor(fabs(atmcoory - ymin)/spacing);
            int ztag = floor(fabs(atmcoorz - zmin)/spacing);
            
            string atmid = rna->getAtom(n)->getAtmName(); //get atom name
            if(atmid == "MG"){
                mgx.push_back(xtag);
                mgy.push_back(ytag);
                mgz.push_back(ztag);
            }
        
        }
        
        // if(mgx.empty()){
        //     continue;
        // }
        
        vector <string> O_AtmName = {"O2", "O2'", "O3'", "O4'", "O5'", "O6'", "OP1'", "OP2'"}; //O2+O2'+O3'+O4+O4'+O5'+O6+OP1+OP2
        
        vector <string> rna_resname = {"A", "U", "G", "C"};
    
        
        //distance cutoff for N, OP1, OP2, others, and HOH
        double N_lower = 1.7; // lower bound, unit Å
        double N_higher = 7; // higher bound
        double OP1_lower = 1.6;
        double OP1_higher = 6.2;
        double OP2_lower = 1.5;
        double OP2_higher = 5.25;
        double others_lower = 1.65;
        double others_higher = 4.4;
        double hoh_lower = 1.525;
        double hoh_higher = 100000;
        
        //loop over atoms rather than grids
        for(unsigned int n=0; n<natom; n++){
            //information of the nth atom
            string name = rna->getAtom(n)->getAtmName();
            string resname = rna->getAtom(n)->getResName();
            double atmcoorx = rna->getAtom(n)->getX();
            double atmcoory = rna->getAtom(n)->getY();
            double atmcoorz = rna->getAtom(n)->getZ();
            
            //the nth atom is in (i,j,k)th grid
            int i = floor((rna->getAtom(n)->getX()-xmin)/spacing);
            int j = floor((rna->getAtom(n)->getY()-ymin)/spacing);
            int k = floor((rna->getAtom(n)->getZ()-zmin)/spacing);
            
            
            //for N, examine the grids around (i,j,k)th grid
            if(count(rna_resname.begin(), rna_resname.end(), resname) && name[0] == 'N'){
                int i_start = i - ceil(N_higher/spacing);
                int i_end = i + ceil(N_higher/spacing);
                int j_start = j - ceil(N_higher/spacing);
                int j_end = j + ceil(N_higher/spacing);
                int k_start = k - ceil(N_higher/spacing);
                int k_end = k + ceil(N_higher/spacing);
                
                if(i_start<0){
                    i_start = 0;
                }
                if(i_end>nofxnode){
                    i_end = nofxnode;
                }
                if(j_start<0){
                    j_start = 0;
                }
                if(j_end>nofynode){
                    j_end = nofynode;
                }
                if(k_start<0){
                    k_start = 0;
                }
                if(k_end>nofznode){
                    k_end = nofznode;
                }
                
                // candidate grids to examine.  x direction: i_start ~ i_end;  y direction: j-radius ~ j+radius;  z direction: k_start ~ k_end
                for(int x=i_start; x<=i_end; x++){
                    for(int y=j_start; y<=j_end; y++){
                        for(int z=k_start; z<=k_end; z++){
                            double gridcoorx = xmin+(x+1)*spacing - spacing/2;
                            double gridcoory = ymin+(y+1)*spacing - spacing/2;
                            double gridcoorz = zmin+(z+1)*spacing - spacing/2;
                            float distance = sqrt(pow(gridcoorx-atmcoorx,2)+pow(gridcoory-atmcoory,2)+pow(gridcoorz-atmcoorz,2)); //calculate distance between the nth atm and the (x,y,z)th grid center
                            int index = (x+1) + y*nofxnode + z*(nofxnode*nofynode);
                            // if it's the first time visit the grid
                            if(dist_dic.find(index)==dist_dic.end()){
                                // construct and store a distance list for Mg-N, OP1, OP2, others and HOH
                                vector <float> distlist = {distance, 9999, 9999, 9999, 9999}; // construct a distlist for the grid
                                dist_dic[index] = distlist;
                            }
                            // if the grid has been visited
                            else{
                                //update value
                                if(distance < dist_dic[index][0]){
                                    dist_dic[index][0] = distance;
                                }
                                else{
                                    // do nothing
                                }
                            }
                        }
                    }
                    
                }
            }
            
            // for OP1, do the same thing
            if(count(rna_resname.begin(), rna_resname.end(), resname) && name == "OP1"){
                int i_start = i - ceil(OP1_higher/spacing);
                int i_end = i + ceil(OP1_higher/spacing);
                int j_start = j - ceil(OP1_higher/spacing);
                int j_end = j + ceil(OP1_higher/spacing);
                int k_start = k - ceil(OP1_higher/spacing);
                int k_end = k + ceil(OP1_higher/spacing);
                
                if(i_start<0){
                    i_start = 0;
                }
                if(i_end>nofxnode){
                    i_end = nofxnode;
                }
                if(j_start<0){
                    j_start = 0;
                }
                if(j_end>nofynode){
                    j_end = nofynode;
                }
                if(k_start<0){
                    k_start = 0;
                }
                if(k_end>nofznode){
                    k_end = nofznode;
                }
                //grids to examine: for x grids, i-radius ~ i+radius; for y grids, j-radius ~ j+radius; for z grids, k-radius ~ k+radius
                for(int x=i_start; x<=i_end; x++){
                    for(int y=j_start; y<=j_end; y++){
                        for(int z=k_start; z<=k_end; z++){
                            double gridcoorx = xmin+(x+1)*spacing - spacing/2;
                            double gridcoory = ymin+(y+1)*spacing - spacing/2;
                            double gridcoorz = zmin+(z+1)*spacing - spacing/2;
                            float distance = sqrt(pow(gridcoorx-atmcoorx,2)+pow(gridcoory-atmcoory,2)+pow(gridcoorz-atmcoorz,2)); //calculate distance between the nth atm and the (x,y,z)th grid center
                            int index = (x+1) + y*nofxnode + z*(nofxnode*nofynode);
                           // if it's the first time visit the grid
                            if(dist_dic.find(index)==dist_dic.end()){
                                // construct and store a distance list for Mg-N, OP1, OP2, others and HOH
                                vector <float> distlist = {9999, distance, 9999, 9999, 9999}; // construct a distlist for the grid
                                dist_dic[index] = distlist;
                            }
                            // if the grid has been visited
                            else{
                                //update value
                                if(distance < dist_dic[index][1]){
                                    dist_dic[index][1] = distance;
                                }
                                else{
                                    // do nothing
                                }
                            }
                    }
                    
                }
                
            }
            }
            
            
            // for OP2, do the same thing
            if(count(rna_resname.begin(), rna_resname.end(), resname) && name == "OP2"){
                int i_start = i - ceil(OP2_higher/spacing);
                int i_end = i + ceil(OP2_higher/spacing);
                int j_start = j - ceil(OP2_higher/spacing);
                int j_end = j + ceil(OP2_higher/spacing);
                int k_start = k - ceil(OP2_higher/spacing);
                int k_end = k + ceil(OP2_higher/spacing);
                
                if(i_start<0){
                    i_start = 0;
                }
                if(i_end>nofxnode){
                    i_end = nofxnode;
                }
                if(j_start<0){
                    j_start = 0;
                }
                if(j_end>nofynode){
                    j_end = nofynode;
                }
                if(k_start<0){
                    k_start = 0;
                }
                if(k_end>nofznode){
                    k_end = nofznode;
                }
                //grids to examine: for x grids, i-radius ~ i+radius; for y grids, j-radius ~ j+radius; for z grids, k-radius ~ k+radius
                for(int x=i_start; x<=i_end; x++){
                    for(int y=j_start; y<=j_end; y++){
                        for(int z=k_start; z<=k_end; z++){
                            double gridcoorx = xmin+(x+1)*spacing - spacing/2;
                            double gridcoory = ymin+(y+1)*spacing - spacing/2;
                            double gridcoorz = zmin+(z+1)*spacing - spacing/2;
                            float distance = sqrt(pow(gridcoorx-atmcoorx,2)+pow(gridcoory-atmcoory,2)+pow(gridcoorz-atmcoorz,2)); //calculate distance between the nth atm and the (x,y,z)th grid center
                            int index = (x+1) + y*nofxnode + z*(nofxnode*nofynode);
                            // if it's the first time visit the grid
                            if(dist_dic.find(index)==dist_dic.end()){
                                // construct and store a distance list for Mg-N, OP1, OP2, others and HOH
                                vector <float> distlist = {9999, 9999, distance, 9999, 9999}; // construct a distlist for the grid
                                dist_dic[index] = distlist;
                            }
                            // if the grid has been visited
                            else{
                                //update value
                                if(distance < dist_dic[index][2]){
                                    dist_dic[index][2] = distance;
                                }
                                else{
                                    // do nothing
                                }
                            }
                    }
                    
                }
            }
            }
            
            
            
            // for other types of atoms, do the same thing
            if(name != "MG" && (name[0] != 'N' || name[0] != 'O')){
                int i_start = i - ceil(others_higher/spacing);
                int i_end = i + ceil(others_higher/spacing);
                int j_start = j - ceil(others_higher/spacing);
                int j_end = j + ceil(others_higher/spacing);
                int k_start = k - ceil(others_higher/spacing);
                int k_end = k + ceil(others_higher/spacing);
                
                if(i_start<0){
                    i_start = 0;
                }
                if(i_end>nofxnode){
                    i_end = nofxnode;
                }
                if(j_start<0){
                    j_start = 0;
                }
                if(j_end>nofynode){
                    j_end = nofynode;
                }
                if(k_start<0){
                    k_start = 0;
                }
                if(k_end>nofznode){
                    k_end = nofznode;
                }
                //grids to examine: for x grids, i-radius ~ i+radius; for y grids, j-radius ~ j+radius; for z grids, k-radius ~ k+radius
                for(int x=i_start; x<=i_end; x++){
                    for(int y=j_start; y<=j_end; y++){
                        for(int z=k_start; z<=k_end; z++){
                            double gridcoorx = xmin+(x+1)*spacing - spacing/2;
                            double gridcoory = ymin+(y+1)*spacing - spacing/2;
                            double gridcoorz = zmin+(z+1)*spacing - spacing/2;
                            float distance = sqrt(pow(gridcoorx-atmcoorx,2)+pow(gridcoory-atmcoory,2)+pow(gridcoorz-atmcoorz,2)); //calculate distance between the nth atm and the (x,y,z)th grid center
                            int index = (x+1) + y*nofxnode + z*(nofxnode*nofynode);
                             // if it's the first time visit the grid
                            if(dist_dic.find(index)==dist_dic.end()){
                                // construct and store a distance list for Mg-N, OP1, OP2, others and HOH
                                vector <float> distlist = {9999, 9999, 9999, distance, 9999}; // construct a distlist for the grid
                                dist_dic[index] = distlist;
                            }
                            // if the grid has been visited
                            else{
                                //update value
                                if(distance < dist_dic[index][3]){
                                    dist_dic[index][3] = distance;
                                }
                                else{
                                    // do nothing
                                }
                            }
                        }
                    }
                    
                }
                
            }
               
           // for O in HOH, do the same thing
           if(name == "O" && resname == "HOH"){
               int i_start = i - ceil(hoh_higher/spacing);
               int i_end = i + ceil(hoh_higher/spacing);
               int j_start = j - ceil(hoh_higher/spacing);
               int j_end = j + ceil(hoh_higher/spacing);
               int k_start = k - ceil(hoh_higher/spacing);
               int k_end = k + ceil(hoh_higher/spacing);
               
               if(i_start<0){
                   i_start = 0;
               }
               if(i_end>nofxnode){
                   i_end = nofxnode;
               }
               if(j_start<0){
                   j_start = 0;
               }
               if(j_end>nofynode){
                   j_end = nofynode;
               }
               if(k_start<0){
                   k_start = 0;
               }
               if(k_end>nofznode){
                   k_end = nofznode;
               }
               //grids to examine: for x grids, i-radius ~ i+radius; for y grids, j-radius ~ j+radius; for z grids, k-radius ~ k+radius
               for(int x=i_start; x<=i_end; x++){
                   for(int y=j_start; y<=j_end; y++){
                       for(int z=k_start; z<=k_end; z++){
                           double gridcoorx = xmin+(x+1)*spacing - spacing/2;
                           double gridcoory = ymin+(y+1)*spacing - spacing/2;
                           double gridcoorz = zmin+(z+1)*spacing - spacing/2;
                           float distance = sqrt(pow(gridcoorx-atmcoorx,2)+pow(gridcoory-atmcoory,2)+pow(gridcoorz-atmcoorz,2)); //calculate distance between the nth atm and the (x,y,z)th grid center
                           int index = (x+1) + y*nofxnode + z*(nofxnode*nofynode);
                           // if it's the first time visit the grid
                           if(dist_dic.find(index)==dist_dic.end()){
                               // construct and store a distance list for Mg-N, OP1, OP2, others and HOH
                               vector <float> distlist = {9999, 9999, 9999, 9999, distance}; // construct a distlist for the grid
                               dist_dic[index] = distlist;
                           }
                           // if the grid has been visited
                           else{
                               //update value
                               if(distance < dist_dic[index][4]){
                                   dist_dic[index][4] = distance;
                               }
                               else{
                                   // do nothing
                               }
                           }
                       }
                   }
                   
               }
           }
        }//bracket for looping over atms
           
        ofstream fout;
        //add dummy atoms
        //J is placed in qualified empty grids
		int atmid = 1;
        fout.open(("./addatom/addatom_"+pdbname+".pdb").c_str());
        for (int i=0; i<nofxnode; i++){
            for (int j=0; j<nofynode; j++){
                for(int k=0; k<nofznode; k++){
                    int index = (i+1) + j*nofxnode + k*(nofxnode*nofynode);
                    bool condition = false;
                    // if the grid hasn't been visited
                    if(dist_dic.find(index)==dist_dic.end()){
                        continue; // examine next one
                    }
                    else{
                        if(dist_dic[index][0]<=N_higher && dist_dic[index][0]>=N_lower && dist_dic[index][1]<=OP1_higher && dist_dic[index][1]>=OP1_lower && dist_dic[index][2]<=OP2_higher && dist_dic[index][2]>=OP2_lower && dist_dic[index][3]<=others_higher && dist_dic[index][3]>=others_lower && dist_dic[index][4]<=hoh_higher && dist_dic[index][4]>=hoh_lower && dist_dic[index][2] <= dist_dic[index][1] + 2.7 && dist_dic[index][2] >= dist_dic[index][1] - 2.7){
                        condition = true;
                    }
                    }
    
                    if(condition){ //only add dummy atoms in empty grids
                        bool isMg = false;
						for(unsigned int mi=0; mi<mgx.size(); mi++){
							if (i == mgx[mi] && j == mgy[mi] && k == mgz[mi]){
                              fout<<"ATOM"<<right<<setw(7)<<atmid<<right<<setw(7)<<"L    L"<<"     "<<left<<setw(4)<<2<<right<<setw(11)<<setiosflags(ios::fixed)<<setprecision(3)<<xmin+(i+1)*spacing - spacing/2<<right<<setw(8)<<ymin+(j+1)*spacing - spacing/2<<right<<setw(8)<<zmin+(k+1)*spacing - spacing/2<<"  1.00  0.00           L"<<endl;//important to output in this format otherwise addatom.pdb won't be read properly
							  isMg = true;
							  mgx.erase(mgx.begin() + mi);
							  mgy.erase(mgy.begin() + mi);
							  mgz.erase(mgz.begin() + mi);
							  break;
							}
						}
						if (! isMg){
						  fout<<"ATOM"<<right<<setw(7)<<atmid<<right<<setw(7)<<"J    J"<<"     "<<left<<setw(4)<<1<<right<<setw(11)<<setiosflags(ios::fixed)<<setprecision(3)<<xmin+(i+1)*spacing - spacing/2<<right<<setw(8)<<ymin+(j+1)*spacing - spacing/2<<right<<setw(8)<<zmin+(k+1)*spacing - spacing/2<<"  1.00  0.00           J"<<endl;//important to output in this format otherwise addatom.pdb won't be read properly
						}
						atmid += 1;
                    }
                }
            }
        }
        // L is placed in Mg2+ occupied grids
        for(unsigned int i=0; i<mgx.size(); i++){
           fout<<"ATOM"<<right<<setw(7)<<atmid<<right<<setw(7)<<"M    M"<<"     "<<left<<setw(4)<<2<<right<<setw(11)<<setiosflags(ios::fixed)<<setprecision(3)<<xmin+(mgx[i]+1)*spacing - spacing/2<<right<<setw(8)<<ymin+(mgy[i]+1)*spacing - spacing/2<<right<<setw(8)<<zmin+(mgz[i]+1)*spacing - spacing/2<<"  1.00  0.00           M"<<endl;//important to output in this format otherwise addatom.pdb won't be read properly
		   atmid += 1;
        }
        
        Molecule *DummyAtms = NULL;
        Atom *a, *atmEntry = new Atom;
        
        Chain *c, *chnEntry;
        Residue *res, *resEntry;
        
        c=NULL;
        res=NULL;
        
        chnEntry=NULL;
        resEntry=NULL;
        
        DummyAtms = DummyAtms->readPDB("./addatom/addatom_"+pdbname+".pdb");
        //cout << " aaa " <<endl;
        
        for (unsigned int i=0; i< DummyAtms->getChnVecSize(); i++){
            c=DummyAtms->getChain(i);
            chnEntry=new Chain;
            
            for (unsigned int j=0; j< c->getResVecSize(); j++){
                res=c->getResidue(j);
                resEntry=new Residue;
                
                for (unsigned int k=0; k< res->getAtmVecSize(); k++){
                    a=res->getAtom(k);
                    //Add each selected atom
                    atmEntry=new Atom;
                    atmEntry->clone(a);
                    rna->addAtom(atmEntry);
                    resEntry->addAtom(atmEntry);
                    chnEntry->addAtom(atmEntry);
                }
                if (resEntry->getAtmVecSize() > 0){
                    rna->addResidue(resEntry);
                    chnEntry->addResidue(resEntry);
                }
                else{
                    delete resEntry;
                }
            }
            if (chnEntry->getResVecSize() > 0){
                rna->addChain(chnEntry);
            }
            else{
                delete chnEntry;
            }
        }
        
        fout.close();
        
        
        //write out dummy atoms added rna structures
        ofstream dmyaddedrna;
        dmyaddedrna.open(("./proba/dmyaddedrna_"+pdbname+".pdb").c_str());
        dmyaddedrna<<rna->writePDB(true,false,false);
        dmyaddedrna.close();
        
        
        stringstream ss;
        
        ss << ".:.J#" << fc+1 << "_.:.L#" << fc+1 << "_.:.M#" << fc+1;
        rna->select(ss.str());
        // rna->select(".:.J_.:.L_.:.M");
		rna = rna->copy(true);
         
        //test.open("./test_featurization.pdb");
        //test<<rna->writePDB(true,false,false);
        //test.close();
         
         
        
        
        //featurization
        cout << "start featurizing" << endl;
        vector<double> etalist;
        vector<vector<double> > features;
        vector<double> feature_tmp;
        feature_tmp.clear();
        
        //create list of eta value
        for (int i = 0; i < numEta; i ++){
            eta = pow(etaBase, i + etaStartPow);
            etalist.push_back(eta);
        }
        
        AtomicFeaturizer *atomicfeature = new AtomicFeaturizer(rna);
        //Metallo *metallo = new Metallo(predictor);
        
        //save space for featurization
        dist_dic.clear();
        
        string addmgfilename="./addatom/add_mg_"+pdbname+".pdb";
        char*p=(char*)addmgfilename.data();
        remove(p);
        
        //ofstream proba;
        //proba.open(("./proba/"+pdbs.at(f).substr(0,4)+".txt").c_str());
        
        atomicfeature->featurizeScalar(fc, etalist, features, outfile, sel_atmname, row_atmname);
        //save space for featurization
        delete rna;
        for (unsigned int i=0; i<features[0].size(); i++){
            //cout << i  << " " << features[i].size() << " " <<endl;
            feature_tmp.clear();
            feature_tmp.resize(features.size());
            for (unsigned int j=0; j<features.size(); j++){
                //cout << i  << " " << features.size() << " " <<endl;
                feature_tmp.at(j) = features[j][i];
            }
            //metallo->Predict(feature_tmp);
            // add logic to take prediction and place a Mg2+ for which the prediction is above a certain threshold
            
            //  if(metallo->Predict(feature_tmp)>threshold ){
            //metallo->addmg(rna, DummyAtms, pdbname, i);
            // }
            
            // cout << i  << " " << metallo->Predict(feature_tmp) << " " <<endl;
            //proba<< i  << " " << metallo->Predict(feature_tmp) << " " <<endl;
        }
        
        string adddmyatmfile = "./addatom/addatom_"+pdbname+".pdb";
        char*q=(char*)adddmyatmfile.data();
        
        if(keep_add_mg == false){
            remove(p);
        }
        if(keep_add_dummy == false){
            remove(q);
        }
        if(keep_addedmg_pdb == true){
            ofstream addedmgpdb;
            addedmgpdb.open(("./addatom/mg_predicted_"+pdbs.at(f).substr(0,4)+".pdb").c_str());
            addedmgpdb<<rna->writePDB();
            addedmgpdb.close();
        }
        
        features.clear();
        if(atomicfeature!=NULL){
            delete atomicfeature;
        }
        delete chnEntry;
        delete resEntry;
        delete a;
        delete atmEntry;
    }
    return 0;
};
