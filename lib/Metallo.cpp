//Aaron T. Frank
//Sean M. Law

    
/*
This file is part of MoleTools.

MoleTools is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MoleTools is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MoleTools.  If not, see <http://www.gnu.org/licenses/>.
*/

//Code generated using: awk '{print "{\""$1":"$2":"$3"\","$4"},"}' larmorD_both.dat | tr '\n' ' '

#include "Metallo.hpp"
#include "Metallo_RF.hpp"
#include "Molecule.hpp"
#include "Misc.hpp"
#include "DTree.hpp"


#include <fstream>
#include <stdlib.h>
#include <iomanip>
#include <iostream>
#include <cstdio>




Metallo::Metallo(unsigned int predictor){
    /* initialize classs */
    this->Load();
    this->predictor = predictor;
}

void Metallo::Load(){
  unsigned int ntree;
  std::vector<std::string> tokens;
  ntree = Metallo_RF::getNTree();
  trees.clear();
  trees.resize(ntree);
  tokens.clear();
    

  for (unsigned int i=0; i< ntree; i++){
    Misc::splitStr(Metallo_RF::getTree(i, this->predictor), " \t", tokens, false);
    trees.at(i) = new DTree;
    trees.at(i)->genDTree(tokens, ":");
  }
}

double Metallo::Predict(const std::vector<double> features){
  /* stores features */
  double p1;
  double ptmp;
  unsigned int ntree;
  std::string s;
    
  // Metallo prediction
  // this prediction scheme is only valid for a binary (2-class) predictor
  p1 = 0.0;
  ptmp = 0.0;
  ntree=Metallo_RF::getNTree();
  for (unsigned int i=0; i< ntree; i++){
    ptmp = atof(trees.at(i)->getDTreeClass(features).c_str());
    if(ptmp != 0.0){
        p1 += 1.0;
    } else {
        p1 += 0.0;
    }
  }
  p1 /= ntree;
  return p1;
}

//add dummy atoms
bool Metallo::addmg(Molecule *mol, Molecule *DummyAtms, string pdbname, unsigned int i){
    Molecule *addatom=NULL;
    Molecule *newaddatom = NULL;
    //string addmgfilename="./addatom/add_mg_"+pdbname+".pdb";
    //char*p=(char*)addmgfilename.data();
    
    ofstream addfile;
    addfile.open("./addatom/add_mg_"+pdbname+".pdb");
    
    if(mol->getAtmVecSize()+i<10000){
    addfile<<"HETATM"<<right<<setw(5)<<mol->getAtmVecSize()+i<<"  M     M M"<<right<<setw(4)<<mol->getAtmVecSize()+i<<right<<setw(12)<<setiosflags(ios::fixed)<<setprecision(3)<<DummyAtms->getAtom(i)->getX()<<right<<setw(8)<<DummyAtms->getAtom(i)->getY()<<right<<setw(8)<<DummyAtms->getAtom(i)->getZ()<<"  1.00  0.00           M"<<endl;
    }
    
    else if(mol->getAtmVecSize()+i > 9999 && mol->getAtmVecSize()+i<99999){
        addfile<<"HETATM"<<right<<setw(5)<<mol->getAtmVecSize()+i<<"  M     M M"<<right<<setw(4)<<mol->getAtmVecSize()+i<<right<<setw(11)<<setiosflags(ios::fixed)<<setprecision(3)<<DummyAtms->getAtom(i)->getX()<<right<<setw(8)<<DummyAtms->getAtom(i)->getY()<<right<<setw(8)<<DummyAtms->getAtom(i)->getZ()<<"  1.00  0.00           M"<<endl;
    }
    
    else if(mol->getAtmVecSize()+i > 99999){
        addfile<<"HETATM"<<right<<setw(5)<<mol->getAtmVecSize()+i<<"  M     M M"<<right<<setw(4)<<mol->getAtmVecSize()+i<<right<<setw(10)<<setiosflags(ios::fixed)<<setprecision(3)<<DummyAtms->getAtom(i)->getX()<<right<<setw(8)<<DummyAtms->getAtom(i)->getY()<<right<<setw(8)<<DummyAtms->getAtom(i)->getZ()<<"  1.00  0.00           M"<<endl;
        
    }
    
    addfile.close();

    

    
    
    addatom = newaddatom->readPDB("./addatom/add_mg_"+pdbname+".pdb");
    
    Atom *a, *atmEntry = new Atom;
    
    Chain *c, *chnEntry;
    Residue *res, *resEntry;
    
    c=NULL;
    res=NULL;
    
    chnEntry=NULL;
    resEntry=NULL;
    
    for (unsigned int i=0; i< addatom->getChnVecSize(); i++){
        c=addatom->getChain(i);
        chnEntry=new Chain;
        
        for (unsigned int j=0; j< c->getResVecSize(); j++){
            res=c->getResidue(j);
            resEntry=new Residue;
            
            for (unsigned int k=0; k< res->getAtmVecSize(); k++){
                a=res->getAtom(k);
                //Add each selected atom
                atmEntry=new Atom;
                atmEntry->clone(a);
                mol->addAtom(atmEntry);
                resEntry->addAtom(atmEntry);
                chnEntry->addAtom(atmEntry);
            }
            if (resEntry->getAtmVecSize() > 0){
                mol->addResidue(resEntry);
                chnEntry->addResidue(resEntry);
            }
            else{
                delete resEntry;
            }
        }
        if (chnEntry->getResVecSize() > 0){
            mol->addChain(chnEntry);
        }
        else{
            delete chnEntry;
        }
    }
    return 0;

}
