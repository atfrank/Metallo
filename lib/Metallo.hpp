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

#include <string>
#include <vector>
#include <map>
#include "Molecule.hpp"
#include "Residue.hpp"
#include "Atom.hpp"
#include "Misc.hpp"
#include "Analyze.hpp"
#include "Trajectory.hpp"
#include "AtomicFeaturizer.hpp"
#include "Chain.hpp"


using namespace std;

//#include "Molecule.hpp"
//#include "DTree.hpp"

/* Forward Declaration  (only valid for pointers and references) */
class Molecule;
class DTree;
class Metallo_RF;

class Metallo {
  private:
    std::vector<DTree *> trees;
    unsigned int predictor;
  public:
    Metallo(unsigned int predictor);
    void Load();
    void SetPrediction(double prediction);
    double Predict(const std::vector<double> features);
    bool addmg(Molecule *mol, Molecule *DummyAtms, string pdbname, unsigned int i);
};

