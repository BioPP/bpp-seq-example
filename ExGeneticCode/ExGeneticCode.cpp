/*
 * File: ExGeneticCode.cpp
 * Created by: Julien Dutheil
 * Created on: Jun Mon 10 21:00 2013
 * Last modified: Jun Mon 10 21:00 2013
 *
 * An example code which displays a given genetic code.
 */

/*----------------------------------------------------------------------------------------------------*/

/*
 * We start by including what we'll need, and sort the inclusions a bit:
 */

/*
 * From the STL:
 */
#include <iostream> /* to be able to output stuff in the terminal. */

/*
 * We'll use the standard template library namespace:
 */
using namespace std;

/*
 * From bpp-seq:
 */
#include <Bpp/Seq/Alphabet.all>
#include <Bpp/Seq/GeneticCode.all>

/*
 * All Bio++ functions are also in a namespace, so we'll use it:
 */
using namespace bpp;

/*----------------------------------------------------------------------------------------------------*/
/*
 * Now starts the real stuff...
 */


int main(int args, char ** argv)
{
  try {
    DNA*           nt = new DNA();
    CodonAlphabet* ca = new CodonAlphabet(nt);
    //GeneticCode*   gc = new StandardGeneticCode(nt); 
    //GeneticCode*   gc = new VertebrateMitochondrialGeneticCode(nt); 
    //GeneticCode*   gc = new YeastMitochondrialGeneticCode(nt); 
    //GeneticCode*   gc = new InvertebrateMitochondrialGeneticCode(nt); 
    //GeneticCode*   gc = new EchinodermMitochondrialGeneticCode(nt); 
    //GeneticCode*   gc = new AscidianMitochondrialGeneticCode(nt); 
    GeneticCode*   gc = new MoldMitochondrialGeneticCode(nt); 
    
    //To keep the same order as NCBI:
    vector<int> x(4);
    x[0] = 3; x[1] = 1; x[2] = 0; x[3] = 2;

    cout << "--------------+---------------+---------------+--------------" << endl;
    for (int pos1 = 0; pos1 < 4; ++pos1) {
      for (int pos3 = 0; pos3 < 4; ++pos3) {
        for (int pos2 = 0; pos2 < 4; ++pos2) {
          string s = nt->intToChar(x[pos1]) + nt->intToChar(x[pos2]) + nt->intToChar(x[pos3]);
          int i = ca->charToInt(s);
          string aa = "";
          if (gc->isStop(s))
            aa = "STP";
          else {
            int j = gc->translate(i);
            aa = gc->getTargetAlphabet()->getAbbr(j);
          }
          if (pos2 > 0)
            cout << "| ";
          cout << s << " (" << (i < 10 ? "0" : "") << i << "): " << aa;
          if (gc->isStart(i))
            cout << "#";
          else if (gc->isAltStart(i))
            cout << "*";
          else
            cout << " ";
        }
        cout << endl;
      }
      cout << "--------------+---------------+---------------+--------------" << endl;
    }

  } catch(Exception& e) {
    cerr << e.what() << endl;
  }

  return 0;
}

