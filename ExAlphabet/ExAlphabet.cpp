/*
 * File: ExAlphabet.cpp
 * Created by: Julien Dutheil
 * Created on: Dec Mon 01 13:33 2008
 * 
 * Introduction to the Alphabet class, and a bit of propaganda on good habits to take!
 *
 * HOW TO USE THAT FILE:
 * - General comments are written using the * * syntax.
 * - Code lines are switched off using '//'. To activate those lines, just remove the '//' characters!
 * - You're welcome to extensively modify that file!
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
#include <Bpp/Seq/Alphabet.all> /* this include all alphabets in one shot */

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
  /*
   * We surround our code with a try-catch block, in case some error occurs:
   */
  try
  {
    cout << "Hello World!" << endl;

    Alphabet* alphabet = new DNA();
    cout << "Playing with alphabet " << alphabet->getAlphabetType() << endl;

    string state = "A";
    bool test = alphabet->isCharInAlphabet("A");
    cout << "Test for state " << state << ": " << (test ? ":)" : ":(") << endl;
    if(test) cout << state << " actually stands for " << alphabet->getName(state) << endl;
    
    cout << "Numbers of characters: " << alphabet->getNumberOfChars() << endl;
    cout << "Numbers of types     : " << alphabet->getNumberOfTypes() << endl;
    cout << "Size of the alphabet : " << alphabet->getSize()          << endl;

    /*
     * Alphabet link character states (human representation) to integer codes (computer representation).
     *
     * Although you'll probably never have to care about that (you're human, not a computer after all),
     * it might be good to sneak a bit in here... try the following piece of code:
     */
    vector<string> states = alphabet->getSupportedChars();
    cout << "Char\tInt" << endl;
    for (unsigned int i = 0; i < states.size(); i++)
    {
      cout << states[i] << "\t" << alphabet->charToInt(states[i]) << endl; 
    }

  }
  catch(Exception& e)
  {
    cout << "Bio++ exception:" << endl;
    cout << e.what() << endl;
    return 1;
  }
  catch(exception& e)
  {
    cout << "Any other exception:" << endl;
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}

