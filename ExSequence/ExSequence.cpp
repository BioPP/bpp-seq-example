/*
 * File: ExSequence.cpp
 * Created by: Julien Dutheil
 * Created on: Dec Tue 03 17:13 2008
 * Last modified: Jun Thu 04 07:09 2009
 * 
 * Introduction to the Sequence class.
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
#include <iostream> //to be able to output stuff in the terminal.

/*
 * We'll use the standard template library namespace:
 */
using namespace std;

/*
 * From Bpp-Seq:
 */
#include <Bpp/Seq/Alphabet.all> /* this include all alphabets in one shot */
#include <Bpp/Seq/Sequence.h> /* this include the definition of the Sequence object */
#include <Bpp/Seq/SequenceTools.h> /* this include some tool sto deal with sequences */
#include <Bpp/Seq/DNAToRNA.h> /* A few translators here... */
#include <Bpp/Seq/NucleicAcidsReplication.h>
#include <Bpp/Seq/GeneticCode/StandardGeneticCode.h>

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

    /*
     * We first create a sequence object from scratch, specifying the alphabet:
     */
    Sequence* sequence = new BasicSequence("My first sequence", "GATTACAATGATTACATGGT", &AlphabetTools::DNA_ALPHABET);
    cout << sequence->getName() << endl;
    cout << sequence->toString() << endl;
    cout << sequence->size() << " positions." << endl;

    /*
     * Accessing the single positions:
     */
    for (unsigned int i = 0; i < sequence->size(); i++)
    {
      cout << sequence->getChar(i) << endl;
    }

    //Writes a dot plot:
    for (unsigned int i = 0; i < sequence->size(); i++) {
      for (unsigned int j = 0; j < sequence->size(); j++) {
        if (sequence->getChar(i) == sequence->getChar(j))
          cout << sequence->getChar(i);
        else
          cout << ".";
      }
      cout << endl;
    }

    /*
     * Sequence objects are derived from a more general structure called SymbolList.
     * Nothing really important here, you just have to know that the Site objects, that we'll
     * meet later, also inherit from SymbolList. Hence all methods presented here will
     * also work with Site objects.
     */
    Sequence* seq1 = sequence->clone();
    Sequence* seq2 = sequence->clone();
    seq2->deleteElement(4);
    Sequence* seq3 = sequence->clone();
    seq3->addElement(4, "T");
    Sequence* seq4 = sequence->clone();
    seq4->setElement(4, "A");

    cout << "Seq1: " << seq1->toString() << endl;
    cout << "Seq2: " << seq2->toString() << endl;
    cout << "Seq3: " << seq3->toString() << endl;
    cout << "Seq4: " << seq4->toString() << endl;

    /*
     * A bit of cleaning...
     */
    delete seq1;
    delete seq2;
    delete seq3;
    delete seq4;
    
    /*
     * Now we'll see what we can do with sequences...
     *
     * Several methods are available in the SequenceTools static class.
     * 'static' means that this class can be used without any instance, you can just call
     * directly any of their methods. In bio++, there are several {*}Tools classes which are all static,
     * and deal with particular data structures. This includes SiteTools and CodonSiteTools for instance.
     * These three classes inherit from the general SymbolListTools class, and add more specialized functions.
     */
    Sequence* subSequence = SequenceTools::subseq(*sequence, 6, 8);
    cout << "SubSeq: " << subSequence->toString() << endl;
    delete subSequence;

    Sequence* reverseSequence = 0;
    try {
      reverseSequence = SequenceTools::reverseTranscript(*sequence);
    } catch(Exception& e) {
      cerr << e.what() << endl;
    }
    
    Sequence* transSequence = SequenceTools::transcript(*sequence);
    cout << "TransSeq: " << transSequence->toString() << endl;

    reverseSequence = SequenceTools::reverseTranscript(*transSequence);
    cout << "RevSeq: " << reverseSequence->toString() << endl;
    delete transSequence;
    delete reverseSequence;

    /*
     * We'll now do a bit of /in silico/ molecular biology.
     * Basically, this means decoding/recoding sequences according to given alphabets.
     *
     * First we need to understand how sequences are coded.
     * A sequence (or a site) is coded as a vector of int codes, and the correspondance between 
     * int code and the actual character string is ensured by the Alphabet object (see previous exercise).
     * So when you create a Sequence object from a string, you are actually *parsing* its content.
     * Try the following:
     */
    cout << "This sequence is coded with a " << sequence->getAlphabet()->getAlphabetType() << endl;
    for (unsigned int i = 0; i < sequence->size(); i++)
    {
      cout << sequence->getChar(i) << "\t" << sequence->getValue(i) << "\t" << (*sequence)[i] << endl;
    }

    /*
     * To change the Alphabet of a sequence, you need to decode and recode it:
     */
    try {
      Sequence* protSequence = new BasicSequence(sequence->getName(), sequence->toString(), &AlphabetTools::PROTEIN_ALPHABET);
      cout << "This sequence is now coded with a " << protSequence->getAlphabet()->getAlphabetType() << endl;
      delete protSequence;
    } catch (Exception& ex) {
      cerr << ex.what() << endl;
    }

    try {
      Sequence* rnaSequence = new BasicSequence(sequence->getName(), sequence->toString(), &AlphabetTools::RNA_ALPHABET);
      cout << "This sequence is now coded with a " << rnaSequence->getAlphabet()->getAlphabetType() << endl;
      delete rnaSequence;
    } catch (Exception& ex) {
      cerr << ex.what() << endl;
    }

    CodonAlphabet* codonAlphabet = new StandardCodonAlphabet(&AlphabetTools::DNA_ALPHABET);
    Sequence* codonSequence = new BasicSequence(sequence->getName(), sequence->toString(), codonAlphabet);
    cout << "This sequence is now coded with a " << codonSequence->getAlphabet()->getAlphabetType() << endl;
    for (unsigned int i = 0; i < codonSequence->size(); i++)
    {
      cout << codonSequence->getChar(i) << "\t" << codonSequence->getValue(i) << endl;
    }
    delete codonSequence;
    delete codonAlphabet;

    /*
     * To make more complexe deparsing/reparsing, you need *Translator* objects.
     * These objects allow you to convert from an alphabet to another in a very general way.
     *
     * Example 1: changing DNA to RNA:
     */

    Transliterator* transliterator = new DNAToRNA();
    Sequence* trSequence = transliterator->translate(*sequence);
    cout << "Original  : " << sequence->toString() << endl;
    cout << "Translated: " << trSequence->toString() << endl;
    delete trSequence;
    delete transliterator;

    /*
     * Example 2: getting the complement sequence, in the same alphabet:
     */
    transliterator = new NucleicAcidsReplication(&AlphabetTools::DNA_ALPHABET, &AlphabetTools::DNA_ALPHABET);
    trSequence = transliterator->translate(*sequence);
    cout << "Original  : " << sequence->toString() << endl;
    cout << "Translated: " << trSequence->toString() << endl;
    delete trSequence;
    delete transliterator;
    
    /*
     * Example 3: The same but with an RNA complement:
     */
    transliterator = new NucleicAcidsReplication(&AlphabetTools::DNA_ALPHABET, &AlphabetTools::RNA_ALPHABET);
    trSequence = transliterator->translate(*sequence);
    cout << "Original  : " << sequence->toString() << endl;
    cout << "Translated: " << trSequence->toString() << endl;
    delete trSequence;
    delete transliterator;
    
  }
  catch (Exception& e)
  {
    cout << "Bio++ exception:" << endl;
    cout << e.what() << endl;
    return 1;
  }
  catch (exception& e)
  {
    cout << "Any other exception:" << endl;
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}

