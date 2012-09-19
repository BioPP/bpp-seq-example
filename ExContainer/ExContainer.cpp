/*
 * File: ExContainer.cpp
 * Created by: Julien Dutheil
 * Created on: Dec Fri 05 11:54 2008
 * Last modified: Jun Fri 05 11:49 2009
 * 
 * Introduction to the Container classes, and reader/writing from sequence files.
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
 * From SeqLib:
 */
#include <Bpp/Seq/Alphabet.all> /* this includes all alphabets in one shot */
#include <Bpp/Seq/Container.all> /* this includes all containers */
#include <Bpp/Seq/Io.all> /* this includes all sequence readers and writers */
#include <Bpp/Seq/StateProperties/DefaultNucleotideScore.h>

/*
 * We'll need a few tools from the Bio++ core library:
 */
#include <Bpp/App/ApplicationTools.h>

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
     * Containers are objects that store a set of sequences.
     * They hence play a major role in the Bio++ libraries.
     *
     * We will first create a simple container from scratch:
     */
    Sequence* seq1 = new BasicSequence("Sequence 1", "GATTACTGATTACATGGT", &AlphabetTools::DNA_ALPHABET);
    Sequence* seq2 = new BasicSequence("Sequence 2", "GATCACAATGTTACGCT", &AlphabetTools::DNA_ALPHABET);
    VectorSequenceContainer* vsc = new VectorSequenceContainer(&AlphabetTools::DNA_ALPHABET);
    vsc->addSequence(*seq1);
    vsc->addSequence(*seq2);
    cout << "This container has " << vsc->getNumberOfSequences() << " sequences." << endl;

    /*
     * The VectorSequenceContainer class is one of the simplest container.
     * It is a mere vector of sequences (it is hence an OrderedSequenceContainer).
     *
     * Sequence can be accessed either by their name or by their position:
     */
    cout << vsc->getSequence("Sequence 1").toString() << endl;
    cout << vsc->getSequence(0).toString() << endl;
    /*
     * Another important point is that sequences are cloned while being added to the container,
     * so that the container has full control over the data. As a consequence:
     * - deleting the sequence object is safe and does not damage the container
     * - deleting the container will not damage your initial sequences.
     * As a consequence, a sequence object can be safely added to several containers and being destroyed after that:
     */
    delete seq1;
    delete seq2;
    cout << vsc->getSequence("Sequence 1").toString() << endl;
    /* still works! */

    /*
     * As before, we have several utilitary methods available.
     */
    cout << "Is that an alignment? " << (SequenceContainerTools::sequencesHaveTheSameLength(*vsc) ? "yes" : "no") << endl;
    SiteContainer* alignedSeq = SiteContainerTools::alignNW(vsc->getSequence(0), vsc->getSequence(1), DefaultNucleotideScore(&AlphabetTools::DNA_ALPHABET), -5);
    cout << alignedSeq->getSequence(0).toString() << endl;
    cout << alignedSeq->getSequence(1).toString() << endl;
    cout << "Is that an alignment? " << (SequenceContainerTools::sequencesHaveTheSameLength(*alignedSeq) ? "yes" : "no") << endl;

    /*
     * Alignement are stored in a particular structure called SiteContainer.
     * A SiteContainer has all the properties of a SequenceContainer + specific methods
     * to access the data per column. Check the documentation of those classes.
     */
    cout << "This container has " << alignedSeq->getNumberOfSites() << " positions." << endl;
    delete alignedSeq;

    /* ----------------
     * QUESTION 1: Look at the SiteContainerTools class, and try to:
     * - remove the gaps from the alignment
     * - change the gaps to unresolve characters
     * - compute the frequencies of nucleotides
     * ----------------
     */

    /*
     * In most cases, containers are created from a file and not from scratch.
     * We hence need to introduce a new kind of objects, called sequence reader (ISequence).
     * ISequence objects have a method 'read' that creates a container from a file and an alphabet.
     */
    Fasta fasReader;
    OrderedSequenceContainer* sequences = fasReader.readSequences("TIMnuc.aln.fasta", &AlphabetTools::DNA_ALPHABET);
    cout << "This container has " << sequences->getNumberOfSequences() << " sequences." << endl;
    cout << "Is that an alignment? " << (SequenceContainerTools::sequencesHaveTheSameLength(*sequences) ? "yes" : "no") << endl;

    /*
     * The Fasta format can store sequences which are aligned or not.
     * It hence returns a SequenceContainer object, not an alignment.
     * Since we know that the sequences are aligned, we will create a SiteContainer
     * object from the SequenceContainer. There are two types of SiteContainer objects in Bio++,
     * (i) the AlignedSequenceContainer, which is a specialization of the VectorSequenceContainer,
     * and (ii) the VectorSiteContainer. They differ by the way they store the data.
     * We will compare them to assess their properties:
     */
    SiteContainer* align = new AlignedSequenceContainer(*sequences);
    SiteContainer* sites = new VectorSiteContainer(*sequences);

    unsigned int nrep = 1000; /* WATCHOUT!!! reduce this number if this is too slow on your computer! */

    ApplicationTools::startTimer();
    for (unsigned int j = 0; j < nrep; j++)
    {
      ApplicationTools::displayGauge(j, nrep-1, '=');
      for (unsigned int i = 0; i < align->getNumberOfSequences(); i++)
        align->getSequence(i);
    }
    ApplicationTools::displayTime("\nTotal time used for sequence access in ASC:");

    ApplicationTools::startTimer();
    for (unsigned int j = 0; j < nrep; j++)
    {
      ApplicationTools::displayGauge(j, nrep-1, '=');
      for (unsigned int i = 0; i < sites->getNumberOfSequences(); i++)
        sites->getSequence(i);
    }
    ApplicationTools::displayTime("\nTotal time used for sequence access in VSC:");

    ApplicationTools::startTimer();
    for (unsigned int j = 0; j < nrep; j++)
    {
      ApplicationTools::displayGauge(j, nrep-1, '=');
      for (unsigned int i = 0; i < align->getNumberOfSites(); i++)
        align->getSite(i);
    }
    ApplicationTools::displayTime("\nTotal time used for site access in ASC:");

    ApplicationTools::startTimer();
    for (unsigned int j = 0; j < nrep; j++)
    {
      ApplicationTools::displayGauge(j, nrep-1, '=');
      for (unsigned int i = 0; i < sites->getNumberOfSites(); i++)
        sites->getSite(i);
    }
    ApplicationTools::displayTime("\nTotal time used for site access in VSC:");

    delete sequences;
    delete align;
    delete sites;

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

