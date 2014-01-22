/*
 * Author:
 *
 *      Vedran Sabadoš
 *
 * Version:
 *
 *      V1.0         26. 12. 2013
 */

package FMIndex;

import javax.swing.JOptionPane;

/**
 * Class ReferenceSequenceIndex contains all data structures constructed
 * upon reference sequence necessary to efficiently count occurrences of
 * query sequence.
 * 
 * @version     V1.0                26.12.2013.
 * @author      Vedran Sabadoš
 */
public class ReferenceSequenceIndex implements java.io.Serializable {

    // <editor-fold desc="Fields">

    /** Actual alphabet dictionary contains one byte for each possible character */
    public AlphabetDictionary alphabetDictionary;

    /** Prefix-sum array */
    public int[] prefixSumTable;

    /** Occurrence Wavelet tree */
    public FMIndexWaveletTreeNode occurrenceWaveletTree;
    // </editor-fold>

    // <editor-fold desc="Methods">

    /**
     * Function C provides total number of characters which are lexicographically
     * smaller than given character. If given character is for 1 greater than the
     * greatest character, than total number of characters in BWT is returned.
     * 
     * @param characterToCount      character to which to count
     * @return                      total number of lexicographically smaller characters
     * @throws ReportedException 
     */
    public int C(byte characterToCount) throws ReportedException {
        if (characterToCount > (alphabetDictionary.getAlphabetSize()) + 1) {
            reportError("Pokušaj dohvata broja nižih znakova za znak veći od dopuštenog!");
            throw new ReportedException();
        }
        return prefixSumTable[characterToCount];
    }

    /**
     * Function getCharacterTotal provides how many times specified character is
     * contained into reference sequence.
     * 
     * @param characterToGet    character whose occurrence is requested
     * @return                  number how many times characterToGet is contained in reference sequence
     * @throws ReportedException 
     */
    public int getCharacterTotal(byte characterToGet) throws ReportedException {
        if (characterToGet > alphabetDictionary.getAlphabetSize()) {
            reportError("Pokušaj dohvata broja znakova za znak koji ne postoji!");
            throw new ReportedException();
        }
        return prefixSumTable[characterToGet + 1] - prefixSumTable[characterToGet];
    }

    /**
     * Function reportError reports execution errors by opening confirmation
     * window with appropriate message.
     * 
     * @param   text    message to be reported in confirmation window
     */
    private void reportError (String text) {
        JOptionPane.showConfirmDialog(null, text, "FMIndexCount", JOptionPane.DEFAULT_OPTION, JOptionPane.ERROR_MESSAGE);
    }
    // </editor-fold>
}
