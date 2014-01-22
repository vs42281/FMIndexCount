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
 * Class AlphabetDictionary provides dictionary for characters used in reference
 * string and methods necessary to convert it to compact binary encoding.
 *
 * @version     V1.0                26.12.2013.
 * @author      Vedran Sabadoš
 */
public class AlphabetDictionary implements java.io.Serializable {
    /*
     * this version supports only characters whose binary value is in range 1 - 127.
     */

    // <editor-fold desc="Fields">

    /** Actual alphabet dictionary contains one byte for each possible character */
    private final byte[] dictionaryArray;

    /** Number of characters in the alphabet */
    private int alphabetSize;

    /** Flag telling whether dictionary is finalized or not */
    private boolean finalized;
    // </editor-fold>

    // <editor-fold desc="Constructor">

    /** Class AlphabetDictionary Constructor */
    public AlphabetDictionary() {
        dictionaryArray = new byte[128];
        finalized = false;
    }
    // </editor-fold>

    // <editor-fold desc="Methods">

    /**
     * Function addCharacter adds new character into not finalized alphabet
     * dictionary.
     * 
     * @param newCharacter      A character to be added to the dictionary
     * @throws ReportedException 
     */
    public void addCharacter(int newCharacter) throws ReportedException {
        if (finalized) {
            reportError ("Pokušaj dodavanja znaka već finaliziranom rječniku!");
            throw new ReportedException();
        }
        if ((newCharacter <= 0) || (newCharacter > 127)) {
            reportError ("Pokušaj dodavanja u rječnik znaka izvan dozvoljenih granica: " + Integer.toString(newCharacter));
            throw new ReportedException();
        }
        dictionaryArray[newCharacter] = 1;
    }

    /**
     * Function getCompactCode provides compact code for raw coded character
     * 
     * @param rawCode       raw code of character the compact code is to be provided for
     * @return              compact code for given character
     * @throws ReportedException 
     */
    public byte getCompactCode(byte rawCode) throws ReportedException {

        if (!finalized) {
            reportError ("Pokušaj dohvata kompaktnog koda iz nefinaliziranog rječnika!");
            throw new ReportedException();
        }
        return dictionaryArray[rawCode];
    }

    /**
     * Function getAlphabetSize provides number of characters in alphabet
     * 
     * @return      number of characters in alphabet
     * @throws ReportedException 
     */
    public int getAlphabetSize() throws ReportedException {

        if (!finalized) {
            reportError ("Pokušaj dohvata veličine abecede iz nefinaliziranog rječnika!");
            throw new ReportedException();
        }
        return alphabetSize;
    }

    /**
     * Function finalizeDictionary finalizes alphabet dictionary, that is,
     * provides compact code for each character from alphabet.
     * 
     * @throws ReportedException 
     */
    public void finalizeDictionary() throws ReportedException {

        byte compactCode = 1;

        if (finalized) {
            reportError ("Pokušaj finaliziranja već finaliziranom rječniku!");
            throw new ReportedException();
        }
        for (int i = 1; i < 128; i++) {
            if (dictionaryArray[i] != 0) {
                dictionaryArray[i] = compactCode;
                compactCode++;
            }
        }
        alphabetSize = (compactCode - 1);
        finalized = true;
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
