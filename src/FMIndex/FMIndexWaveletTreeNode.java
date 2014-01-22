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
 * Class FMIndexWaveletTreeNode provides functionality of wavelet tree node used
 * in conjunction with FM index occurrence wavelet tree.
 *
 * @version     V1.0                26.12.2013.
 * @author      Vedran Sabadoš
 */
public class FMIndexWaveletTreeNode implements java.io.Serializable {

    // <editor-fold desc="Fields">

    /** Wavelet tree node underlaying bit string */
    private FMIndexBitString waveletNodeBitString;

    /** Left child wavelet tree node */
    public FMIndexWaveletTreeNode leftChild;

    /** Right child wavelet tree node */
    public FMIndexWaveletTreeNode rightChild;

    /** Lowest character number stored in this node */
    private int lowCharacterNumber;

    /** Highest character number stored in this node */
    private int highCharacterNumber;

    /** First character number stored in this node as binary 1 */
    private int thresholdCharacterNumber;

    /** Wavelet tree node size in bits */
    private int nodeSize;

    /** Wavelet tree node sequential access index */
    private int sequentialIndex;
    // </editor-fold>

    // <editor-fold desc="Constructors">

    /**
     * Class FMIndexWaveletNode constructor.
     * 
     * @param characterString   character string wavelet tree is to be constructed upon
     * @param prefixSumTable    prefix-sum table for given character string
     * @throws ReportedException 
     */
    public FMIndexWaveletTreeNode(  byte[] characterString,
                                    int[] prefixSumTable    ) throws ReportedException {
        this(characterString, prefixSumTable, 0, (prefixSumTable.length - 2));
    }

    /**
     * Class FMIndexWaveletNode constructor.
     * 
     * @param characterString   character string wavelet tree is to be constructed upon
     * @param prefixSumTable    prefix-sum table for given character string
     * @param lowChar           lowest character for this node
     * @param highChar          highest character for this node
     * @throws ReportedException 
     */
    public FMIndexWaveletTreeNode(  byte[] characterString,
                                    int[] prefixSumTable,
                                    int lowChar,
                                    int highChar            ) throws ReportedException {

        int midpoint;

        lowCharacterNumber = lowChar;
        highCharacterNumber = highChar;

        /* calculate threshold character */
        if ((highCharacterNumber - lowCharacterNumber) < 1) {
            reportError ("Class FMIndexWaveletNode construction error (highChar/lowCharacterNumber difference is less than 1)!");
            throw new ReportedException();
        }
        midpoint = ((prefixSumTable[lowCharacterNumber] + prefixSumTable[highCharacterNumber + 1]) / 2);
        thresholdCharacterNumber = lowCharacterNumber + 1;
        while (prefixSumTable[thresholdCharacterNumber] < midpoint) {
            thresholdCharacterNumber++;
        }
        if (    (prefixSumTable[thresholdCharacterNumber] - midpoint)
                > (midpoint - prefixSumTable[thresholdCharacterNumber - 1]) ) {
            thresholdCharacterNumber--;
        }

        /* create underlaying bitstring */
        nodeSize = (prefixSumTable[highCharacterNumber + 1] - prefixSumTable[lowCharacterNumber]);
        waveletNodeBitString = new FMIndexBitString(nodeSize + 1);

        /* create left child node */
        if ((thresholdCharacterNumber - lowCharacterNumber) >= 2) {
            leftChild = new FMIndexWaveletTreeNode( null,
                                                    prefixSumTable,
                                                    lowCharacterNumber,
                                                    (thresholdCharacterNumber - 1)  );
        } else {
            leftChild = null;
        }

        /* create right child node */
        if ((highCharacterNumber - thresholdCharacterNumber) >= 1) {
            rightChild = new FMIndexWaveletTreeNode(    null,
                                                        prefixSumTable,
                                                        thresholdCharacterNumber,
                                                        highCharacterNumber         );
        } else {
            rightChild = null;
        }

        /* populate wavelet tree */
        if (characterString != null) {
            resetSequentialIndex();
            for (int i = 0; i < nodeSize; i++) {
                insertCharacter(characterString[i]);
            }
            refreshBuckets();
        }
    }
    // </editor-fold>

    // <editor-fold desc="Methods">

    /**
     * Function insertCharacter inserts a character to this wavelet tree node
     * and all the nodes below, at the place defined by sequential index.
     * Sequential index is incremented after this operation.
     * 
     * @param characterToInsert     character to insert
     * @throws ReportedException 
     */
    final public void insertCharacter(byte characterToInsert) throws ReportedException {

        /* check if character and sequential index are valid */
        if (    (characterToInsert < lowCharacterNumber)
                || (characterToInsert > highCharacterNumber)    ) {
            reportError ("Character to insert into wavelet tree node is invalid!");
            throw new ReportedException();
        }
        if (sequentialIndex >= nodeSize) {
            reportError ("Traying to insert too many characters into wavelet tree node!");
            throw new ReportedException();
        }

        if (characterToInsert < thresholdCharacterNumber) {
            waveletNodeBitString.resetBitNoBucket(sequentialIndex++);
            if (leftChild != null) {
                leftChild.insertCharacter(characterToInsert);
            }
        } else {
            waveletNodeBitString.setBitNoBucket(sequentialIndex++);
            if (rightChild != null) {
                rightChild.insertCharacter(characterToInsert);
            }
        }
    }

    /**
     * Function resetSequentialIndex resets sequential access index for this node
     * and all the nodes below.
     */
    final public void resetSequentialIndex() {
        sequentialIndex = 0;
        if (leftChild != null) {
            leftChild.resetSequentialIndex();
        }
        if (rightChild != null) {
            rightChild.resetSequentialIndex();
        }
    }

    /**
     * Function refreshBuckets refreshes bucket and superbucket arrays in 
     * underlaying bit strings for this node and all the nodes below.
     */
    final public void refreshBuckets() {
        waveletNodeBitString.refreshBuckets();
        if (leftChild != null) {
            leftChild.refreshBuckets();
        }
        if (rightChild != null) {
            rightChild.refreshBuckets();
        }
    }

    /**
     * Function Occ counts number of occurrences of given character in BW
     * transform of reference sequence before given index. Character at given
     * index position is not included in the sum.
     * 
     * @param countedCharacter  character to count its occurrence
     * @param index             index to count characters to
     * @return                  number of occurrence of countedCharacter up to index position
     * @throws ReportedException 
     */
    public int Occ(byte countedCharacter, int index) throws ReportedException {

        int childIndex;

        if ((countedCharacter < lowCharacterNumber) || (countedCharacter > highCharacterNumber)) {
            reportError("Pokušaj dohvata broja znakova za znak koji ne postoji!");
            throw new ReportedException();
        }
        if (index > nodeSize) {
            reportError("Pokušaj dohvata broja znakova za prevelik indeks!");
            throw new ReportedException();
        }

        if (countedCharacter < thresholdCharacterNumber) {
            childIndex = (index - waveletNodeBitString.rank1(index));
            if (leftChild != null) {
                return leftChild.Occ(countedCharacter, childIndex);
            } else {
                return childIndex;
            }
        } else {
            childIndex = waveletNodeBitString.rank1(index);
            if (rightChild != null) {
                return rightChild.Occ(countedCharacter, childIndex);
            } else {
                return childIndex;
            }
        }
    }

    /**
     * Function getMemoryConsumption estimates memory consumption of this node
     * and all the nodes below.
     * 
     * @return returns      estimated memory consumption in bytes
     */
    final public int getMemoryConsumption() {

        int sum;

        sum = waveletNodeBitString.getMemoryConsumption();
        if (leftChild != null) {
            sum += leftChild.getMemoryConsumption();
        }
        if (rightChild != null) {
            sum += rightChild.getMemoryConsumption();
        }
        return sum;
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
