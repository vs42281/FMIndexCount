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
 * Class FMIndexBitString provides functionality of binary string used in
 * conjunction with FM index occurrence wavelet tree.
 *
 * @version     V1.0                26.12.2013.
 * @author      Vedran Sabadoš
 */
public class FMIndexBitString implements java.io.Serializable {

    // <editor-fold desc="Constants">

    final public int BUCKET_SIZE = 64;                              // bucket size
    final public int BUCKET_SIZE_SHIFT = 6;                         // bucktt size multiplication/divsion shift
    final public int BUCKET_REMINDER_BITMASK = (0x0000003f);        // bucket size multiplication/divsion shift
    final public int SUPER_BUCKET_SIZE = 4096;                      // super bucket size
    final public int SUPER_BUCKET_SIZE_SHIFT = 12;                  // super bucket size multiplication/divsion shift
    final public int SUPER_BUCKET_REMINDER_BITMASK = (0x00000fff);  // super bucket size multiplication/divsion shift
    // </editor-fold>

    // <editor-fold desc="Fields">

    /** Bit string storage */
    private final long[] bitString;

    /** Bucket storage */
    private final short[] bucketArray;

    /** Super bucket storage */
    private final int[] superBucketArray;

    /** Number of buckets */
    private final int noBuckets;

    /** Number of  superbuckets */
    private final int noSuperBuckets;
    // </editor-fold>

    // <editor-fold desc="Constructor">

    /**
     *  Class FMIndexBitString Constructor 
     * 
     * @param size      size of bitstring
     */
    public FMIndexBitString(int size) {
        noSuperBuckets = ((size + (SUPER_BUCKET_SIZE - 1)) / SUPER_BUCKET_SIZE);
        noBuckets = (noSuperBuckets * BUCKET_SIZE);
        superBucketArray = new int[noSuperBuckets];
        bucketArray = new short[noBuckets];
        bitString = new long[noBuckets];
    }
    // </editor-fold>

    // <editor-fold desc="Methods">

    /**
     * Function setBit sets bit at specified index.
     * 
     * @param index     index to set bit
     * @throws ReportedException 
     */
    public void setBit(int index) throws ReportedException {

        int bitStringIndex;
        int bucketLimitIndex;
        int superBucketIndex;
        long bitMask;
        long bits;

        bitStringIndex = (index >>> BUCKET_SIZE_SHIFT);
        bitMask = (1L << (index & BUCKET_REMINDER_BITMASK));
        if (bitStringIndex >= noBuckets) {
            reportError ("Character to insert into wavelet tree node is invalid!");
            throw new ReportedException();
        }
        bits = bitString[bitStringIndex];
        if ((bits & bitMask) == 0) {
            bitString[bitStringIndex] = (bits | bitMask);
            superBucketIndex = ((bitStringIndex >>> BUCKET_SIZE_SHIFT) + 1);
            bitStringIndex++;
            bucketLimitIndex = (superBucketIndex << BUCKET_SIZE_SHIFT);
            while (bitStringIndex < bucketLimitIndex) {
                bucketArray[bitStringIndex++]++;
            }
            while (superBucketIndex < noSuperBuckets) {
                superBucketArray[superBucketIndex++]++;
            }
        }
    }

    /**
     * Function setBitNoBucket sets bit at specified index. This function does
     * not affects bucket and superbucket array.
     * 
     * @param index     index to set bit
     * @throws ReportedException 
     */
    public void setBitNoBucket(int index) throws ReportedException {

        int bitStringIndex;
        long bitMask;
        long bits;

        bitStringIndex = (index >>> BUCKET_SIZE_SHIFT);
        bitMask = (1L << (index & BUCKET_REMINDER_BITMASK));
        if (bitStringIndex >= noBuckets) {
            reportError ("Character to insert into wavelet tree node is invalid!");
            throw new ReportedException();
        }
        bits = bitString[bitStringIndex];
        if ((bits & bitMask) == 0) {
            bitString[bitStringIndex] = (bits | bitMask);
        }
    }

    /**
     * Function resetBit resets bit at specified index.
     * 
     * @param index     index to reset bit
     * @throws ReportedException 
     */
    public void resetBit(int index) throws ReportedException {

        int bitStringIndex;
        int bucketLimitIndex;
        int superBucketIndex;
        long bitMask;
        long bits;

        bitStringIndex = (index >>> BUCKET_SIZE_SHIFT);
        bitMask = (1L << (index & BUCKET_REMINDER_BITMASK));
        if (bitStringIndex >= noBuckets) {
            reportError ("Character to insert into wavelet tree node is invalid!");
            throw new ReportedException();
        }
        bits = bitString[bitStringIndex];
        if ((bits & bitMask) != 0) {
            bitString[bitStringIndex] = (bits & (~bitMask));
            superBucketIndex = ((bitStringIndex >>> BUCKET_SIZE_SHIFT) + 1);
            bitStringIndex++;
            bucketLimitIndex = (superBucketIndex << BUCKET_SIZE_SHIFT);
            while (bitStringIndex < bucketLimitIndex) {
                bucketArray[bitStringIndex++]--;
            }
            while (superBucketIndex < noSuperBuckets) {
                superBucketArray[superBucketIndex++]--;
            }
        }
    }

    /**
     * Function resetBitNoBucket resets bit at specified index. This function
     * does not affects bucket and superbucket array.
     * 
     * @param index     index to reset bit
     * @throws ReportedException 
     */
    public void resetBitNoBucket(int index) throws ReportedException {

        int bitStringIndex;
        long bitMask;
        long bits;

        bitStringIndex = (index >>> BUCKET_SIZE_SHIFT);
        bitMask = (1L << (index & BUCKET_REMINDER_BITMASK));
        if (bitStringIndex >= noBuckets) {
            reportError ("Character to insert into wavelet tree node is invalid!");
            throw new ReportedException();
        }
        bits = bitString[bitStringIndex];
        if ((bits & bitMask) != 0) {
            bitString[bitStringIndex] = (bits & (~bitMask));
        }
    }

    /**
     * Function refreshBuckets refreshes bucket and superbucket array to reflect
     * bitstring array condition.
     */
    public void refreshBuckets() {

        int bucketIndex;
        int superSum;
        short sum;

        bucketIndex = -1;
        superSum = 0;
        superBucketArray[0] = 0;
        for (int i = 0; i < noSuperBuckets; i++) {
            sum = 0;
            bucketArray[++bucketIndex] = 0;
            for (int j = 0; j < BUCKET_SIZE; j++) {
                sum += Long.bitCount(bitString[bucketIndex]);
                if (j < (BUCKET_SIZE - 1)) {
                    bucketArray[++bucketIndex] = sum;
                }
            }
            if (i < (noSuperBuckets - 1)) {
                superSum += sum;
                superBucketArray[i + 1] = superSum;
            }
        }
    }

    /**
     * Function rank1 provides number of ones in bitstring before specified index.
     * The bit at specified index is not included in the sum.
     * 
     * @param index     index to count ones to
     * @return          number of ones up to index
     * @throws ReportedException 
     */
    public int rank1(int index) throws ReportedException {

        int inx;
        int reminder;
        long bitmask;
        long bitStrigValue;
        int sum;

        if (index >= (noBuckets << BUCKET_SIZE_SHIFT)) {
            reportError("Pokušaj dohvata broja jedinica do indeksa koji je prevelik!");
            throw new ReportedException();
        }

        inx = (index >>> SUPER_BUCKET_SIZE_SHIFT);
        sum = superBucketArray[inx];
        inx = (index >>> BUCKET_SIZE_SHIFT);
        sum += bucketArray[inx];
        reminder = (index & BUCKET_REMINDER_BITMASK);
        bitmask = ((1L << reminder) - 1);
        bitStrigValue = bitString[inx];
        bitStrigValue &= bitmask;
        sum += Long.bitCount(bitStrigValue);
        return sum;
    }

    /**
     * Function getMemoryConsumption estimates memory consumption of this object.
     * 
     * @return returns      estimated memory consumption in bytes
     */
    public int getMemoryConsumption() {
        return ((noSuperBuckets * 4) + (noBuckets * 10));
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
