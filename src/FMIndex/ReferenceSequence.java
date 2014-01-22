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

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import static java.lang.StrictMath.min;
import static java.lang.StrictMath.pow;
import static java.lang.StrictMath.round;
import java.nio.charset.Charset;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.LinkedList;
import javax.swing.JOptionPane;

/**
 * Class ReferenceSequence represents reference sequence to which all query
 * sequences are compared to find number of occurrences.
 * 
 * @version     V1.0                26.12.2013.
 * @author      Vedran Sabadoš
 */
public class ReferenceSequence {

    // <editor-fold desc="Constants">

    final public int NEW_LINE_CHARACTER = 0x0a;                 // \n character in UTF8
    final public int SPACE_CHARACTER = 0x20;                    // space character in UTF8
    final public int SEQUENCE_NAME_START_CHARACTER = 0x3e;      // > character in UTF8
    final public int LOWRCASE_A_CHARACTER = 0x61;               // a character in UTF8
    final public int LOWRCASE_Z_CHARACTER = 0x7a;               // z character in UTF8

    final public int INITIAL_START_SEQUENCE_LENGTH_5 = 3;       // initial start sequence for alphabet of no more than 5 characters
    final public int INITIAL_START_SEQUENCE_LENGTH_10 = 3;      // initial start sequence for alphabet of no more than 10 characters
    final public int INITIAL_START_SEQUENCE_LENGTH_10_PLUS = 2; // initial start sequence for alphabet of more than 10 characters

    final public int READ_BUFFER_SIZE = 1000;                   // size of read buffer
    final public int LINKED_BUFFER_SIZE = 1000;                 // size of linked buffer
    final public int START_SEQUENCE_BUFFER_SIZE = 1000;         // size of start sequence buffer
    final public int ROTATION_BUFFER_SIZE = 10000000;            // size of rotation buffer
    final public int LONG_RUN_SEQUENCE_THRESHOLD = 100;         // size of rotation buffer
    final public int TEST_LONG_RUN_SEQUENCE_AFTER = 10;         // test if long run sequence after specified number of passes when comaring two rotations
    // </editor-fold>

    // <editor-fold desc="Fields">

    /** Object containing all data structures necessary to efficiently count
     * occurrences of query sequence.
     */
    private ReferenceSequenceIndex referenceSequenceIndex;

    /** Parent frame reference */
    FMIndexCountFrame parentFrame;

    /** Compact coded reference sequence buffer */
    private byte[] compactCodeBuffer;

    /** Reference sequence Burrows-Wheeler transform buffer */
    private byte[] bwtBuffer;

    /** Reference sequence Burrows-Wheeler transform buffer index */
    int bwtBufferIndex;

    /** Reference sequence file read stream */
    FileInputStream referenceSequenceStream;

    /** Reference sequence file read intermediate buffer */
    private byte[] intermediateReadBuffer;

    /** Reference sequence file read intermediate buffer index */
    private int intermediateReadBufferIndex;

    /** Number of still not read characters in reference sequence file */
    long referenceSequenceNotRead;

    /** Current start sequence buffer */
    private byte[] startSequenceBuffer;

    /** Array with count of long run sequences per character */
    private int[] longRunCount;

    /** Array with start rotations of long run sequences per character */
    private int[][] longRunStart;

    /** Array with lengths of long run sequences per character */
    private int[][] longRunLength;

    /** Objects for reference sequence statistics */
    private SimpleDateFormat sdf;
    private String startTimeBuffer;
    private String fileReadTimeBuffer;
    private String bwtCreationTimeBuffer;
    private String indexCreationTimeBuffer;
    private int noStartSequenceCombinations;
    private int currentStartSequenceCombination;
    private long maximalPreprocessingMemoryConsumption;
    private long currentPreprocessingMemoryConsumption;

//    private int[] testBwtBuffer;
    // </editor-fold>

    // <editor-fold desc="Constructor">

    /**
     * Class ReferenceSequence Constructor
     * 
     * @param filename      Reference sequence filename
     * @param format        Reference sequence file format
     * @param sequenceName  Reference sequence name (in case of Fasta format)
     * @param parent        Reference for parent frame
     * @throws ReportedException 
     */
    @SuppressWarnings({"BroadCatchBlock", "TooBroadCatch"})
    public ReferenceSequence(   String filename,
                                ReferenceSequenceFormat format,
                                String sequenceName,
                                FMIndexCountFrame parent        ) throws ReportedException {

        parentFrame = parent;
        sdf = new SimpleDateFormat("HH:mm:ss.SSS");
        startTimeBuffer = sdf.format(Calendar.getInstance().getTime());
        System.out.println("Start time: " + startTimeBuffer);

        /* initialize memory occupation fields */
        maximalPreprocessingMemoryConsumption = 0;
        currentPreprocessingMemoryConsumption = 0;

        if (format == ReferenceSequenceFormat.FASTA) {
            try {
                preprocessFastaReferenceSequence(filename, sequenceName);
            } catch (ReportedException ex) {
                throw ex;
            } catch (Exception ex) {
                reportError(    "Neočekivana pogreška tijekom predprocesiranja referentnog slijeda u fasta formatu ("
                                + getExceptionType(ex)
                                + "):\n\n    "
                                + ex.getMessage()                                                                       );
                throw new ReportedException();
            }
        } else if (format == ReferenceSequenceFormat.TEXT) {
            try {
                preprocessTextReferenceSequence(filename);
            } catch (ReportedException ex) {
                throw ex;
            } catch (Exception ex) {
                reportError(    "Neočekivana pogreška tijekom predprocesiranja referentnog slijeda u tekst formatu ("
                                + getExceptionType(ex)
                                + "):\n\n    "
                                + ex.getMessage()                                                                       );
                throw new ReportedException();
            }
        } else {
            try {
                readPreprocessedReferenceSequence(filename);
            } catch (ReportedException ex) {
                throw ex;
            } catch (Exception ex) {
                reportError(    "Neočekivana pogreška tijekom učitavanja predprocesiranog referentnog slijeda ("
                                + getExceptionType(ex)
                                + "):\n\n    "
                                + ex.getMessage()                                                                       );
                throw new ReportedException();
            }
        }

        /* report statistics */
        reportStatistics();

        /* update maximal memory consumption field */
        maximalPreprocessingMemoryConsumption = currentPreprocessingMemoryConsumption;
    }
    // </editor-fold>

    // <editor-fold desc="Methods">

    /**
     * Function preprocessFastaReferenceSequence preprocesses reference sequence
     * from file in Fasta format.
     * 
     * @param filename      Reference sequence filename
     * @param sequenceName  Reference sequence name (in case of Fasta format)
     * @throws ReportedException 
     */
    private void preprocessFastaReferenceSequence(  String filename,
                                                    String sequenceName ) throws ReportedException {

        File referenceSequenceFile;
        Charset utf8Charset;
        byte[] sequenceNameBuffer;
        LinkedList<byte[]> linkedBufferList;
        byte[] linkedBuffer;
        long sequenceLength;
        int linkedBufferIndex;
        int nextChar;
        int compactCodeBufferIndex;

        /* Create reference sequence index */
        referenceSequenceIndex = new ReferenceSequenceIndex();
        referenceSequenceIndex.alphabetDictionary = new AlphabetDictionary();
        currentPreprocessingMemoryConsumption += 128;
        if (currentPreprocessingMemoryConsumption > maximalPreprocessingMemoryConsumption) {
            maximalPreprocessingMemoryConsumption = currentPreprocessingMemoryConsumption;
        }

        /* Set progress bar to indeterminate */
        java.awt.EventQueue.invokeLater(new Runnable() {
            @Override
            public void run() {
                parentFrame.setProgressBarIndeterminate(true);
            }
        });

        /* prepare sequence name buffer */
        utf8Charset = Charset.forName("UTF-8");
        sequenceNameBuffer = (">" + sequenceName).getBytes(utf8Charset);
        utf8Charset = null;

        /* open reference sequence file */
        try {
            referenceSequenceFile = new File(filename);
            referenceSequenceNotRead = referenceSequenceFile.length();
            referenceSequenceStream = null;
            referenceSequenceStream = new FileInputStream(referenceSequenceFile);
            intermediateReadBuffer = new byte[READ_BUFFER_SIZE];
            intermediateReadBufferIndex = READ_BUFFER_SIZE;
            currentPreprocessingMemoryConsumption += READ_BUFFER_SIZE;
            if (currentPreprocessingMemoryConsumption > maximalPreprocessingMemoryConsumption) {
                maximalPreprocessingMemoryConsumption = currentPreprocessingMemoryConsumption;
            }

            /* find start of given sequence */
            while (!isSequenceName(sequenceNameBuffer)) {
                checkCancel();
                if (isEndOfFile()) {
                    break;
                }
            }
            if (isEndOfFile()) {
                reportError(    "Sekvenca \""
                                + sequenceName
                                + "\" nije sadržana u datoteci:\n\n"
                                + filename                          );
                throw new ReportedException();
            }

            /* read given sequence */
            linkedBufferList = new LinkedList();
            linkedBuffer = new byte[LINKED_BUFFER_SIZE];
            sequenceLength = 0;
            linkedBufferIndex = 0;
            while ((nextChar = getCharacter()) >= 0) {
                if (nextChar == NEW_LINE_CHARACTER) {
                    continue;
                }
                if (nextChar == SEQUENCE_NAME_START_CHARACTER) {
                    break;
                }
                nextChar = toUpperCase(nextChar);
                referenceSequenceIndex.alphabetDictionary.addCharacter(nextChar);
                if (linkedBufferIndex >= LINKED_BUFFER_SIZE) {
                    checkCancel();
                    linkedBufferList.add(linkedBuffer);
                    sequenceLength += LINKED_BUFFER_SIZE;
                    linkedBufferIndex = 0;
                    linkedBuffer = new byte[LINKED_BUFFER_SIZE];
                }
                linkedBuffer[linkedBufferIndex++] = ((byte) nextChar);
            }
        } catch (ReportedException ex) {
            throw ex;
        } catch (Exception ex) {
            reportError(    "Greška tijekom otvaranja datoteke referentnog slijeda ("
                            + getExceptionType(ex)
                            + "):\n\n    "
                            + ex.getMessage()                                           );
            throw new ReportedException();
        } finally {
            if (referenceSequenceStream != null) {
                try {
                    referenceSequenceStream.close();
                } catch (IOException ex) {
                }
            }
        }
        sequenceLength += linkedBufferIndex;
        if ((sequenceLength >= 0x7fffffff) || (sequenceLength < 2)) {     //leave place for final character
            reportError(    "Referentni slijed je prevelik: "
                            + Long.toString(sequenceLength)
                            + " znakova!"                       );
            throw new ReportedException();
        }

        System.out.println("Reference sequence readout completed at: " + sdf.format(Calendar.getInstance().getTime()));


        /* create compact code reference sequence buffer */
        currentPreprocessingMemoryConsumption -= READ_BUFFER_SIZE;
        referenceSequenceStream = null;
        intermediateReadBuffer = null;
        referenceSequenceIndex.alphabetDictionary.finalizeDictionary();
        compactCodeBuffer = new byte[(int) (sequenceLength + 1)];
        currentPreprocessingMemoryConsumption += (sequenceLength + 1);
        currentPreprocessingMemoryConsumption += (linkedBufferList.size() * 8);
        currentPreprocessingMemoryConsumption += (   (linkedBufferList.size() + 1)
                                                    * LINKED_BUFFER_SIZE            );
        if (currentPreprocessingMemoryConsumption > maximalPreprocessingMemoryConsumption) {
            maximalPreprocessingMemoryConsumption = currentPreprocessingMemoryConsumption;
        }

        /* fill compact code reference sequence buffer */
        compactCodeBufferIndex = 0;
        for (byte[] buf: linkedBufferList) {
            for (byte i: buf) {
                compactCodeBuffer[compactCodeBufferIndex++] =
                        referenceSequenceIndex.alphabetDictionary.getCompactCode(i);
            }
        }
        for (int i = 0; i < linkedBufferIndex; i++) {
            compactCodeBuffer[compactCodeBufferIndex++] =
                    referenceSequenceIndex.alphabetDictionary.getCompactCode(linkedBuffer[i]);
        }
        compactCodeBuffer[compactCodeBufferIndex++] = 0;

        /* release unused memory */
        currentPreprocessingMemoryConsumption -= (linkedBufferList.size() * 8);
        currentPreprocessingMemoryConsumption -= (   (linkedBufferList.size() + 1)
                                                    * LINKED_BUFFER_SIZE            );
        linkedBufferList = null;
        linkedBuffer = null;
        checkCancel();

        fileReadTimeBuffer = sdf.format(Calendar.getInstance().getTime());
        System.out.println("Reference sequence compact code buffer creation completed at: " + fileReadTimeBuffer);

        preprocessCompactReferenceSequence();
    }


    /**
     * Function preprocessTextReferenceSequence preprocesses reference sequence
     * from file in Text format.
     * 
     * @param filename      Reference sequence filename
     * @throws ReportedException 
     */
    private void preprocessTextReferenceSequence(String filename) throws ReportedException {

        File referenceSequenceFile;
        LinkedList<byte[]> linkedBufferList;
        byte[] linkedBuffer;
        long sequenceLength;
        int linkedBufferIndex;
        int nextChar;
        int compactCodeBufferIndex;

        /* Create reference sequence index */
        referenceSequenceIndex = new ReferenceSequenceIndex();
        referenceSequenceIndex.alphabetDictionary = new AlphabetDictionary();
        currentPreprocessingMemoryConsumption += 128;
        if (currentPreprocessingMemoryConsumption > maximalPreprocessingMemoryConsumption) {
            maximalPreprocessingMemoryConsumption = currentPreprocessingMemoryConsumption;
        }

        /* Set progress bar to indeterminate */
        java.awt.EventQueue.invokeLater(new Runnable() {
            @Override
            public void run() {
                parentFrame.setProgressBarIndeterminate(true);
            }
        });

        /* open reference sequence file */
        try {
            referenceSequenceFile = new File(filename);
            referenceSequenceNotRead = referenceSequenceFile.length();
            referenceSequenceStream = null;
            referenceSequenceStream = new FileInputStream(referenceSequenceFile);
            intermediateReadBuffer = new byte[READ_BUFFER_SIZE];
            intermediateReadBufferIndex = READ_BUFFER_SIZE;
            currentPreprocessingMemoryConsumption += READ_BUFFER_SIZE;
            if (currentPreprocessingMemoryConsumption > maximalPreprocessingMemoryConsumption) {
                maximalPreprocessingMemoryConsumption = currentPreprocessingMemoryConsumption;
            }
            if (isEndOfFile()) {
                reportError(    "Datoteka \""
                                + filename
                                + "\" ne sadrži niti jedan znak\n\n");
                throw new ReportedException();
            }

            /* read given sequence */
            linkedBufferList = new LinkedList();
            linkedBuffer = new byte[LINKED_BUFFER_SIZE];
            sequenceLength = 0;
            linkedBufferIndex = 0;
            while ((nextChar = getCharacter()) >= 0) {
                if ((nextChar < 32) || (nextChar > 126)) {
                    continue;
                }
                nextChar = toUpperCase(nextChar);
                referenceSequenceIndex.alphabetDictionary.addCharacter(nextChar);
                if (linkedBufferIndex >= LINKED_BUFFER_SIZE) {
                    checkCancel();
                    linkedBufferList.add(linkedBuffer);
                    sequenceLength += LINKED_BUFFER_SIZE;
                    linkedBufferIndex = 0;
                    linkedBuffer = new byte[LINKED_BUFFER_SIZE];
                }
                linkedBuffer[linkedBufferIndex++] = ((byte) nextChar);
            }
        } catch (ReportedException ex) {
            throw ex;
        } catch (Exception ex) {
            reportError(    "Greška tijekom otvaranja datoteke referentnog slijeda ("
                            + getExceptionType(ex)
                            + "):\n\n    "
                            + ex.getMessage()                                           );
            throw new ReportedException();
        } finally {
            if (referenceSequenceStream != null) {
                try {
                    referenceSequenceStream.close();
                } catch (IOException ex) {
                }
            }
        }
        sequenceLength += linkedBufferIndex;
        if ((sequenceLength >= 0x7fffffff) || (sequenceLength < 2)) {     //leave place for final character
            reportError(    "Referentni slijed je prevelik: "
                            + Long.toString(sequenceLength)
                            + " znakova!"                       );
            throw new ReportedException();
        }

        System.out.println("Reference sequence readout completed at: " + sdf.format(Calendar.getInstance().getTime()));


        /* create compact code reference sequence buffer */
        currentPreprocessingMemoryConsumption -= READ_BUFFER_SIZE;
        referenceSequenceStream = null;
        intermediateReadBuffer = null;
        referenceSequenceIndex.alphabetDictionary.finalizeDictionary();
        compactCodeBuffer = new byte[(int) (sequenceLength + 1)];
        currentPreprocessingMemoryConsumption += (sequenceLength + 1);
        currentPreprocessingMemoryConsumption += (linkedBufferList.size() * 8);
        currentPreprocessingMemoryConsumption += (   (linkedBufferList.size() + 1)
                                                    * LINKED_BUFFER_SIZE            );
        if (currentPreprocessingMemoryConsumption > maximalPreprocessingMemoryConsumption) {
            maximalPreprocessingMemoryConsumption = currentPreprocessingMemoryConsumption;
        }

        /* fill compact code reference sequence buffer */
        compactCodeBufferIndex = 0;
        for (byte[] buf: linkedBufferList) {
            for (byte i: buf) {
                compactCodeBuffer[compactCodeBufferIndex++] =
                        referenceSequenceIndex.alphabetDictionary.getCompactCode(i);
            }
        }
        for (int i = 0; i < linkedBufferIndex; i++) {
            compactCodeBuffer[compactCodeBufferIndex++] =
                    referenceSequenceIndex.alphabetDictionary.getCompactCode(linkedBuffer[i]);
        }
        compactCodeBuffer[compactCodeBufferIndex++] = 0;

        /* release unused memory */
        currentPreprocessingMemoryConsumption -= (linkedBufferList.size() * 8);
        currentPreprocessingMemoryConsumption -= (   (linkedBufferList.size() + 1)
                                                    * LINKED_BUFFER_SIZE            );
        linkedBufferList = null;
        linkedBuffer = null;
        checkCancel();

        System.out.println("Reference sequence compact code buffer creation completed at: " + sdf.format(Calendar.getInstance().getTime()));

        preprocessCompactReferenceSequence();
    }

    /**
     * Function readPreprocessedReferenceSequence reads preprocessed and saved
     * reference sequence index from file.
     * 
     * @param filename      Filename of file containing saved reference sequence index
     * @throws ReportedException 
     */
    private void readPreprocessedReferenceSequence(String filename) throws ReportedException {

        ObjectInputStream in = null;

        try {
            in = new ObjectInputStream(new FileInputStream(filename));
            try {
                referenceSequenceIndex = new ReferenceSequenceIndex();
                referenceSequenceIndex = (ReferenceSequenceIndex) in.readObject();

                /* estimate occurence wavelet node memory consumption */
                currentPreprocessingMemoryConsumption = 128;
                currentPreprocessingMemoryConsumption += (referenceSequenceIndex.prefixSumTable.length * 4);
                currentPreprocessingMemoryConsumption += referenceSequenceIndex.occurrenceWaveletTree.getMemoryConsumption();
                maximalPreprocessingMemoryConsumption = currentPreprocessingMemoryConsumption;
            } catch (Exception ex) {
                reportError(    "Greška tijekom učitavanja predprocesirane datoteke referentnog slijeda ("
                                + getExceptionType(ex)
                                + "):\n\n    "
                                + ex.getMessage()                                           );
                throw new ReportedException();
            }
        } catch (ReportedException ex) {
            throw ex;
        } catch (Exception ex) {
            reportError(    "Greška tijekom otvaranja datoteke za učitavanje predprocesiranog referentnog slijeda ("
                            + getExceptionType(ex)
                            + "):\n\n    "
                            + ex.getMessage()                                           );
            throw new ReportedException();
        } finally {
            if (in != null) {
                try {
                    in.close();
                } catch (IOException ex) {
                }
            }
        }

        indexCreationTimeBuffer = sdf.format(Calendar.getInstance().getTime());
        System.out.println("Reference sequence index creation completed at: " + indexCreationTimeBuffer);
    }

    /**
     * Function preprocessCompactReferenceSequence preprocesses reference sequence
     * coded with compact code and stored in compactCodeBuffer. Alphabet is
     * defined by alphabetDictionary object.
     */
    private void preprocessCompactReferenceSequence() throws ReportedException {
        createBwt();
        createIndex();

        /* estimate memory consumption */
        currentPreprocessingMemoryConsumption = 128;
        currentPreprocessingMemoryConsumption += (referenceSequenceIndex.prefixSumTable.length * 4);
        currentPreprocessingMemoryConsumption += referenceSequenceIndex.occurrenceWaveletTree.getMemoryConsumption();
    }

    /**
     * Function createBwt creates reference sequence BW transform.
     * 
     * @throws ReportedException 
     */
    private void createBwt() throws ReportedException {

        int[][] rotationBuffers;
        int compactCodeBufferSize;
        int alphabetSize;
        int initial_start_sequence_length;
        int changeDepth;
        String startSequenceString;



//        testBwtBuffer = new int[compactCodeBuffer.length];


        compactCodeBufferSize = compactCodeBuffer.length;
        bwtBuffer = new byte[compactCodeBufferSize];
        bwtBuffer[0] = compactCodeBuffer[compactCodeBufferSize - 2];
        bwtBufferIndex = 1;



//        testBwtBuffer[0] = compactCodeBufferSize - 1;



        /* initialize prefix-sum table */
        alphabetSize = referenceSequenceIndex.alphabetDictionary.getAlphabetSize();
        referenceSequenceIndex.prefixSumTable = new int[alphabetSize + 2];
        referenceSequenceIndex.prefixSumTable[0] = 0;
        referenceSequenceIndex.prefixSumTable[alphabetSize + 1] = compactCodeBufferSize;    // last entry contains reference sequence length

        /* prepare start sequence buffer */
        startSequenceBuffer = new byte[START_SEQUENCE_BUFFER_SIZE];
        Arrays.fill(startSequenceBuffer, (byte) 1);

        /* prepare rotation buffers array */
        rotationBuffers = new int[START_SEQUENCE_BUFFER_SIZE + 1][];

        /* detect long runs of same character */
        detectLongRuns();

        /* calculate initial start sequence length */
        initial_start_sequence_length = INITIAL_START_SEQUENCE_LENGTH_5;
        if (alphabetSize > 5) {
            initial_start_sequence_length = INITIAL_START_SEQUENCE_LENGTH_10;
        }
        if (alphabetSize > 10) {
            initial_start_sequence_length = INITIAL_START_SEQUENCE_LENGTH_10_PLUS;
        }

        /* update progress bar */
        noStartSequenceCombinations = ((int) round(pow( alphabetSize,
                                                        initial_start_sequence_length   )));
        currentStartSequenceCombination = 0;
        java.awt.EventQueue.invokeLater(new Runnable() {
            @Override
            public void run() {
                parentFrame.setProgressBarMaximum(noStartSequenceCombinations);
                parentFrame.setProgressBarValue(currentStartSequenceCombination);
                parentFrame.setProgressBarIndeterminate(false);
            }
        });

        /* update statistics */
        currentPreprocessingMemoryConsumption += compactCodeBufferSize;
        currentPreprocessingMemoryConsumption += START_SEQUENCE_BUFFER_SIZE;
        currentPreprocessingMemoryConsumption += ((START_SEQUENCE_BUFFER_SIZE + 1) * 8);
        if (currentPreprocessingMemoryConsumption > maximalPreprocessingMemoryConsumption) {
            maximalPreprocessingMemoryConsumption = currentPreprocessingMemoryConsumption;
        }

        /* process each start sequence */
        changeDepth = 1;
        do {

            /* fill prefix-sum table */
            if (changeDepth == 1) {
                referenceSequenceIndex.prefixSumTable[startSequenceBuffer[0]] = bwtBufferIndex;
            }

            /*
             * process sequences that are shorter than initial sequence length.
             * They are located at the end of reference sequence.
             */
            for (; changeDepth < initial_start_sequence_length; changeDepth++) {
                checkTail(changeDepth);
            }


            startSequenceString = Byte.toString(startSequenceBuffer[0]);
            for (int i = 1; i < initial_start_sequence_length; i++) {
                startSequenceString += (", " + Byte.toString(startSequenceBuffer[i]));
            }
            System.out.println("Starting sequence \"" + startSequenceString + "\" at: " + sdf.format(Calendar.getInstance().getTime()));


            /* proces all rotations starting with sequence from start sequence buffer */
            findRotations(initial_start_sequence_length, 0, rotationBuffers);

            /* update progress bar */
            currentStartSequenceCombination++;
            java.awt.EventQueue.invokeLater(new Runnable() {
                @Override
                public void run() {
                    parentFrame.setProgressBarValue(currentStartSequenceCombination);
                }
            });

            /* calculate next start sequence */
            for (changeDepth = initial_start_sequence_length; changeDepth > 0; changeDepth--) {
                startSequenceBuffer[changeDepth - 1]++;
                if (startSequenceBuffer[changeDepth - 1] <= alphabetSize) {
                    break;
                }
                startSequenceBuffer[changeDepth - 1] = 1;
            }
        } while (changeDepth > 0);



//        System.out.println("BWT test started at: " + sdf.format(Calendar.getInstance().getTime()));
//        testBwt();



        /* release unused objects */
        releaseLongRunMemory();
        compactCodeBuffer = null;
        startSequenceBuffer = null;
        rotationBuffers = null;
        currentPreprocessingMemoryConsumption -= compactCodeBufferSize;
        currentPreprocessingMemoryConsumption -= START_SEQUENCE_BUFFER_SIZE;
        currentPreprocessingMemoryConsumption -= ((START_SEQUENCE_BUFFER_SIZE + 1) * 8);

        bwtCreationTimeBuffer = sdf.format(Calendar.getInstance().getTime());
        System.out.println("Reference sequence BWT creation completed at: " + bwtCreationTimeBuffer);
    }

    /**
     * Function createIndex creates reference sequence index upon its
     * BW transform.
     * 
     * @throws ReportedException 
     */
    private void createIndex() throws ReportedException {
        referenceSequenceIndex.occurrenceWaveletTree = new FMIndexWaveletTreeNode(bwtBuffer, referenceSequenceIndex.prefixSumTable);

        /* estimate occurence wavelet node memory consumption */
        currentPreprocessingMemoryConsumption += referenceSequenceIndex.occurrenceWaveletTree.getMemoryConsumption();
        if (currentPreprocessingMemoryConsumption > maximalPreprocessingMemoryConsumption) {
            maximalPreprocessingMemoryConsumption = currentPreprocessingMemoryConsumption;
        }

        /* release unused objects */
        currentPreprocessingMemoryConsumption -= bwtBuffer.length;
        bwtBuffer = null;

        indexCreationTimeBuffer = sdf.format(Calendar.getInstance().getTime());
        System.out.println("Reference sequence index creation completed at: " + indexCreationTimeBuffer);
    }

    /**
     * Function findRotations collects all reference sequence rotations beginning
     * with start sequence from startSequenceBuffer, sorts them and fills BWT
     * buffer. Function can be called recursively to increase start sequence
     * depth and to further reduce memory occupation by lowering number of
     * rotations to sort in one step.
     * 
     * @param startSequenceLength       depth of start sequence
     * @param startRotation             rotation to start search with
     * @param rotationBuffers           rotation buffers from previous recursions
     * @throws ReportedException 
     */
    private void findRotations( int startSequenceLength,
                                int startRotation,
                                int[][] rotationBuffers     ) throws ReportedException {

        int[] rotationBuffer;
        int rotationBufferIndex;
        int[] previousRotationBuffer;
        int previousRotationBufferSize;
        int referenceSeqIndex;
        int referenceSequenceSize;
        int totalRotations;
        int rotation;

        /* initialize rotation buffer */
        checkCancel();
        rotationBuffer = new int[ROTATION_BUFFER_SIZE];
        rotationBuffers[startSequenceLength] = rotationBuffer;
        rotationBufferIndex = 0;
        currentPreprocessingMemoryConsumption += (ROTATION_BUFFER_SIZE * 4);
        if (currentPreprocessingMemoryConsumption > maximalPreprocessingMemoryConsumption) {
            maximalPreprocessingMemoryConsumption = currentPreprocessingMemoryConsumption;
        }

        /* check the rest of reference sequence for matching start sequences */
        referenceSequenceSize = (compactCodeBuffer.length - startSequenceLength);
        referenceSeqIndex = startRotation;
        while ((referenceSeqIndex < referenceSequenceSize) && (rotationBufferIndex < ROTATION_BUFFER_SIZE)) {
            if (rotationEquals(referenceSeqIndex, startSequenceLength)) {
                rotationBuffer[rotationBufferIndex++] = referenceSeqIndex;
            }
            referenceSeqIndex++;
        }
        checkCancel();

        /* check if all matching sequences from reference sequence are now in rotation buffer */
        if (referenceSeqIndex >= referenceSequenceSize) {

            /* all matching rotations are taken from reference sequence buffer */

            /* count matching rotations from current level and all previous levels */
            totalRotations = rotationBufferIndex;
            for (int i = 1; i < startSequenceLength; i++) {
                previousRotationBuffer = rotationBuffers[i];
                if (previousRotationBuffer != null) {
                    previousRotationBufferSize = previousRotationBuffer.length;
                    for (int j = 0; j < previousRotationBufferSize; j++) {
                        rotation = previousRotationBuffer[j];
                        if (rotation >= 0) {
                            if (rotationEquals(rotation, startSequenceLength)) {
                                totalRotations++;
                            }
                        }
                    }
                }
            }
            checkCancel();

            /* check if there are matching rotations in previous levels */
            if (totalRotations != rotationBufferIndex) {

                /* there are matching rotations in previous levels */

                int[] sortBuffer;
                int sortBufferIndex;
                
                /* move matching rotations from all previous levels to sort buffer */
                sortBuffer = new int[totalRotations];
                currentPreprocessingMemoryConsumption += (totalRotations * 4);
                if (currentPreprocessingMemoryConsumption > maximalPreprocessingMemoryConsumption) {
                    maximalPreprocessingMemoryConsumption = currentPreprocessingMemoryConsumption;
                }
                sortBufferIndex = 0;
                for (int i = 1; i < startSequenceLength; i++) {
                    previousRotationBuffer = rotationBuffers[i];
                    if (previousRotationBuffer != null) {
                        previousRotationBufferSize = previousRotationBuffer.length;
                        for (int j = 0; j < previousRotationBufferSize; j++) {
                            rotation = previousRotationBuffer[j];
                            if (rotation >= 0) {
                                if (rotationEquals(rotation, startSequenceLength)) {
                                    sortBuffer[sortBufferIndex++] = rotation;
                                    previousRotationBuffer[j] = -1;
                                }
                            }
                        }
                    }
                }
                checkCancel();

                /* move rotations from current rotation buffer to sort buffer */
                for (int i = 0; i < rotationBufferIndex; i++) {
                    sortBuffer[sortBufferIndex++] = rotationBuffer[i];
                }

                /* all matching rotations from reference sequence are now in sort buffer */
                /* sort rotations */
                sortRotations(sortBuffer, totalRotations, startSequenceLength);

                /* insert sorted rotations into BWT buffer */
                insertBwt(sortBuffer, totalRotations);

                /* release sort buffer */
                sortBuffer = null;
                currentPreprocessingMemoryConsumption -= (totalRotations * 4);
            } else {

                /* there is no matching rotations in previous levels */

                /* all matching rotations from reference sequence are in rotation buffer */
                if (rotationBufferIndex > 0) {

                    /* sort rotations */
                    if (rotationBufferIndex > 1) {
                        sortRotations(rotationBuffer, rotationBufferIndex, startSequenceLength);
                    }

                    /* insert sorted rotations into BWT buffer */
                    insertBwt(rotationBuffer, rotationBufferIndex);
                }
            }
        } else {
 
            /* not all matching rotations from reference sequence are taken */
            /* increase start sequence length and do recursion for each possible extension */
            if (startSequenceLength >= START_SEQUENCE_BUFFER_SIZE) {
                reportError("Referentna sekvenca je prevelika!");
                throw new ReportedException();
            }
            checkTail(startSequenceLength);
            for (byte i = 1; i <= referenceSequenceIndex.alphabetDictionary.getAlphabetSize(); i++) {
                startSequenceBuffer[startSequenceLength] = i;
                findRotations(startSequenceLength + 1, referenceSeqIndex, rotationBuffers);
            }
        }

        /* release unused objects and update statistics */
        rotationBuffers[startSequenceLength] = null;
        currentPreprocessingMemoryConsumption -= (ROTATION_BUFFER_SIZE * 4);
    }

    /**
     * Function sortRotations sorts given rotations using heap sort algorithm
     * 
     * @param rotationBuffer        buffer containing rotations to sort
     * @param noRotations           number of rotations to sort
     * @param startSequenceLength   number of characters guaranteed to be the same
     *                              at the beginning of each rotation 
     * @throws ReportedException 
     */
    private void sortRotations( int[] rotationBuffer,
                                int noRotations,
                                int startSequenceLength ) throws ReportedException {

        int Child;
        int Parent;
        int CurrentParent;
        int PreviousParent;
        int LeftChild;
        int RightChild;
        int tmp;
        int cnt;

        /* heapify rotations */
        Parent = 0;
        Child = 1;
        while(Child < noRotations) {
            checkCancel();
            if (isRotationGreater(rotationBuffer[Child], rotationBuffer[Parent], startSequenceLength)) {
                tmp = rotationBuffer[Child];
                PreviousParent = Child;
                CurrentParent = Parent;

                do {
                    rotationBuffer[PreviousParent] = rotationBuffer[CurrentParent];
                    PreviousParent = CurrentParent;
                    if (CurrentParent < 1) {
                        break;
                    }
                    CurrentParent = ((PreviousParent - 1) / 2);
                } while (isRotationGreater(tmp, rotationBuffer[CurrentParent], startSequenceLength));

                rotationBuffer[PreviousParent] = tmp;
            }

            Child++;
            Parent = ((Child - 1) / 2);
        }

        /* sort rotations */
        cnt = noRotations - 1;
        while(cnt > 0) {
            checkCancel();
            tmp = rotationBuffer[cnt];
            rotationBuffer[cnt] = rotationBuffer[0];

            Parent = 0;
            LeftChild = 1;
            RightChild = 2;

            while (true) {
                if (RightChild < cnt) {
                    if (isRotationGreater(rotationBuffer[LeftChild], rotationBuffer[RightChild], startSequenceLength)) {
                        if (isRotationGreater(rotationBuffer[LeftChild], tmp, startSequenceLength)) {
                            rotationBuffer[Parent] = rotationBuffer[LeftChild];
                            Parent = LeftChild;
                        } else {
                            rotationBuffer[Parent] = tmp;
                            break;
                        }
                    } else {
                        if (isRotationGreater(rotationBuffer[RightChild], tmp, startSequenceLength)) {
                            rotationBuffer[Parent] = rotationBuffer[RightChild];
                            Parent = RightChild;
                        } else {
                            rotationBuffer[Parent] = tmp;
                            break;
                        }
                    }
                } else {
                    if (LeftChild < cnt) {
                        if (isRotationGreater(rotationBuffer[LeftChild], tmp, startSequenceLength)) {
                            rotationBuffer[Parent] = rotationBuffer[LeftChild];
                            rotationBuffer[LeftChild] = tmp;
                            break;
                        } else {
                            rotationBuffer[Parent] = tmp;
                            break;
                        }
                    } else {
                        rotationBuffer[Parent] = tmp;
                        break;
                    }
                }

                LeftChild = (Parent * 2) + 1;
                RightChild = LeftChild + 1;
            }
            cnt--;
        }
    }

    /**
     * Function insertBwt converts rotation number to BWT character and inserts it into BWT buffer.
     * 
     * @param sortBuffer    buffer containing sorted rotations
     */
    private void insertBwt(int[] sortBuffer, int noRotations) {

        int rotation;

        for (int i = 0; i < noRotations; i++) {
            rotation = sortBuffer[i];



//            testBwtBuffer[bwtBufferIndex] = rotation;



            if (rotation == 0) {
                bwtBuffer[bwtBufferIndex++] = 0;
            } else {
                bwtBuffer[bwtBufferIndex++] = compactCodeBuffer[rotation - 1];
            }
        }
    }

    /**
     * Function checkTail checks if sequence in start sequence buffer matches
     * sequence at the end of reference sequence buffer, just before terminating
     * character. If so, this is the next rotation in BW transformation and the
     * associated character is placed into BWT buffer.
     * 
     * @param tailSize      the length of sequence in start sequence buffer
     */
    private void checkTail(int tailSize) {

        int startSequenceIndex;
        int referenceSeqIndex;

        startSequenceIndex = 0;
        referenceSeqIndex = ((compactCodeBuffer.length - 1) - tailSize);
        for (int i = tailSize; i > 0; i--) {
            if (startSequenceBuffer[startSequenceIndex++] != compactCodeBuffer[referenceSeqIndex++]) {
                return;
            }
        }
        referenceSeqIndex = ((compactCodeBuffer.length - 2) - tailSize);
        bwtBuffer[bwtBufferIndex++] = compactCodeBuffer[referenceSeqIndex];



//        testBwtBuffer[bwtBufferIndex - 1] = ((compactCodeBuffer.length - 1) - tailSize);
    }

    /**
     * Function detectLongRuns detects long sequences of same character in reference
     * sequence and prepares data structures to efficiently sort reference sequence
     * rotations.
     * 
     * @throws ReportedException 
     */
    private void detectLongRuns() throws ReportedException {

        ArrayList<Integer>[] start;
        ArrayList<Integer>[] length;
        int referenceSequenceLength;
        byte startCharacter;
        int startRotation;
        int noChar;
        int count;

        /* create empty data structures */
        noChar = referenceSequenceIndex.alphabetDictionary.getAlphabetSize();
        longRunCount = new int[noChar + 1];
        longRunStart = new int[noChar + 1][];
        longRunLength = new int[noChar + 1][];
        currentPreprocessingMemoryConsumption += ((noChar + 1) * 12);
        if (currentPreprocessingMemoryConsumption > maximalPreprocessingMemoryConsumption) {
            maximalPreprocessingMemoryConsumption = currentPreprocessingMemoryConsumption;
        }
        start = new ArrayList[noChar + 1];
        length = new ArrayList[noChar + 1];
        for (int i = 1; i <= noChar; i++) {
            start[i] = new ArrayList();
            length[i] = new ArrayList();
        }

        /* find long run sequences */
        referenceSequenceLength = compactCodeBuffer.length;
        startCharacter = 0;
        startRotation = 0;
        count = 1;
        for (int i = 0; i < referenceSequenceLength; i++) {
            if (compactCodeBuffer[i] == startCharacter) {
                count++;
            } else {
                if (count >= LONG_RUN_SEQUENCE_THRESHOLD) {
                    start[startCharacter].ensureCapacity(longRunCount[startCharacter] + 1);
                    start[startCharacter].add(startRotation);
                    length[startCharacter].ensureCapacity(longRunCount[startCharacter] + 1);
                    length[startCharacter].add(count);
                    longRunCount[startCharacter]++;
                }
                startCharacter = compactCodeBuffer[i];
                startRotation = i;
                count = 1;
            }
        }
        if (count >= LONG_RUN_SEQUENCE_THRESHOLD) {
            start[startCharacter].ensureCapacity(longRunCount[startCharacter] + 1);
            start[startCharacter].add(startRotation);
            length[startCharacter].ensureCapacity(longRunCount[startCharacter] + 1);
            length[startCharacter].add(count);
            longRunCount[startCharacter]++;
        }

        /* move found long run sequences to reference sequence object fields */
        for (int i = 1; i <= noChar; i++) {
            if (longRunCount[i] > 0) {
                longRunStart[i] = new int[longRunCount[i]];
                longRunLength[i] = new int[longRunCount[i]];
                currentPreprocessingMemoryConsumption += (longRunCount[i] * 8);
                if (currentPreprocessingMemoryConsumption > maximalPreprocessingMemoryConsumption) {
                    maximalPreprocessingMemoryConsumption = currentPreprocessingMemoryConsumption;
                }
                for (int j = 0; j < longRunCount[i]; j++) {
                    longRunStart[i][j] = start[i].get(j);
                    longRunLength[i][j] = length[i].get(j);
                }
            }
        }
    }

    /**
     * Function releaseLongRunMemory releases memory used by long run data
     * structures.
     * 
     * @throws ReportedException 
     */
    private void releaseLongRunMemory() throws ReportedException {

        int noChar;

        noChar = referenceSequenceIndex.alphabetDictionary.getAlphabetSize();
        for (int i = 1; i <= noChar; i++) {
            if (longRunCount[i] > 0) {
                longRunStart[i] = null;
                longRunLength[i] = null;
                currentPreprocessingMemoryConsumption -= (longRunCount[i] * 8);
            }
        }
    }

    /**
     * Function rotationEquals checks if given rotation of reference sequence
     * begins with sequence in start sequence buffer.
     * 
     * @param rotation      rotation to check against start sequence buffer
     * @param length        length of sequence
     * @return              true if given rotation matches the rotation from start sequence buffer
     */
    private boolean rotationEquals(int rotation, int length) {

        int startSequenceIndex;

        for (startSequenceIndex = 0; startSequenceIndex < length; startSequenceIndex++) {
            if (startSequenceBuffer[startSequenceIndex] != compactCodeBuffer[rotation++]) {
                return false;
            }
        }
        return true;
    }

    /**
     * Function isRotationGreater checks if rotation A of reference sequence is
     * greater than rotation B. Rotations are guaranteed to be equal at first
     * startSequenceLength places.
     * 
     * @param rotationA             rotation A
     * @param rotationB             rotation B
     * @param startSequenceLength   number of places at the beginning that are guaranteed to be equal
     * @return                      true if rotation A is greater than rotation B
     */
    private boolean isRotationGreater(int rotationA, int rotationB,  int startSequenceLength) {

        int charactersToSkip;

        rotationA += startSequenceLength;
        rotationB += startSequenceLength;
        for (int i = 0; i < TEST_LONG_RUN_SEQUENCE_AFTER; i++) {
            if (compactCodeBuffer[rotationA++] != compactCodeBuffer[rotationB++]) {
                if (compactCodeBuffer[--rotationA] > compactCodeBuffer[--rotationB]) {
                    return true;
                }
                return false;
            }
        }

        /* check if comparing two long run rotations */
        charactersToSkip = checkLongRunSequence(rotationA, rotationB);
        rotationA += charactersToSkip;
        rotationB += charactersToSkip;

        while (compactCodeBuffer[rotationA++] == compactCodeBuffer[rotationB++]) {
        }
        if (compactCodeBuffer[--rotationA] > compactCodeBuffer[--rotationB]) {
            return true;
        }
        return false;
    }

    /**
     * Function checkLongRunSequence checks if both input rotations are inside
     * long run sequence for the same character, and calculate for how many
     * following characters the rotations will stay equal.
     * 
     * @param rotationA     rotation A
     * @param rotationB     rotation B
     * @return              number of characters rotation A and B will stay equal
     */
    private int checkLongRunSequence(int rotationA, int rotationB) {

        int characterCode;
        int rotationALength;
        int rotationBLength;
        int sequenceCount;
        int sequenceIndex;
        int sequenceIndexLow;
        int sequenceIndexHigh;
        int longRunSequenceStart;
        int longRunSequenceLength;

        /* check if characters are equal and if character has long run sequences */
        characterCode = compactCodeBuffer[rotationA];
        if (characterCode != compactCodeBuffer[rotationB]) {
            return 0;
        }
        sequenceCount = longRunCount[characterCode];
        if (sequenceCount <= 0) {
            return 0;
        }

        /* checks rotation A */
        sequenceIndexHigh = (sequenceCount - 1);
        sequenceIndexLow = 0;
        while (true) {
            sequenceIndex = ((sequenceIndexHigh + sequenceIndexLow) / 2);
            longRunSequenceStart = longRunStart[characterCode][sequenceIndex];
            longRunSequenceLength = longRunLength[characterCode][sequenceIndex];
            if (longRunSequenceStart > rotationA) {
                sequenceIndexHigh = (sequenceIndex - 1);
                if (sequenceIndexHigh < sequenceIndexLow) {
                    return 0;
                }
            } else {
                if ((longRunSequenceStart + longRunSequenceLength) > rotationA) {
                    break;
                } else {
                    sequenceIndexLow = (sequenceIndex + 1);
                    if (sequenceIndexHigh < sequenceIndexLow) {
                        return 0;
                    }
                }
            }
        }
        rotationALength = ((longRunSequenceLength) - (rotationA - longRunSequenceStart));

        /* checks rotation B */
        sequenceIndexHigh = (sequenceCount - 1);
        sequenceIndexLow = 0;
        while (true) {
            sequenceIndex = ((sequenceIndexHigh + sequenceIndexLow) / 2);
            longRunSequenceStart = longRunStart[characterCode][sequenceIndex];
            longRunSequenceLength = longRunLength[characterCode][sequenceIndex];
            if (longRunSequenceStart > rotationB) {
                sequenceIndexHigh = (sequenceIndex - 1);
                if (sequenceIndexHigh < sequenceIndexLow) {
                    return 0;
                }
            } else {
                if ((longRunSequenceStart + longRunSequenceLength) > rotationB) {
                    break;
                } else {
                    sequenceIndexLow = (sequenceIndex + 1);
                    if (sequenceIndexHigh < sequenceIndexLow) {
                        return 0;
                    }
                }
            }
        }
        rotationBLength = ((longRunSequenceLength) - (rotationB - longRunSequenceStart));

        return min(rotationALength, rotationBLength);
    }

    /**
     * Function isSequenceName reads next line from reference sequence Fasta
     * file and checks if it is sequence start line for sequence given by
     * sequenceNameBuffer parameter.
     * 
     * @param inStream              FileInputStream object associated with reference sequence file.
     * @param sequenceNameBuffer    sequence name buffer.
     * @param fileLength            Length of the rest of the reference sequence file.
     * @return                      returns true if current line is given sequence start line.
     * @throws ReportedException 
     */
    private boolean isSequenceName(byte[] sequenceNameBuffer) throws ReportedException {

        int nextChar;
        int sequenceNameLength = sequenceNameBuffer.length;

        for (int i = 0; i < sequenceNameLength; i++) {
            nextChar = getCharacter();
            if (nextChar != sequenceNameBuffer[i]) {
                if (nextChar != NEW_LINE_CHARACTER) {
                    discardLine();
                }
                return false;
            }
        }
        nextChar = getCharacter();
        if (nextChar == NEW_LINE_CHARACTER) {
            return true;
        }
        discardLine();
        if (nextChar != SPACE_CHARACTER) {
            return false;
        }
        return true;
    }

    /**
     * Function discardLine discards all characters of current line in reference
     * sequence file
     *
     * @throws ReportedException 
     */
    private void discardLine() throws ReportedException {

        int nextChar;

        do {
            nextChar = getCharacter();
            if (nextChar == NEW_LINE_CHARACTER) {
                return;
            }
        } while (nextChar >= 0);
    }

    /**
     * Function getCharacter provides next character from reference sequence
     * intermediate buffer. If intermediate buffer is empty, it is filled up
     * from reference sequence file before the character is returned.
     * 
     * @return      return next character or -1 if whole file is red.
     * @throws ReportedException 
     */
    private int getCharacter() throws ReportedException {

        int charactersRead;
        int charactersRequested;

        if (intermediateReadBufferIndex >= READ_BUFFER_SIZE) {
            if (referenceSequenceNotRead <= 0) {
                return -1;
            }
            intermediateReadBufferIndex = 0;
            if (referenceSequenceNotRead < READ_BUFFER_SIZE) {
                intermediateReadBufferIndex = ((int) (READ_BUFFER_SIZE - referenceSequenceNotRead));
            }
            charactersRequested = (READ_BUFFER_SIZE - intermediateReadBufferIndex);
            try {
                charactersRead = referenceSequenceStream.read(  intermediateReadBuffer,
                                                                intermediateReadBufferIndex,
                                                                charactersRequested             );
            } catch (Exception ex) {
                reportError(    "Greška tijekom čitanja datoteke referentnog slijeda ("
                                + getExceptionType(ex)
                                + "):\n\n    "
                                + ex.getMessage()                                           );
                throw new ReportedException();
            }
            if (charactersRead != charactersRequested) {
                reportError("Nepoznata greška tijekom čitanja datoteke referentnog slijeda!");
                throw new ReportedException();
            }
            referenceSequenceNotRead -= charactersRead;
        }

        return (intermediateReadBuffer[intermediateReadBufferIndex++] & 0x000000ff);
    }

    /**
     * Function isEndOfFile checks if end of reference sequence file is reached.
     * 
     * @return      returns true if end of file is reached
     */
    private boolean isEndOfFile() {
        if (intermediateReadBufferIndex >= READ_BUFFER_SIZE) {
            if (referenceSequenceNotRead <= 0) {
                return true;
            }
        }
        return false;
    }

    /**
     * Function toUpperCase converts given character to upper case
     * 
     * @param inChar        character to convert to upper case
     * @return              upper case of given character
     */
    private int toUpperCase (int inChar) {
        if ((inChar >= LOWRCASE_A_CHARACTER) && (inChar <= LOWRCASE_Z_CHARACTER)) {
            return (inChar & (~0x20));
        }
        return inChar;
    }

    /**
     * Function getExceptionType provides the type of exception
     * 
     * @param   ex    exception which type is to be determined
     */
    private String getExceptionType (Exception ex) {

        String exceptionFullName;
        int lastDotIndex;

        exceptionFullName = ex.getClass().toString();
        lastDotIndex = exceptionFullName.lastIndexOf('.');
        if (lastDotIndex < 0) {
            lastDotIndex = 0;
        } else {
            lastDotIndex++;
        }
        return exceptionFullName.substring(lastDotIndex);
    }

    /**
     * Function checkCancel checks if reference sequence creation cancellation
     * is requested. If so, confirmation window is opened and creation is
     * interrupted.
     * 
     * @throws ReportedException 
     */
    private void checkCancel() throws ReportedException {
        if (Thread.interrupted()) {
            reportError("Kreiranje indeksa referentnog slijeda prekinuto na zahtjev korisnika!");
            throw new ReportedException();
        }
    }

    /**
     * Function saveReferenceSequenceIndex saves preprocessed reference sequence
     * index.
     * 
     * @param filename      filename to save reference sequence index to
     * @throws ReportedException 
     */
    public void saveReferenceSequenceIndex(String filename) throws ReportedException {

        ObjectOutputStream out = null;

        try {
            out = new ObjectOutputStream(new FileOutputStream(filename));
            try {
                out.writeObject(referenceSequenceIndex);
            } catch (Exception ex) {
                reportError(    "Greška tijekom spremanja predprocesirane datoteke referentnog slijeda ("
                                + getExceptionType(ex)
                                + "):\n\n    "
                                + ex.getMessage()                                           );
                throw new ReportedException();
            }
        } catch (ReportedException ex) {
            throw ex;
        } catch (Exception ex) {
            reportError(    "Greška tijekom otvaranja datoteke za spremanje predprocesiranog referentnog slijeda ("
                            + getExceptionType(ex)
                            + "):\n\n    "
                            + ex.getMessage()                                           );
            throw new ReportedException();
        } finally {
            if (out != null) {
                try {
                    out.close();
                } catch (IOException ex) {
                }
            }
        }

        JOptionPane.showConfirmDialog(null, "Referentni slijed je uspješno spremljen!", "FMIndexCount", JOptionPane.DEFAULT_OPTION, JOptionPane.INFORMATION_MESSAGE);
    }

    /**
     * Function countQuerySequenceOccurrences counts occurrences of given query
     * sequence in reference sequence.
     * 
     * @param qurySeq       string representing query sequence
     * @return              number of occurrences of query sequence in reference sequence
     * @throws ReportedException 
     */
    public int countQuerySequenceOccurrences (String qurySeq) throws ReportedException {

        Charset utf8Charset;
        Calendar calendar;
        String queryStartTimeBuffer;
        String queryEndTimeBuffer;
        String statistics;
        byte[] querySequenceBuffer;
        long countDuration;
        int queryCharacter;
        int queryIndex;
        int queryLength;
        int lowIndex;
        int lowIndexShift;
        int highIndex;
        int highIndexShift;
        byte characterToProcess;

        /* get query start time */
        calendar = Calendar.getInstance();
        countDuration = calendar.getTimeInMillis();
        queryStartTimeBuffer = sdf.format(calendar.getTime());
        System.out.println("Query start time: " + queryStartTimeBuffer);

        /* prepare query buffer */
        queryLength = qurySeq.length();
        if (queryLength <= 0) {
            reportError("Prazan upitni slijed!");
            throw new ReportedException();
        }
        querySequenceBuffer = new byte[queryLength];
        utf8Charset = Charset.forName("UTF-8");
        querySequenceBuffer = qurySeq.getBytes(utf8Charset);
        utf8Charset = null;
        queryIndex = 0;
        for (int i = 0; i < queryLength; i++) {
            queryCharacter = querySequenceBuffer[i];
            if ((queryCharacter >= 32) && (queryCharacter <= 126)) {
                queryCharacter = toUpperCase(queryCharacter);
                queryCharacter = referenceSequenceIndex.alphabetDictionary.getCompactCode((byte) queryCharacter);
                if (queryCharacter == 0) {
                    reportError("Znak (" + qurySeq.substring(i, i + 1) + ") iz upitnog slijeda ne nalazi se u referentnom slijedu!");
                    throw new ReportedException();
                }
                querySequenceBuffer[queryIndex++] = ((byte) queryCharacter);
            }
        }
        queryLength = queryIndex;
        if (queryLength <= 0) {
            reportError("Prazan upitni slijed!");
            throw new ReportedException();
        }

        /* count matches */
        characterToProcess = querySequenceBuffer[queryLength - 1];
        lowIndex =  referenceSequenceIndex.C(characterToProcess);
        highIndex = (   referenceSequenceIndex.C((byte) (characterToProcess + 1)) - 1);
        for (int i = (queryLength - 2); i >= 0; i--) {
            characterToProcess = querySequenceBuffer[i];

            lowIndexShift = referenceSequenceIndex.occurrenceWaveletTree.Occ(characterToProcess, lowIndex);
            lowIndex =  referenceSequenceIndex.C(characterToProcess);
            lowIndex += lowIndexShift;

            highIndexShift = referenceSequenceIndex.getCharacterTotal(characterToProcess);
            highIndexShift -= referenceSequenceIndex.occurrenceWaveletTree.Occ(characterToProcess, (highIndex + 1));
            highIndex = (   referenceSequenceIndex.C((byte) (characterToProcess + 1)) - 1);
            highIndex -= highIndexShift;

            if (highIndex < lowIndex) {
                break;
            }
        }
        highIndex -= lowIndex;
        highIndex++;

        /* get query end time */
        calendar = Calendar.getInstance();
        countDuration = (calendar.getTimeInMillis() - countDuration);
        queryEndTimeBuffer = sdf.format(calendar.getTime());
        System.out.println("Query start time: " + queryEndTimeBuffer);

        /* report statistics */
        statistics = "Statistika brojanja ponavljanja upita u referentnoj sekvenci:\n";
        statistics += (  "\n    Vrijeme početka brojanja:                           " + queryStartTimeBuffer);
        statistics += (  "\n    Vrijeme završetka brojanja:                         " + queryEndTimeBuffer);
        statistics += (  "\n    Ukupno vrijeme brojanja:                            " + Long.toString(countDuration) + " ms");

        statistics += ("\n\n    Dužina referentnog slijede:                         " + Integer.toString(referenceSequenceIndex.prefixSumTable[referenceSequenceIndex.prefixSumTable.length - 1] - 1));
        statistics += (  "\n    Dužina upitnog slijede:                             " + Integer.toString(queryLength));

        statistics += ("\n\n    Broj pronađenih podudaranja:                        " + Integer.toString(highIndex));

        statistics += ("\n\n    Maksimalno zauzeće memorije tijekom brojanja:       " + Integer.toString((int) ((maximalPreprocessingMemoryConsumption + querySequenceBuffer.length) / 1000)) + " KB");

        JOptionPane.showConfirmDialog(null, statistics, "FMIndexCount", JOptionPane.DEFAULT_OPTION, JOptionPane.INFORMATION_MESSAGE);

        if (highIndex <= 0) {
            return 0;
        } else {
            return highIndex;
        }
    }
    
    /** Function reportStatistics reports reference sequence creation statistics */
    private void reportStatistics() {

        String statistics;

        statistics = "Statistika kreiranja indeksa referentne sekvence:\n";
        statistics += (     "\n    Vrijeme početka kreiranja:                                    " + startTimeBuffer);
        if (fileReadTimeBuffer != null) {
            statistics += ( "\n    Vrijeme završetka učitavanja datoteke referentnog slijeda:    " + fileReadTimeBuffer);
        }
        if (bwtCreationTimeBuffer != null) {
            statistics += ( "\n    Vrijeme završetka kreiranja BW transformacije:                " + bwtCreationTimeBuffer);
        }
        statistics += (     "\n    Vrijeme završetka kreiranja:                                  " + indexCreationTimeBuffer);

        statistics += (   "\n\n    Dužina referentnog slijede:                                   " + Integer.toString(referenceSequenceIndex.prefixSumTable[referenceSequenceIndex.prefixSumTable.length - 1] - 1));

        statistics += (   "\n\n    Maksimalno zauzeće memorije tijekom kreiranja:                " + Integer.toString((int) (maximalPreprocessingMemoryConsumption / 1000)) + " KB");

        JOptionPane.showConfirmDialog(null, statistics, "FMIndexCount", JOptionPane.DEFAULT_OPTION, JOptionPane.INFORMATION_MESSAGE);
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

    
    
    // <editor-fold desc="Test functions">
    
    private void testBwt() {
//
//        int bwtBufferSize;
//        int rotation;
//        byte bwt;
//
//        bwtBufferSize = bwtBuffer.length;
//        for (int i = 1; i < bwtBufferSize; i++) {
//            if (isRotationGreater(testBwtBuffer[i - 1], testBwtBuffer[i],  0)) {
////            if (testIsRotationGreater(testBwtBuffer[i - 1], testBwtBuffer[i])) {
//                reportError("Rotation sort error!");
//            }
//
//            rotation = testBwtBuffer[i];
//            if (rotation == 0) {
//                bwt = 0;
//            } else {
//                bwt = compactCodeBuffer[rotation - 1];
//            }
//            if (bwtBuffer[i] != bwt) {
//                reportError("BWT sort error!");
//            }
//        }
//        
    }

    private boolean testIsRotationGreater(int rotationA, int rotationB) {
        while (compactCodeBuffer[rotationA++] == compactCodeBuffer[rotationB++]) {
        }
        if (compactCodeBuffer[--rotationA] > compactCodeBuffer[--rotationB]) {
            return true;
        }
        return false;
    }
    // </editor-fold>



    // <editor-fold desc="Enumerations">
    public enum ReferenceSequenceFormat {
		FASTA,
		TEXT,
		PREPROCESSED
	}
    // </editor-fold>
}
