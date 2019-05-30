package edu.ucr.lxu051;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

public class MultipleSequenceAlignment {

    public String alpha, beta, gamma;
    public String alphaR, betaR, gammaR;
    public StringBuilder a, b, c;
    public int[][] scoreFunction;

    
    public MultipleSequenceAlignment(String alpha, String beta, String gamma) {
        this.a = new StringBuilder();
        this.b = new StringBuilder();
        this.c = new StringBuilder();
        this.alpha = alpha;
        this.beta = beta;
        this.gamma = gamma;
        this.alphaR = reverseString(alpha);
        this.betaR = reverseString(beta);
        this.gammaR = reverseString(gamma);
        this.scoreFunction = genBlastScoreFunction();
    }

    public void align() {
        //dp(alpha, beta, gamma, scoreFunction);

        if (alpha.length() < 2 || beta.length() < 2 || gamma.length() < 2) {
            dp(alpha, beta, gamma);
        } else {
            path(0, 0, 0, alpha.length(), beta.length(), gamma.length());
        }

//        System.out.println(a);
//        System.out.println(b);
//        System.out.println(c);
        int[] stats = calFinalScore();
        System.out.println(stats[0]);
        System.out.println(stats[1]);
    }

    public static void main(String[] args) throws FileNotFoundException {
        Scanner sc = new Scanner(new File("test1.txt"));
        String alpha = sc.nextLine();
        String beta = sc.nextLine();
        String gamma = sc.nextLine();
        long startTime = System.currentTimeMillis();
        MultipleSequenceAlignment msa = new MultipleSequenceAlignment(alpha, beta, gamma);
        msa.align();
        long endTime = System.currentTimeMillis();
        long time = (endTime - startTime) / 1000;
        System.out.println("Running time: " + time/60 + " min & " + time % 60 + " sec.");
        System.out.println(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());
    }

    public void dp(String alpha, String beta, String gamma) {
        int[][][] matrix = new int[alpha.length()+1][beta.length()+1][gamma.length()+1];
        for (int i = 0; i < alpha.length()+1; i++) {
            for (int j = 0; j < beta.length()+1; j++) {
                for (int k = 0; k < gamma.length()+1; k++) {
                    if (i == 0 && j == 0 && k == 0) {
                        matrix[i][j][k] = 0;
                    } else {
                        int one, two, three, four, five, six, seven;
                        try {
                            one = matrix[i - 1][j][k] + sigma(alpha.charAt(i - 1), '-', '-', scoreFunction);
                        } catch (IndexOutOfBoundsException e) {
                            one = -Integer.MAX_VALUE;
                        }
                        try {
                            two = matrix[i][j - 1][k] + sigma('-', beta.charAt(j - 1), '-', scoreFunction);
                        } catch (IndexOutOfBoundsException e) {
                            two = -Integer.MAX_VALUE;
                        }
                        try {
                            three = matrix[i][j][k - 1] + sigma('-', '-', gamma.charAt(k - 1), scoreFunction);
                        } catch (IndexOutOfBoundsException e) {
                            three = -Integer.MAX_VALUE;
                        }
                        try {
                            four = matrix[i - 1][j - 1][k] + sigma(alpha.charAt(i - 1), beta.charAt(j - 1), '-', scoreFunction);
                        } catch (IndexOutOfBoundsException e) {
                            four = -Integer.MAX_VALUE;
                        }
                        try {
                            five = matrix[i - 1][j][k - 1] + sigma(alpha.charAt(i - 1), '-', gamma.charAt(k - 1), scoreFunction);
                        } catch (IndexOutOfBoundsException e) {
                            five = -Integer.MAX_VALUE;
                        }
                        try {
                            six = matrix[i][j - 1][k - 1] + sigma('-', beta.charAt(j - 1), gamma.charAt(k - 1), scoreFunction);
                        } catch (IndexOutOfBoundsException e) {
                            six = -Integer.MAX_VALUE;
                        }
                        try {
                            seven = matrix[i - 1][j - 1][k - 1] + sigma(alpha.charAt(i - 1), beta.charAt(j - 1), gamma.charAt(k - 1), scoreFunction);
                        } catch (IndexOutOfBoundsException e) {
                            seven = -Integer.MAX_VALUE;
                        }
                        int[] scores = {one, two, three, four, five, six, seven};
                        matrix[i][j][k] = max(scores);
                    }
//                    System.out.print(matrix[i][j][k] + " ");
                }
//                System.out.println();
            }
//            System.out.println();
        }
        backtracing(alpha, beta, gamma, matrix);
    }

    public void backtracing(String alpha, String beta, String gamma, int[][][] matrix) {
        int i = matrix.length - 1;
        int j = matrix[0].length - 1;
        int k = matrix[0][0].length - 1;
        StringBuilder tempA = new StringBuilder();
        StringBuilder tempB = new StringBuilder();
        StringBuilder tempC = new StringBuilder();
        while (i > 0 || j > 0 || k > 0) {
            int cur = matrix[i][j][k];
            if ((i>0) && matrix[i-1][j][k] + sigma(alpha.charAt(i-1), '-', '-', scoreFunction) == cur) {
                tempA.append(alpha.charAt(i-1)); tempB.append('-'); tempC.append('-');
                i--;
            } else if ((j>0) && matrix[i][j-1][k] + sigma('-', beta.charAt(j-1), '-', scoreFunction) == cur) {
                tempA.append('-'); tempB.append(beta.charAt(j-1)); tempC.append('-');
                j--;
            } else if ((k>0) && matrix[i][j][k-1] + sigma('-', '-', gamma.charAt(k-1), scoreFunction) == cur) {
                tempA.append('-'); tempB.append('-'); tempC.append(gamma.charAt(k-1));
                k--;
            } else if ((i>0 && j>0) && matrix[i-1][j-1][k] + sigma(alpha.charAt(i-1), beta.charAt(j-1), '-', scoreFunction) == cur) {
                tempA.append(alpha.charAt(i-1)); tempB.append(beta.charAt(j-1)); tempC.append('-');
                i--; j--;
            } else if ((i>0 && k>0) && matrix[i-1][j][k-1] + sigma(alpha.charAt(i-1), '-', gamma.charAt(k-1), scoreFunction) == cur) {
                tempA.append(alpha.charAt(i-1)); tempB.append('-'); tempC.append(gamma.charAt(k-1));
                i--; k--;
            } else if ((j>0 && k>0) && matrix[i][j-1][k-1] + sigma('-', beta.charAt(j-1), gamma.charAt(k-1), scoreFunction) == cur) {
                tempA.append('-'); tempB.append(beta.charAt(j-1)); tempC.append(gamma.charAt(k-1));
                j--; k--;
            } else if ((i>0 && j>0 && k>0) && matrix[i-1][j-1][k-1] + sigma(alpha.charAt(i-1), beta.charAt(j-1), gamma.charAt(k-1), scoreFunction) == cur) {
                tempA.append(alpha.charAt(i-1)); tempB.append(beta.charAt(j-1)); tempC.append(gamma.charAt(k-1));
                i--; j--; k--;
            } else {
                System.out.println("Error");
            }
        }
        a.append(tempA.reverse());
        b.append(tempB.reverse());
        c.append(tempC.reverse());
    }

    public void path(int a1, int b1, int c1, int a2, int b2, int c2) {
//        System.out.println("A1 = " + a1 + "; B1 = " + b1 + "; C1 = " + c1 + "; A2 = " + a2 + "; B2 = " + b2 + "; C2 = " + c2);
        if (a1 + 1 >= a2 || b1 + 1 >= b2 || c1 + 1 >= c2) {
            dp(alpha.substring(a1, a2), beta.substring(b1, b2), gamma.substring(c1, c2));
            //System.out.println(alpha.substring(a1, a2) + " " + beta.substring(b1, b2) + " " + gamma.substring(c1, c2));
        } else {
            int mid = Math.floorDiv(c1 + c2, 2);
            int[][] forwardScores = maxPathScore(a1, b1, c1, a2, b2, mid, "f");
            int[][] backwardScores = maxPathScore(a1, b1, mid, a2, b2, c2, "b");

//            printMatrix(forwardScores);
//            System.out.println();
//            printMatrix(backwardScores);

            int j1 = 0, j2 = 0;
            int max = -Integer.MAX_VALUE;
            for (int i = 0; i < forwardScores.length; i++) {
                for (int j = 0; j < forwardScores[0].length; j++) {
                    int temp = forwardScores[i][j] + backwardScores[i][j];
                    if (temp > max) {
                        max = temp;
                        j1 = i;
                        j2 = j;
                    }
                }
            }
            //System.out.println("i = " + j1 + "; j = " + j2);
            path(a1, b1, c1, a1+j1, b1+j2, mid);
            path(a1+j1, b1+j2, mid, a2, b2, c2);
        }
    }

    public int[][] maxPathScore(int a1, int b1, int c1, int a2, int b2, int c2, String direction) {
        String alpha, beta, gamma;
        if (direction.equals("f")) {
            alpha = this.alpha;
            beta = this.beta;
            gamma = this.gamma;
        } else {
            alpha = this.alphaR;
            beta = this.betaR;
            gamma = this.gammaR;
            int tempA1 = this.alpha.length() - a2;
            int tempA2 = this.alpha.length() - a1;
            int tempB1 = this.beta.length() - b2;
            int tempB2 = this.beta.length() - b1;
            int tempC1 = this.gamma.length() - c2;
            int tempC2 = this.gamma.length() - c1;
            a1 = tempA1; a2 = tempA2; b1 = tempB1; b2 = tempB2; c1 = tempC1; c2 = tempC2;
        }
        int alphaSize = a2 - a1 + 1, betaSize = b2 - b1 + 1;
        int[][] tempScore = new int[alphaSize][betaSize];
        int[][] finalScore = new int[alphaSize][betaSize];
        boolean init = true;

        for (int i = 0; i < alphaSize; i++) {
            for (int j = 0; j < betaSize; j++) {
                if (i == 0 && j == 0 && init) {
                    tempScore[i][j] = 0;
                    init = false;
                } else {
                    int one, two, three;
                    try {
                        one = tempScore[i - 1][j] + sigma(alpha.charAt(a1 + i - 1), '-', '-', scoreFunction);
                    } catch (IndexOutOfBoundsException e) {
                        one = -Integer.MAX_VALUE;
                    }
                    try {
                        two = tempScore[i][j - 1] + sigma('-', beta.charAt(b1 + j - 1), '-', scoreFunction);
                    } catch (IndexOutOfBoundsException e) {
                        two = -Integer.MAX_VALUE;
                    }
                    try {
                        three = tempScore[i - 1][j - 1] + sigma(alpha.charAt(a1 + i - 1), beta.charAt(b1 + j - 1), '-', scoreFunction);
                    } catch (IndexOutOfBoundsException e) {
                        three = -Integer.MAX_VALUE;
                    }
                    int[] scores = {one, two, three};
                    tempScore[i][j] = max(scores);
                }
            }
        }

//        printMatrix(tempScore);

        for (int k = 0; k < c2 - c1 + 1; k++) {
            finalScore = new int[alphaSize][betaSize];
            for (int i = 0; i < alphaSize; i++) {
                for (int j = 0; j < betaSize; j++) {
                    if (k == 0) {
                        finalScore[i][j] = tempScore[i][j];
                    } else {
                        int one, two, three, four, five, six, seven;
                        try {
                            one = finalScore[i - 1][j] + sigma(alpha.charAt(a1 + i - 1), '-', '-', scoreFunction);
                        } catch (IndexOutOfBoundsException e) {
                            one = -Integer.MAX_VALUE;
                        }
                        try {
                            two = finalScore[i][j - 1] + sigma('-', beta.charAt(b1 + j - 1), '-', scoreFunction);
                        } catch (IndexOutOfBoundsException e) {
                            two = -Integer.MAX_VALUE;
                        }
                        try {
                            three = finalScore[i - 1][j - 1] + sigma(alpha.charAt(a1 + i - 1), beta.charAt(b1 + j - 1), '-', scoreFunction);
                        } catch (IndexOutOfBoundsException e) {
                            three = -Integer.MAX_VALUE;
                        }
                        try {
                            four = tempScore[i][j] + sigma('-', '-', gamma.charAt(c1 + k - 1), scoreFunction);
                        } catch (IndexOutOfBoundsException e) {
                            four = -Integer.MAX_VALUE;
                        }
                        try {
                            five = tempScore[i - 1][j] + sigma(alpha.charAt(a1 + i - 1), '-', gamma.charAt(c1 + k - 1), scoreFunction);
                        } catch (IndexOutOfBoundsException e) {
                            five = -Integer.MAX_VALUE;
                        }
                        try {
                            six = tempScore[i][j - 1] + sigma('-', beta.charAt(b1 + j - 1), gamma.charAt(c1 + k - 1), scoreFunction);
                        } catch (IndexOutOfBoundsException e) {
                            six = -Integer.MAX_VALUE;
                        }
                        try {
                            seven = tempScore[i - 1][j - 1] + sigma(alpha.charAt(a1 + i - 1), beta.charAt(b1 + j - 1), gamma.charAt(c1 + k - 1), scoreFunction);
                        } catch (IndexOutOfBoundsException e) {
                            seven = -Integer.MAX_VALUE;
                        }
                        int[] scores = {one, two, three, four, five, six, seven};
                        finalScore[i][j] = max(scores);
                    }
                }
            }
            for (int i = 0; i < tempScore.length; i++) {
                for (int j = 0; j < tempScore[0].length; j++) {
                    tempScore[i][j] = finalScore[i][j]; // soft copy
                }
            }
        }

        if (direction.equals("b")) {
            return reverseMatrix(finalScore);
        } else {
            return finalScore;
        }
    }

    public int max(int[] scores) {
        int max = -Integer.MAX_VALUE;
        for (int s : scores) {
            max = (s > max) ? s : max;
        }
        return max;
    }

    public int sigma(char a, char b, char c, int[][] scoreFunction) { //'A', 'C', 'G', 'T', '-'
        int indexA = getIndex(a);
        int indexB = getIndex(b);
        int indexC = getIndex(c);
        return scoreFunction[indexA][indexB] + scoreFunction[indexA][indexC] + scoreFunction[indexB][indexC];
    }

    public int getIndex(char c) {
        if (c == 'A' || c == 'a') {
            return 0;
        } else if (c == 'C' || c == 'c') {
            return 1;
        } else if (c == 'G' || c == 'g') {
            return 2;
        } else if (c == 'T' || c == 't') {
            return 3;
        } else {
            return 4;
        }
    }

    public String reverseString(String s) {
        StringBuilder reverse = new StringBuilder(s);
        return reverse.reverse().toString();
    }

    public int[][] reverseMatrix(int[][] matrix) {
        int[][] reversed = new int[matrix.length][matrix[0].length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                reversed[i][j] = matrix[matrix.length - i - 1][matrix[0].length - j - 1];
            }
        }
        return reversed;
    }

    // scoreFunction = 'A', 'C', 'G', 'T', '-'
    public int[][] genBlastScoreFunction() {
        int[][] scoreFunction = new int[5][5];
        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 5; j++) {
                if (i == 4 && j == 4) {
                    scoreFunction[i][j] = 0;
                } else if (i == 4 || j == 4) {
                    scoreFunction[i][j] = -8;
                } else if (i == j) {
                    scoreFunction[i][j] = 5;
                } else {
                    scoreFunction[i][j] = -4;
                }
            }
        }
        return scoreFunction;
    }

    public int[] calFinalScore() {
        int[] stats = new int[2];
        int score = 0;
        int matchCount = 0;
        for (int i = 0; i < a.length(); i++) {
            score += scoreFunction[getIndex(a.charAt(i))][getIndex(b.charAt(i))] +
                     scoreFunction[getIndex(a.charAt(i))][getIndex(c.charAt(i))] +
                     scoreFunction[getIndex(b.charAt(i))][getIndex(c.charAt(i))];
            if (a.charAt(i) == b.charAt(i) && a.charAt(i) == c.charAt(i)) {
                matchCount++;
            }
        }
        stats[0] = matchCount;
        stats[1] = score;
        return stats;
    }

    public void printMatrix(int[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                System.out.print(matrix[i][j] + " ");
            }
            System.out.println();
        }
    }
}
