using System;
using System.Collections.Generic;
using System.Text;
using BioInformaticsConsoleApp;

namespace BioInformaticsConsoleApp
{
    class HirschbergAlgorithm
    {
        private string LCS;
        private enum DATA_TYPE { ALPHA, NUMERIC};

        public int[] algorithmB(int mlength, int nlength, String a, String b)
        {
//            int[][] c = new int[2][nlength + 1];
            int[][] c = new int[2][];
//            c[0] = new int[2];
            c[0] = new int[nlength + 1];
            c[1] = new int[nlength + 1];

            for (int j = 0; j < nlength; j++)
            {
                c[1][j] = 0;
            }

            for (int i = 1; i <= mlength; i++)
            {

                for (int j = 0; j <= nlength; j++)
                {
                    c[0][j] = c[1][j];
                }


                for (int j = 1; j <= nlength; j++)
                {
                    if (a[i - 1] == b[j - 1])
                    {
                        c[1][j] = c[0][j - 1] + 1;
                    }
                    else
                    {
                        c[1][j] = Program.Max(c[1][j - 1], c[0][j]);
                    }
                }
            }

            return c[1];
        }


        public String LinearSpaceAlignment(int mlength, int nlength, string A, string B)
        {
            int ipos = 0;
            int j = 0;
            String c = "";


            if (nlength == 0)
            {   // For trivial problems initialize to empty string if the length is zero.
                c = "";
            }
            else if (mlength == 1)
            {  // If the length is one
                c = "";
                for (j = 0; j < nlength; j++)
                { // Iterate through n
                    if (A[0] == B[j])
                    { //if A[i] = B[j]
                        c = "" + A[0]; //Initialize c with char at the zeroth position
                        break;
                    }
                }
            }
            else
            { // else split the problems
                ipos = (int)Math.Floor(((double)mlength) / 2);

                // Evaluation of L1 and L2 termed as length 1 and length 2
                int[] length1 = algorithmB(ipos, nlength, A.Substring(0, ipos), B);
                int[] length2 = algorithmB(mlength - ipos, nlength, SequenceAlignment.ReverseString(A.Substring(ipos)), SequenceAlignment.ReverseString(B));
//                int[] length2 = algorithmB(mlength - ipos, nlength, reverse_String(A.substring(ipos)), reverse_String(B));

                int k = findC(length1, length2, nlength);

                // Find the best match through recursive calls                        
                String c1 = LinearSpaceAlignment(ipos, k, A.Substring(0, ipos), B.Substring(0, k));
                String c2 = LinearSpaceAlignment(mlength - ipos, nlength - k, A.Substring(ipos), B.Substring(k));

                c = c1 + c2;// Concatenation of two strings
            }

            return c; // return the LCS
        }

        private string reverse_String(object p)
        {
            throw new NotImplementedException();
        }

        public int findC(int[] length1, int[] length2, int n)
        {
            int m = 0;
            int c = 0;

            for (int j = 0; j <= n; j++)
            {
                if (m < (length1[j] + length2[n - j]))
                {
                    m = length1[j] + length2[n - j];
                    c = j;
                }
            }

            return c;
        }


        /*
         # GLOBAL ALIGNMENT IN LINEAR SPACE PROBLEM

# Find the highest-scoring alignment between two strings using a scoring matrix in linear space.

# Given:  Two long amino acid strings (of length approximately 10,000).
# Return: The maximum alignment score of these strings, followed by an alignment achieving this maximum score. Use the BLOSUM62 scoring matrix and indel penalty = 5.

import sys

class BLOSUM62( object ):
    def __init__(self):
        with open('BLOSUM62.txt') as input:
            items = [line.strip().split() for line in input.readlines()]
            self.m = {(item[0], item[1]):int(item[2]) for item in items}

    def __getitem__(self, pair):
        return self.m[pair[0], pair[1]]

def linearSpaceGlobalAlignment( v, w, matrix, sigma ):

    def middleColumnScore( v, w, matrix, sigma ):    
        S = [[i*j*sigma for j in xrange(-1, 1)] for i in xrange(len(v)+1)]
        S[0][1] = -sigma
        backtrack = [0]*(len(v)+1)
        for j in xrange(1, len(w)/2+1):
            for i in xrange(0, len(v)+1):
                if i == 0:
                    S[i][1] = -j*sigma
                else:
                    scores = [S[i-1][0] + matrix[v[i-1], w[j-1]], S[i][0] - sigma, S[i-1][1] - sigma]
                    S[i][1] = max(scores)
                    backtrack[i] = scores.index(S[i][1])
            if j != len(w)/2:
                S = [[row[1]]*2 for row in S]
        return [row[1] for row in S], backtrack
    
    def middleEdge( v, w, matrix, sigma ):
        sourceToMiddle = middleColumnScore(v, w, matrix, sigma)[0]
        middleToSink, backtrack = map(lambda l: l[::-1], middleColumnScore(v[::-1], w[::-1]+['', '$'][len(w) % 2 == 1 and len(w) > 1], matrix, sigma))
        scores = map(sum, zip(sourceToMiddle, middleToSink))
        maxMiddle = max(xrange(len(scores)), key=lambda i: scores[i])
        if maxMiddle == len(scores) - 1:
            nextNode = (maxMiddle, len(w)/2 + 1)
        else:
            nextNode = [(maxMiddle + 1, len(w)/2 + 1), (maxMiddle, len(w)/2 + 1), (maxMiddle + 1, len(w)/2),][backtrack[maxMiddle]]
        return (maxMiddle, len(w)/2), nextNode

    def globalAlignment( v, w, matrix, sigma ):
        S = [[0 for repeat_j in xrange(len(w)+1)] for repeat_i in xrange(len(v)+1)]
        backtrack = [[0 for repeat_j in xrange(len(w)+1)] for repeat_i in xrange(len(v)+1)]
        for i in xrange(1, len(v)+1):
            S[i][0] = -i*sigma
        for j in xrange(1, len(w)+1):
            S[0][j] = -j*sigma
        for i in xrange(1, len(v)+1):
            for j in xrange(1, len(w)+1):
                scores = [S[i-1][j] - sigma, S[i][j-1] - sigma, S[i-1][j-1] + matrix[v[i-1], w[j-1]]]
                S[i][j] = max(scores)
                backtrack[i][j] = scores.index(S[i][j])
        insert_indel = lambda word, i: word[:i] + '-' + word[i:]
        vAligned, wAligned = v, w
        i, j = len(v), len(w)
        maxScore = str(S[i][j])
        while i*j != 0:
            if backtrack[i][j] == 0:
                i -= 1
                wAligned = insert_indel(wAligned, j)
            elif backtrack[i][j] == 1:
                j -= 1
                vAligned = insert_indel(vAligned, i)
            else:
                i -= 1
                j -= 1
        for repeat in xrange(i):
            wAligned = insert_indel(wAligned, 0)
        for repeat in xrange(j):
            vAligned = insert_indel(vAligned, 0)
        return maxScore, vAligned, wAligned

    def linearSpaceAlignment( top, bottom, left, right ):
        if left == right:
            return [v[top:bottom], '-'*(bottom - top)]
        elif top == bottom:
            return ['-'*(right - left), w[left:right]]
        elif bottom - top == 1 or right - left == 1:
            return globalAlignment(v[top:bottom], w[left:right], matrix, sigma)[1:]
        else:
            midNode, nextNode = middleEdge(v[top:bottom], w[left:right], matrix, sigma)
            midNode = tuple(map(sum, zip(midNode, [top, left])))
            nextNode = tuple(map(sum, zip(nextNode, [top, left])))
            current = [['-', v[midNode[0] % len(v)]][nextNode[0] - midNode[0]], ['-', w[midNode[1] % len(w)]][nextNode[1] - midNode[1]]]
            A = linearSpaceAlignment(top, midNode[0], left, midNode[1])
            B = linearSpaceAlignment(nextNode[0], bottom, nextNode[1], right)
            return [A[i] + current[i] + B[i] for i in xrange(2)]
    vAligned, wAligned = linearSpaceAlignment(0, len(v), 0, len(w))
    score = sum([-sigma if '-' in pair else matrix[pair] for pair in zip(vAligned, wAligned)])
    return str(score), vAligned, wAligned

if __name__ == '__main__':
    with open(sys.argv[1]) as input:
        word1, word2 = [line.strip() for line in input.readlines()]
    alignment = linearSpaceGlobalAlignment(word1, word2, BLOSUM62(), 5)
    print '\n'.join(alignment)
    */
    }
}
