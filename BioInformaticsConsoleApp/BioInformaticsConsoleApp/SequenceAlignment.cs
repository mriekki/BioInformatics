using System;
using System.Collections.Generic;
using System.Text;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;

namespace BioInformaticsConsoleApp
{
    public class SequenceAlignment
    {
        const int Size = 21; // BLOSUM62 matrix
        readonly private List<string> linesMatrix;
        readonly List<LifeForm> lifeForms = new List<LifeForm>();
        readonly private int[,] matrix = new int[Size, Size];
        readonly Dictionary<string, int> lookup = new Dictionary<string, int>();
        const int DefaultMismatchPenalty = -5; // if not defined, use this
        static readonly string NL = Environment.NewLine;

        // trace back
        const string DONE = @"¤";
        //const string DIAG = @"\";
        const string DIAG = "\u2196";
        //const string UP = @"^";
        const string UP = "\u2191";
        //const string LEFT = @"<";
        const string LEFT = "\u2190";
        const string SKIP = @"S";
        const int ZERO = 0;

        const int UPi = 1;
        const int LEFTi = 2;
        const int DIAGi = 3;

        // print alignment
        const string GAP = @"-";

        public SequenceAlignment()
        {
            linesMatrix = ReadFileToList("..\\..\\..\\Data Files\\BLOSUM62.txt");
            //linesMatrix = ReadFileToList("..\\..\\..\\Data Files\\PAM250.txt");
        }


        // BLOSUM62.txt
        void ParseMatrixFile()
        {
            const int startLine = 7; // # lines before line 7
            for (int i = 0; i < linesMatrix.Count; i++)
            {
                if (i >= startLine)
                {
                    var line = linesMatrix.ElementAt(i);
                    String[] row = Regex.Split(line, @"\s+"); // white space
                    lookup.Add(row[0], i - startLine);

                    for (int j = 1; j < row.Length; j++)
                    {
                        var s = row[j];
                        if (string.IsNullOrEmpty(s))
                            continue;
                        matrix[i - startLine, j - 1] = Int32.Parse(s);
                    }
                }
            }
        }

        static void PrintMatrix<T>(T[,] A, int I, int J)
        {
            for (int i = 0; i < I; i++)
            {
                for (int j = 0; j < J; j++)
                {
                    var v = A[i, j];
                    P(v + " ");
                }
                PL();
            }
        }


        public static string ReverseString(string s)
        {
            char[] arr = s.ToCharArray();
            Array.Reverse(arr);
            return new string(arr);
        }

        public void Run(string DNA1 = "", string DNA2 = "")
        {
            // put your data here
            lifeForms.Add(new LifeForm() { Name = "Sphinx", DNA = "KQRK" });
            lifeForms.Add(new LifeForm() { Name = "Bandersnatch", DNA = "KAK" });
            lifeForms.Add(new LifeForm() { Name = "Snark", DNA = "KQRIKAAKABK" });

            ParseMatrixFile();

            //var sequenceAlign = SequenceAlign(DNA1, DNA2);

            var sequenceAlign = LocalAlignment(DNA1, DNA2);


            /*            for (int i = 0; i < lifeForms.Count; i++)
                        {
                            for (int j = 0; j < i; j++)
                            {
                                if (i == j)
                                    continue;

                                var l1 = lifeForms.ElementAt(i);
                                var l2 = lifeForms.ElementAt(j);
                                var sequenceAlign = SequenceAlign(l1.DNA, l2.DNA).ToString();
                                Console.WriteLine(string.Format("Matching species: {0} and {1}, DNA: {2} and {3}\n{4}"
                                    , l1.Name, l2.Name, l1.DNA, l2.DNA, sequenceAlign));
                            }
                        }   */
        }

        public string[] multipleSequenceAlignment(string r, string s, string t)
        {
            int[,,] DpTable = new int[r.Length + 1, s.Length + 1, t.Length + 1];
            int[,,] backtrack = new int[r.Length + 1, s.Length + 1, t.Length + 1];

            for (int i = 1; i < r.Length + 1; i++)
            {
                for (int j = 1; j < s.Length + 1; j++)
                {
                    for (int k = 1; k < t.Length + 1; k++)
                    {
                        int opt1 = DpTable[i - 1, j, k];            // (r_i, -,   -)
                        int opt2 = DpTable[i, j - 1, k];            // (-,   s_j, -)
                        int opt3 = DpTable[i, j, k - 1];            // (-,   -,   t_k)
                        int opt4 = DpTable[i - 1, j - 1, k];        // (r_i, s_j, -)
                        int opt5 = DpTable[i - 1, j, k - 1];        // (r_i, -,   t_k)
                        int opt6 = DpTable[i, j - 1, k - 1];        // (-,   s_j, t_k)
                        int opt7 = DpTable[i - 1, j - 1, k - 1];    // (r_i, s_j, t_k)

                        if (r[i - 1] == s[j - 1] && s[j - 1] == t[k - 1])
                        {
                            opt7++;
                        }

                        DpTable[i, j, k] = opt1;
                        backtrack[i, j, k] = 1;

                        if (opt2 > DpTable[i, j, k])
                        {
                            DpTable[i, j, k] = opt2;
                            backtrack[i, j, k] = 2;
                        }
                        if (opt3 > DpTable[i, j, k])
                        {
                            DpTable[i, j, k] = opt3;
                            backtrack[i, j, k] = 3;
                        }
                        if (opt4 > DpTable[i, j, k])
                        {
                            DpTable[i, j, k] = opt4;
                            backtrack[i, j, k] = 4;
                        }
                        if (opt5 > DpTable[i, j, k])
                        {
                            DpTable[i, j, k] = opt5;
                            backtrack[i, j, k] = 5;
                        }
                        if (opt6 > DpTable[i, j, k])
                        {
                            DpTable[i, j, k] = opt6;
                            backtrack[i, j, k] = 6;
                        }
                        if (opt7 > DpTable[i, j, k])
                        {
                            DpTable[i, j, k] = opt7;
                            backtrack[i, j, k] = 7;
                        }
                    }
                }
            }
            String[] output = new String[4];
            int ii = r.Length;
            int jj = s.Length;
            int kk = t.Length;

            output[0] = "" + DpTable[ii, jj, kk];
            output[1] = "";
            output[2] = "";
            output[3] = "";
            int backT = 0;

            while(ii > 0 && jj > 0 && kk > 0) {
                backT = backtrack[ii, jj, kk];
                switch (backT) {
                    case 1:
                        output[1] = r[ii-- - 1] + output[1];
                        output[2] = '-' + output[2];
                        output[3] = '-' + output[3];
                            break;
                    case 2:
                        output[1] = '-' + output[1];
                        output[2] = s[jj-- - 1] + output[2];
                        output[3] = '-' + output[3];
                            break;
                    case 3:
                        output[1] = '-' + output[1];
                        output[2] = '-' + output[2];
                        output[3] = t[kk-- - 1] + output[3];
                            break;
                    case 4:
                        output[1] = r[ii-- - 1] + output[1];
                        output[2] = s[jj-- - 1] + output[2];
                        output[3] = '-' + output[3];
                            break;
                    case 5:
                        output[1] = r[ii-- - 1] + output[1];
                        output[2] = '-' + output[2];
                        output[3] = t[kk-- - 1] + output[3];
                            break;
                    case 6:
                        output[1] = '-' + output[1];
                        output[2] = s[jj-- - 1] + output[2];
                        output[3] = t[kk-- - 1] + output[3];
                            break;
                    case 7:
                        output[1] = r[ii-- - 1] + output[1];
                        output[2] = s[jj-- - 1] + output[2];
                        output[3] = t[kk-- - 1] + output[3];
                            break;
                }
            }

            while(ii > 0) {
                output[1] = r[ii-- - 1] + output[1];
                output[2] = '-' + output[2];
                output[3] = '-' + output[3];
            }

            while(jj > 0) {
                    output[1] = '-' + output[1];
                    output[2] = s[jj-- - 1] + output[2];
                    output[3] = '-' + output[3];
            }

            while(kk > 0) {
                    output[1] = '-' + output[1];
                    output[2] = '-' + output[2];
                    output[3] = t[kk-- - 1] + output[3];
            }

            return output;
    }
        public string FindMiddleEdge(string xs, string ys)
        {
            string result = "";

            const int p = -5; //gap penalty, knowledge by looking at matrix file
            int m = xs.Length;
            decimal nn = ys.Length;
            nn = nn / 2;
            decimal n = Math.Ceiling(nn);

            ParseMatrixFile();

            // init the matrix
            var DpTable = new int[m + 1, (int)n + 1]; // dynamic programming buttom up memory table

            for (int i = 0; i < m + 1; i++)
                DpTable[i, 0] = i * p;
            for (int j = 0; j < n + 1; j++)
                DpTable[0, j] = j * p;


            int MaxScore = -9999;
            int iEdge = 0;
            int iPrevEdge = 0;
            int jEdge = 0;
            int jPrevEdge = 0;

            // calc
            for (int i = 1; i < m + 1; i++)
            {
                for (int j = 1; j < n + 1; j++)
                {
                    char vi = xs.ElementAt(i - 1);
                    char wj = ys.ElementAt(j - 1);

                    var alpha = Alpha(xs.ElementAt(i - 1).ToString(), ys.ElementAt(j - 1).ToString());

                    var diag = alpha + DpTable[i - 1, j - 1];
                    var up = p + DpTable[i - 1, j];
                    var left = p + DpTable[i, j - 1];

                    var max = Max(diag, up, left);
                    DpTable[i, j] = max;

                    if (j == n)     // middle column
                    {
                        if (max > MaxScore)
                        {
                            MaxScore = max;
                            iEdge = i;
                            jEdge = j;

                            if (max == up)
                            {
                                iPrevEdge = iEdge - 1;
                                jPrevEdge = jEdge;
                            }
                            else if (max == left)
                            {
                                iPrevEdge = iEdge;
                                jPrevEdge = jEdge -1;
                            }
                            else
                            {
                                iPrevEdge = iEdge - 1;
                                jPrevEdge = jEdge - 1;
                            }

                        }
                    }
                }
            }

            result = "(" + iPrevEdge.ToString() + ", " + jPrevEdge.ToString() + ") " + "(" + iEdge.ToString() + ", " + jEdge.ToString() + ")";

            return result;
        }

        public List<string> AlignmentWithAffineGapPenalties(string s, string t)
        {
            List<string> result = new List<string>();
            int m = s.Length + 1;
            int n = t.Length + 1;
            int ge = 1;     // gap extenstion penalty
            int go = 11;    // gap opending penalty

            ParseMatrixFile();

            // Uses BLOSUM62. PASS go AND ge AS POSITIVE INTEGER, NOT NEGATIVE!
            int[,,] DpTable = new int[3, m, n];       // S[0] = lower, S[1] = middle, S[2] = upper
            int[,,] backtrack = new int[3, m, n];    // 1 = right, 2 = down, 3 = diag

            for (int i = 1; i < m; i++)
            {
                DpTable[0, i, 0] = -1 * (go + (i - 1) * ge);
                DpTable[1, i, 0] = -1 * (go + (i - 1) * ge);
                DpTable[2, i, 0] = -9999 * go;
                backtrack[0, i, 0] = 0;
                backtrack[1, i, 0] = 0;
                backtrack[2, i, 0] = 0;
            }

            for (int j = 1; j < n; j++)
            {
                DpTable[0, 0, j] = -9999 * go;
                DpTable[1, 0, j] = -1 * (go + (j - 1) * ge);
                DpTable[2, 0, j] = -1 * (go + (j - 1) * ge);
                backtrack[0, 0, j] = 1;
                backtrack[1, 0, j] = 1;
                backtrack[2, 0, j] = 1;
            }

            for (int i = 1; i < m; i++)
            {
                for (int j = 1; j < n; j++)
                {
                    int low1 = DpTable[0, i - 1, j] - ge;
                    int low2 = DpTable[1, i - 1, j] - go;
                    if (low1 > low2)
                    {
                        DpTable[0, i, j] = low1;
                        backtrack[0, i, j] = 0;
                    }
                    else
                    {
                        DpTable[0, i, j] = low2;
                        backtrack[0, i, j] = 1;
                    }

                    int up1 = DpTable[2, i, j - 1] - ge;
                    int up2 = DpTable[1, i, j - 1] - go;
                    if (up1 > up2)
                    {
                        DpTable[2, i, j] = up1;
                        backtrack[2, i, j] = 0;
                    }
                    else
                    {
                        DpTable[2, i, j] = up2;
                        backtrack[2, i, j] = 1;
                    }

                    var alpha = Alpha(s.ElementAt(i - 1).ToString(), t.ElementAt(j - 1).ToString());

                    int opt1 = DpTable[0, i, j];
                    int opt2 = DpTable[1, i - 1, j - 1] + alpha;
                    int opt3 = DpTable[2, i, j];

                    DpTable[1, i, j] = opt1;
                    backtrack[1, i, j] = 0;

                    if (opt2 > DpTable[1, i, j])
                    {
                        DpTable[1, i, j] = opt2;
                        backtrack[1, i, j] = 1;
                    }
                    if (opt3 > DpTable[1, i, j])
                    {
                        DpTable[1, i, j] = opt3;
                        backtrack[1, i, j] = 2;
                    }
                }
            }

            int index = s.Length;
            int jj = t.Length;
            int bestSIJ = 0;
            int best = DpTable[0, index, jj];

            if (DpTable[1, index, jj] > best)
            {
                best = DpTable[1, index, jj];
                bestSIJ = 1;
            }
            if (DpTable[2, index, jj] > best)
            {
                best = DpTable[2, index, jj];
                bestSIJ = 2;
            }

            String[] output = new String[3];
            output[0] = "" + best;
            output[1] = "";
            output[2] = "";

            while (index > 0 && jj > 0)
            {
                if (bestSIJ == 0)
                {
                    if (backtrack[0, index, jj] == 1)
                    {
                        bestSIJ = 1;
                    }
                    output[1] = s[index-- - 1] + output[1];
                    output[2] = '-' + output[2];
                }
                else if (bestSIJ == 1)
                {
                    if (backtrack[1, index, jj] == 0)
                    {
                        bestSIJ = 0;
                    }
                    else if (backtrack[1, index, jj] == 2)
                    {
                        bestSIJ = 2;
                    }
                    else
                    {
                        output[1] = s[index-- - 1] + output[1];
                        output[2] = t[jj-- - 1] + output[2];
                    }
                }
                else
                {
                    if (backtrack[2, index, jj] == 1)
                    {
                        bestSIJ = 1;
                    }
                    output[1] = '-' + output[1];
                    output[2] = t[jj-- - 1] + output[2];
                }
            }

            if (index > 0)
            {
                output[1] = s.Substring(0, index) + output[1];
                String add = "";
                for (int x = 0; x < index; x++)
                {
                    add += '-';
                }
                output[2] = add + output[2];
            }

            if (jj > 0)
            {
                output[2] = t.Substring(0, jj) + output[2];
                String add = "";
                for (int x = 0; x < jj; x++)
                {
                    add += '-';
                }
                output[1] = add + output[1];
            }

            foreach (string str in output)
                result.Add(str);

            return result;
        }


        Sequence LocalAlignment(string xs, string ys)
        {
            const int p = -5; //gap penalty, knowledge by looking at matrix file
            int m = xs.Length;
            int n = ys.Length;
            int overallMaxScore = -9999;
            int maxI = m;
            int maxJ = n;

            // init the matrix
            var DpTable = new int[m + 1, n + 1];      // dynamic programming buttom up memory tabl
            var backtrack = new string[m + 1, n + 1];   // trace back

            for (int i = 0; i < m + 1; i++)
                DpTable[i, 0] = i * p;
            for (int j = 0; j < n + 1; j++)
                DpTable[0, j] = j * p;

            backtrack[0, 0] = DONE;
            for (int i = 1; i < m + 1; i++)
                backtrack[i, 0] = UP;
            for (int j = 1; j < n + 1; j++)
                backtrack[0, j] = LEFT;

            // calc
            for (int i = 1; i < m + 1; i++)
            {
                for (int j = 1; j < n + 1; j++)
                {
                    char vi = xs.ElementAt(i - 1);
                    char wj = ys.ElementAt(j - 1);

                    var alpha = Alpha(xs.ElementAt(i - 1).ToString(), ys.ElementAt(j - 1).ToString());

                    var diag = alpha + DpTable[i - 1, j - 1];
                    var up = p + DpTable[i - 1, j];
                    var left = p + DpTable[i, j - 1];
                    var skip = ZERO;

                    var max = Max(diag, up, left, skip);
                    DpTable[i, j] = max;

                    if (max > overallMaxScore)
                    {
                        overallMaxScore = max;
                        maxI = i;
                        maxJ = j;
                    }

                    if (max == ZERO)
                        backtrack[i, j] = SKIP;
                    else if (max == up)
                        backtrack[i, j] = UP;
                    else if (max == left)
                        backtrack[i, j] = LEFT;
                    else
                        backtrack[i, j] = DIAG;
                }
            }

            var traceBack = ParseTraceBack2(backtrack, maxI + 1, maxJ + 1, true);

            string[] vals = BuildAlignment2(traceBack, xs, ys, maxI, maxJ);

            Console.OutputEncoding = System.Text.Encoding.UTF8;
            //            PL("\nScore table");
            //            PrintMatrix(M, m + 1, n + 1);
            //            PL("\nTraceBack");
            //            PrintMatrix(T, m + 1, n + 1);
            //            PL();

            var sequence = new Sequence() { Score = overallMaxScore, Path = traceBack, One = vals[0], Two = vals[1] };
            return sequence;
        }

        public int EditDistance(string s, string t)
        {
            int distance = 0;
            int m = s.Length + 1;
            int n = t.Length + 1;
            int[,] DpTable = new int[m, n];
            int substituationCost = 0;

            for (int i = 0; i < m; i++)
                for (int j = 0; j < n; j++)
                    DpTable[i, j] = 0;

            for (int i = 0; i < m; i++)
                DpTable[i, 0] = i;

            for (int j = 0; j < n; j++)
                DpTable[0, j] = j;

            for (int j = 1; j < n; j++)
            {
                for (int i = 1; i < m; i++)
                {
                    if (s[i-1] == t[j-1])
                        substituationCost = 0;
                    else
                        substituationCost = 1;

                    DpTable[i, j] = Min(DpTable[i - 1, j] + 1,
                                        DpTable[i, j - 1] + 1,
                                        DpTable[i - 1, j - 1] + substituationCost);
                }
            }

            distance = DpTable[m-1, n-1];

            return distance;
        }
        public Sequence SequenceAlign(string xs, string ys)
        {
            const int p = -5; //gap penalty, knowledge by looking at matrix file
            int m = xs.Length;
            int n = ys.Length;

            ParseMatrixFile();

            // init the matrix
            var DpTable = new int[m + 1, n + 1]; // dynamic programming buttom up memory table
            var backtrack = new string[m + 1, n + 1]; // trace back

            for (int i = 0; i < m + 1; i++)
                DpTable[i, 0] = i * p;
            for (int j = 0; j < n + 1; j++)
                DpTable[0, j] = j * p;

            backtrack[0, 0] = DONE;
            for (int i = 1; i < m + 1; i++)
                backtrack[i, 0] = UP;
            for (int j = 1; j < n + 1; j++)
                backtrack[0, j] = LEFT;

            // calc
            for (int i = 1; i < m + 1; i++)
            {
                for (int j = 1; j < n + 1; j++)
                {
                    char vi = xs.ElementAt(i-1);
                    char wj = ys.ElementAt(j-1);

                    var alpha = Alpha(xs.ElementAt(i - 1).ToString(), ys.ElementAt(j - 1).ToString());

                    var diag = alpha + DpTable[i - 1, j - 1];
                    var up = p + DpTable[i - 1, j];
                    var left = p + DpTable[i, j - 1];

                    var max = Max(diag, up, left);
                    DpTable[i, j] = max;

                    if (max == up)
                        backtrack[i, j] = UP;
                    else if (max == left)
                        backtrack[i, j] = LEFT;
                    else
                        backtrack[i, j] = DIAG; 
                }
            }

            var traceBack = ParseTraceBack(backtrack, m + 1, n + 1);

            string[] vals = BuildAlignment(traceBack, xs, ys);

/*            var sb = new StringBuilder();
            string first, second;

            if (xs.Length != ys.Length)
            {
                string s;
                if (xs.Length > ys.Length)
                {
                    s = ys;
                    first = xs;
                }
                else
                {
                    s = xs;
                    first = ys;
                }

                int i = 0;
                foreach (var trace in traceBack)
                {
                    if (trace.ToString() == DIAG)
                        sb.Append(s.ElementAt(i++).ToString());
                    else
                        sb.Append(GAP);
                }

                second = sb.ToString();
            }
            else
            {
                first = xs;
                second = ys;
            }       */

            Console.OutputEncoding = System.Text.Encoding.UTF8;
//            PL("\nScore table");
//            PrintMatrix(M, m + 1, n + 1);
//            PL("\nTraceBack");
//            PrintMatrix(T, m + 1, n + 1);
//            PL();

            var sequence = new Sequence() { Score = DpTable[m, n], Path = traceBack, One = vals[0], Two = vals[1] };
            return sequence;
        }

        public List<string> OverlapAlignment(string s, string t)
        {
            List<string> result = new List<string>();
            int m = s.Length + 1;
            int n = t.Length + 1;

            int[,] DpTable = new int[m, n];
            int[,] backtrack = new int[m, n]; // 1 = right, 2 = down, 3 = diag, 4 = STOP
            int bi = -1;
            int bj = -1;
            int max = -9999;
            int indel = 2;

            for (int i = 1; i < m; i++)
            {
                for (int j = 1; j < n; j++)
                {
                    int opt1 = DpTable[i, j - 1] - indel; // right
                    int opt2 = DpTable[i - 1, j] - indel; // down
                    int opt3 = DpTable[i - 1, j - 1]; // diag

                    if (s[i - 1] == t[j - 1])
                    {
                        opt3++;
                    }
                    else
                    {
                        opt3 -= 2;
                    }

                    DpTable[i, j] = opt1;
                    backtrack[i,j] = 1;

                    if (opt2 > DpTable[i, j])
                    {
                        DpTable[i, j] = opt2;
                        backtrack[i, j] = 2;
                    }
                    if (opt3 > DpTable[i, j])
                    {
                        DpTable[i, j] = opt3;
                        backtrack[i, j] = 3;
                    }

                    if (i == s.Length || j == t.Length)
                    {
                        if (DpTable[i, j] > max)
                        {
                            max = DpTable[i, j];
                            bi = i;
                            bj = j;
                        }
                    }
                }
            }

            string[] output = new string[3];
            output[1] = "";
            output[2] = "";
            output[0] = "" + max;

            while (bi > 0 && bj > 0) 
            {
                if(backtrack[bi, bj] == 1)  // right
                { 
                    output[1] = '-' + output[1];
                    output[2] = t[bj-- - 1] + output[2];
                }
                else if(backtrack[bi, bj] == 2)   // down
                { 
                    output[1] = s[bi-- - 1] + output[1];
                    output[2] = '-' + output[2];
                }
                else if(backtrack[bi, bj] == 3)   // diag
                { 
                    output[1] = s[bi-- - 1] + output[1];
                    output[2] = t[bj-- - 1] + output[2];
                }
            }   
                                
            foreach (string str in output)
                result.Add(str);

            return result;
        }

        public List<string> FittingAlignment(string s1, string s2)
        {
            List<string> result = new List<string>();
            string s = "";
            string t = "";
            
            if (s1.Length > s2.Length)
            {
                s = s1;
                t = s2;
            }
            else
            {
                s = s2;
                t = s1;
            }
            int m = s.Length + 1;
            int n = t.Length + 1;

            int[,] DpTable = new int[m, n];
            int[,] backtrack = new int[m, n];


            for (int i = 1; i < m; i++)
            {
                for (int j = 1; j < n; j++)
                {
                    int opt1 = DpTable[i - 1, j] - 1;    // up
                    int opt2 = DpTable[i, j - 1] - 1;    // left
                    int opt3 = DpTable[i - 1, j - 1];    // diag

                    if (s[i - 1] == t[j - 1])
                        opt3++;
                    else
                        opt3--;

                    DpTable[i, j] = opt1;
                    backtrack[i, j] = UPi;

                    if (opt2 > DpTable[i, j])
                    {
                        DpTable[i, j] = opt2;
                        backtrack[i, j] = LEFTi;
                    }

                    if (opt3 > DpTable[i, j])
                    {
                        DpTable[i, j] = opt3;
                        backtrack[i, j] = DIAGi;
                    }
                }
            }

            int jj = t.Length;
            int max = -9999;
            int index = -1;

            for (int x = t.Length; x < m; x++)
            {
                if (DpTable[x,jj] > max)
                {
                    max = DpTable[x, jj];
                    index = x;
                }
            }

            string[] output = new string[3];

            output[0] = "" + max;
            output[1] = "";
            output[2] = "";

            while(jj > 0) 
            {
                if(backtrack[index, jj] == UPi) 
                {
                    output[1] = s[index-- - 1] + output[1];
                    output[2] = '-' + output[2];
                }
                else if(backtrack[index, jj] == LEFTi) 
                {
                    output[1] = '-' + output[1];
                    output[2] = t[jj-- - 1] + output[2];
                }
                else if(backtrack[index, jj] == DIAGi) 
                {
                    output[1] = s[index-- - 1] + output[1];
                    output[2] = t[jj-- - 1] + output[2];
                }
            }

            foreach (string str in output)
                result.Add(str);

            return result;
        }

        static string[] BuildAlignment(string traceback, string v, string w)
        {
            string[] result = new string[2];
            string str1 = "";
            string str2 = "";

            int i1 = 0;
            int i2 = 0;

            foreach (char c in traceback)
            {
                if (c.ToString() == DIAG)
                {
                    str1 += v.ElementAt(i1++).ToString();
                    str2 += w.ElementAt(i2++).ToString();
                }
                else if (c.ToString() == LEFT)
                {
                    str1 += GAP.ToString();
                    str2 += w.ElementAt(i2++).ToString();
                }
                else if (c.ToString() == UP)
                {
                    str1 += v.ElementAt(i1++);
                    str2 += GAP.ToString();
                }
            }

            result[0] = str1;
            result[1] = str2;

            return result;

        }

        static string[] BuildAlignment2(string traceback, string v, string w, int i, int j)
        {
            string[] result = new string[2];
            string str1 = "";
            string str2 = "";

            int i1 = i-1;
            int i2 = j-1;

            for (int s = traceback.Length - 1; s >= 0; s--)
            {
                char c = traceback[s];

                if (c.ToString() == SKIP)
                {
                    i1--;
                    i2--;
                }
                else if (c.ToString() == DIAG)
                {
                    str1 += v.ElementAt(i1--).ToString();
                    str2 += w.ElementAt(i2--).ToString();
                }
                else if (c.ToString() == LEFT)
                {
                    str1 += GAP.ToString();
                    str2 += w.ElementAt(i2--).ToString();
                }
                else if (c.ToString() == UP)
                {
                    str1 += v.ElementAt(i1--);
                    str2 += GAP.ToString();
                }
            }

            result[0] = ReverseString(str1);
            result[1] = ReverseString(str2);

            return result;

        }

        static string ParseTraceBack(string[,] T, int I, int J)
        {
            var sb = new StringBuilder();
            int i = I - 1;
            int j = J - 1;
            string path = T[i, j];

            while (path != DONE)
            {
                sb.Append(path);

                if (path == DIAG)
                {
                    i--;
                    j--;
                }
                else if (path == UP)
                    i--;
                else if (path == LEFT)
                    j--;
                else
                {
                    int kk = 0;
                }

                path = T[i, j];
            }

            return ReverseString(sb.ToString());
        }

        static string ParseTraceBack2(string[,] T, int I, int J, bool localAlign=false)
        {
            var sb = new StringBuilder();
            int i = I - 1;
            int j = J - 1;
            string path = T[i, j];

            while (path != DONE)
            {
                sb.Append(path);

                if (!localAlign)
                {
                    if (path == DIAG)
                    {
                        i--;
                        j--;
                    }
                    else if (path == UP)
                        i--;
                    else if (path == LEFT)
                        j--;
                    else
                    {
                        int kk = 0;
                    }
                }
                else
                {
                    if (path == SKIP)   // finished, exit loop
                    {
                        break;
                    }
                    else if (path == UP)
                        i--;
                    else if (path == LEFT)
                    {
                        j--;
                    }
                    else
                    {
                        i--;
                        j--;
                    }
                }

                path = T[i, j];
            }

            return ReverseString(sb.ToString());
        }


        static int Max(params int[] numbers)
        {
            return numbers.Max();
        }

        static int Min(params int[] numbers)
        {
            return numbers.Min();
        }

        int Alpha(string x, string y)
        {
            if (lookup.ContainsKey(x) && lookup.ContainsKey(y))
            {
                var i = lookup[x];
                var j = lookup[y];
                return matrix[i, j];
            }
            else if (x == y)
                return 1; // matrix file match * with *

            return DefaultMismatchPenalty; // default mismatch penalty
        }

        class LifeForm
        {
            public string Name { get; set; }
            public string DNA { get; set; }
        }

        public class Sequence
        {
            public int Score { get; set; }
            public string Path { get; set; }
            public string One { get; set; }
            public string Two { get; set; }
            public new string ToString()
            {
                var s = string.Format("score = {0}{1}one = {2}{3}two = {4}\n\n", Score, NL, One, NL, Two);
                return s;
            }
        }

        public static List<string> ReadFileToList(string path)
        {
            var lines = new List<string>();
            StreamReader streamReader = null;
            try
            {
                streamReader = File.OpenText(path);
                string line = streamReader.ReadLine();
                while (line != null)
                {
                    lines.Add(line);
                    line = streamReader.ReadLine();
                }
            }
            catch (FileNotFoundException ex)
            {
                PL(string.Format("Error cannot find file {0}", path));
                PL(ex.Message + Environment.NewLine + ex.StackTrace);
            }
            catch (Exception ex)
            {
                PL(string.Format("Error reading file {0}", path));
                PL(ex.Message + Environment.NewLine + ex.StackTrace);
            }
            finally
            {
                if (streamReader != null)
                    streamReader.Close();
            }
            return lines;
        }

        public static void PL(object o) { Console.WriteLine(o); } //alias
        public static void PL() { Console.WriteLine(); } //alias
        public static void P(object o) { Console.Write(o); } //alias
    }
}
