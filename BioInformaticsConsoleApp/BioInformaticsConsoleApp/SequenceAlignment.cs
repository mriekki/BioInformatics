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
            //linesMatrix = ReadFileToList("..\\..\\..\\Data Files\\BLOSUM62.txt");
            linesMatrix = ReadFileToList("..\\..\\..\\Data Files\\PAM250.txt");
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


        static string ReverseString(string s)
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
        Sequence SequenceAlign(string xs, string ys)
        {
            const int p = -5; //gap penalty, knowledge by looking at matrix file
            int m = xs.Length;
            int n = ys.Length;

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

        class Sequence
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
