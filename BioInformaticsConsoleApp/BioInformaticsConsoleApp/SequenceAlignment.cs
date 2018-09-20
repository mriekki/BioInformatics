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
            var M = new int[m + 1, n + 1];      // dynamic programming buttom up memory table
            var T = new string[m + 1, n + 1];   // trace back

            for (int i = 0; i < m + 1; i++)
                M[i, 0] = i * p;
            for (int j = 0; j < n + 1; j++)
                M[0, j] = j * p;

            T[0, 0] = DONE;
            for (int i = 1; i < m + 1; i++)
                T[i, 0] = UP;
            for (int j = 1; j < n + 1; j++)
                T[0, j] = LEFT;

            // calc
            for (int i = 1; i < m + 1; i++)
            {
                for (int j = 1; j < n + 1; j++)
                {
                    char vi = xs.ElementAt(i - 1);
                    char wj = ys.ElementAt(j - 1);

                    var alpha = Alpha(xs.ElementAt(i - 1).ToString(), ys.ElementAt(j - 1).ToString());

                    var diag = alpha + M[i - 1, j - 1];
                    var up = p + M[i - 1, j];
                    var left = p + M[i, j - 1];
                    var skip = ZERO;

                    var max = Max(diag, up, left, skip);
                    M[i, j] = max;

                    if (max > overallMaxScore)
                    {
                        overallMaxScore = max;
                        maxI = i;
                        maxJ = j;
                    }

                    if (max == ZERO)
                        T[i, j] = SKIP;
                    else if (max == up)
                        T[i, j] = UP;
                    else if (max == left)
                        T[i, j] = LEFT;
                    else
                        T[i, j] = DIAG;
                }
            }

            var traceBack = ParseTraceBack2(T, maxI + 1, maxJ + 1, true);

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

        Sequence SequenceAlign(string xs, string ys)
        {
            const int p = -5; //gap penalty, knowledge by looking at matrix file
            int m = xs.Length;
            int n = ys.Length;

            // init the matrix
            var M = new int[m + 1, n + 1]; // dynamic programming buttom up memory table
            var T = new string[m + 1, n + 1]; // trace back

            for (int i = 0; i < m + 1; i++)
                M[i, 0] = i * p;
            for (int j = 0; j < n + 1; j++)
                M[0, j] = j * p;

            T[0, 0] = DONE;
            for (int i = 1; i < m + 1; i++)
                T[i, 0] = UP;
            for (int j = 1; j < n + 1; j++)
                T[0, j] = LEFT;

            // calc
            for (int i = 1; i < m + 1; i++)
            {
                for (int j = 1; j < n + 1; j++)
                {
                    char vi = xs.ElementAt(i-1);
                    char wj = ys.ElementAt(j-1);

                    var alpha = Alpha(xs.ElementAt(i - 1).ToString(), ys.ElementAt(j - 1).ToString());

                    var diag = alpha + M[i - 1, j - 1];
                    var up = p + M[i - 1, j];
                    var left = p + M[i, j - 1];

                    var max = Max(diag, up, left);
                    M[i, j] = max;

                    if (max == up)
                        T[i, j] = UP;
                    else if (max == left)
                        T[i, j] = LEFT;
                    else
                        T[i, j] = DIAG; 
                }
            }

            var traceBack = ParseTraceBack(T, m + 1, n + 1);

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

            var sequence = new Sequence() { Score = M[m, n], Path = traceBack, One = vals[0], Two = vals[1] };
            return sequence;
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
