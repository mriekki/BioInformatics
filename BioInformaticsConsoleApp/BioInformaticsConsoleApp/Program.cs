using System;
using System.IO;
using System.Collections.Generic;
using System.Text;
using System.Collections;
using System.Linq;

namespace BioInformaticsConsoleApp
{
    class Program
    {
        private const int NucleotideSize = 4;
        private const int INF = -9999;
        private static char[] zero_ones = { '0', '1' };
        private static string right_arrow = "\u2192";
        private static string down_arrow = "\u2193";
        private static string diag_arrow = "\u2198";
        private static string blank = ".";

        private static Dictionary<string, string> condonTable = new Dictionary<string, string> {
        {"UUU","F"}, {"UUC","F"}, {"UUA","L"}, {"UUG","L"},
        {"UCU","S"}, {"UCC","S" }, {"UCA","S" }, {"UCG","S" },
        { "UAU","Y" }, {"UAC","Y" }, {"UAA","STOP" }, {"UAG","STOP" },
        { "UGU","C" }, {"UGC","C" }, {"UGA","STOP" }, {"UGG","W" },
        { "CUU","L" }, {"CUC","L" }, {"CUA","L" }, {"CUG","L" },
        { "CCU","P" }, {"CCC","P" }, {"CCA","P" }, {"CCG","P" },
        { "CAU","H" }, {"CAC","H" }, {"CAA","Q" }, {"CAG","Q" },
        { "CGU","R" }, {"CGC","R" }, {"CGA","R" }, {"CGG","R" },
        { "AUU","I" }, {"AUC","I" }, {"AUA","I" }, {"AUG","M" },
        { "ACU","T" }, {"ACC","T" }, {"ACA","T" }, {"ACG","T" },
        { "AAU","N" }, {"AAC","N" }, {"AAA","K" }, {"AAG","K" },
        { "AGU","S" }, {"AGC","S" }, {"AGA","R" }, {"AGG","R" },
        { "GUU","V" }, {"GUC","V" }, {"GUA","V" }, {"GUG","V" },
        { "GCU","A" }, {"GCC","A" }, {"GCA","A" }, {"GCG","A" },
        { "GAU","D" }, {"GAC","D" }, {"GAA","E" }, {"GAG","E" },
        { "GGU","G" }, {"GGC","G" }, {"GGA","G" }, {"GGG","G"} };

        private static Dictionary<string, int> integerMass = new Dictionary<string, int>
        {
            { "G", 57 }, {"A", 71 }, {"S", 87 }, {"P", 97 },
            { "V", 99 }, {"T", 101 }, {"C", 103 }, {"I", 113 },
            {"L", 113 }, {"N", 114 }, {"D", 115 }, {"K", 128 },
            {"Q", 128 }, {"E", 129 }, {"M", 131 }, {"H", 137 },
            {"F", 147 }, {"R", 156 }, {"Y", 163 }, {"W", 186 } };


        private static List<Int64> peptideMass = new List<Int64>
        { 
              { 57 }, { 71 }, {87 }, { 97 }, 
              { 99 }, {101 }, {103 }, {113 }, 
              {114 }, {115 }, {128 }, 
              {129 }, {131 }, {137 }, 
              {147 }, {156 }, {163 }, {186 } };

        //private const string inputFile = "..\\..\\..\\Data Files\\dataset_289_5.txt";
        private const string inputFile = "..\\..\\..\\Data Files\\MyData.txt";
        private const string method = "ColoredEdges";

        public static List<string> sharedKmers(int k, string x, string y)
        {
            Dictionary<string, List<int>> xKmers = new Dictionary<string, List<int>>();
            for (int i = 0; i < x.Length - k + 1; i++)
            {
                string sub = x.Substring(i, k);
                string rev = reverseComplement(sub);
                List<int> entry = new List<int>();

                if (!xKmers.ContainsKey(sub))
                {
                    List<int> temp = new List<int>();
                    temp.Add(i);
                    xKmers.Add(sub, temp);
                }
                else
                {
                    entry = xKmers[sub];
                    entry.Add(i);
                }

                if (!xKmers.ContainsKey(rev))
                {
                    List<int> temp = new List<int>();
                    temp.Add(i);
                    xKmers.Add(rev, temp);
                }
                else
                {
                    entry = xKmers[rev];
                    entry.Add(i);
                }
            }

            List<int[]> pairs = new List<int[]>();
            for (int i = 0; i < y.Length - k + 1; i++)
            {
                string sub = y.Substring(i, k);

                List<int> xInds = new List<int>();
                if (xKmers.ContainsKey(sub))
                {
                    xInds = xKmers[sub];
                    for (int j = 0; j < xInds.Count(); j++)
                    {
                        int[] pair = new int[2];
                        int xInd = xInds[j];
                        pair[0] = xInd;
                        pair[1] = i;
                        pairs.Add(pair);
                    }
                }
            }
            xKmers = null;

            pairs.Sort(arrayComparator);
            string val = "";

            List<string> output = new List<string>();

            for (int i = 0; i < pairs.Count(); ++i) 
            {
                val = "(" + pairs[i][0] + ", " + pairs[i][1] + ")";
                output.Add(val);
            }   

            return output;
        }



    public static string reverseComplement(string dna)
    {
        string output = "";
        for (int i = dna.Length - 1; i >= 0; i--)
        {
            char curr = dna[i];
            if (curr == 'A')
                    output += 'T';
                else if (curr == 'T')
                    output += 'A';
                else if (curr == 'C')
                    output += 'G';
                else if (curr == 'G')
                    output += 'C';
        }
        return output;
    }

public static List<Tuple<int, int>> MergeEdgeLists(List<Tuple<int, int>> pList, List<Tuple<int, int>> qList)
        {
            List<Tuple<int, int>> combinedEdgeList = new List<Tuple<int, int>>();

            bool onP = true;

            while (pList.Count > 0 || qList.Count > 0)
            {
                if (onP)
                {
                    combinedEdgeList.Add(Tuple.Create(pList[0].Item1, pList[0].Item2));
                    pList.RemoveAt(0);

                    onP = false;
                }
                else
                {
                    combinedEdgeList.Add(Tuple.Create(qList[0].Item1, qList[0].Item2));
                    qList.RemoveAt(0);
                    onP = true;
                }
            }

            return combinedEdgeList;
        }

        public static List<Tuple<int, int>> CreateEdgeList(string P)
        {
            List<Tuple<int, int>> edgeList = new List<Tuple<int, int>>();

            string[] pEdge = P.Split("), ");

            foreach (string str in pEdge)
            {
                string val = str.Substring(1, str.Length - 1);

                if (val[val.Length - 1] == ')')
                    val = val.Substring(0, val.Length - 1);

                string[] items = val.Split(", ");
                edgeList.Add(Tuple.Create(Int32.Parse(items[0]), Int32.Parse(items[1])));
            }
            return edgeList;
        }

        public static string TwoBreakOnGenome(string P, int i1, int i2, int i3, int i4)
        {
            string result = "";
            string g = ColoredEdges(P);

            g = TwoBreakOnGenomeGraph(g, i1, i2, i3, i4);

            result = GraphToGenome(g);

            return result;
        }

        public static string TwoBreakOnGenomeGraph(string GenomeGraph, int i1, int i2, int i3, int i4)
        {
            string result = "";
            List<Tuple<int, int>> breakpointGraph = new List<Tuple<int, int>>();
            List<Tuple<int, int>> coloredList = new List<Tuple<int, int>>();
            List<Tuple<int, int>> indicesList = new List<Tuple<int, int>>();

            string[] cycle = GenomeGraph.Split("), ");
            foreach (string str in cycle)
            {
                string val = str.Substring(1, str.Length - 1);

                if (val[val.Length - 1] == ')')
                    val = val.Substring(0, val.Length - 1);

                string[] items = val.Split(", ");
                coloredList.Add(Tuple.Create(Int32.Parse(items[0]), Int32.Parse(items[1])));
            }

            indicesList.Add(Tuple.Create(i1, i2));
            indicesList.Add(Tuple.Create(i2, i1));
            indicesList.Add(Tuple.Create(i3, i4));
            indicesList.Add(Tuple.Create(i4, i3));

            foreach (var t in coloredList)
            {
                bool found = false;

                foreach (var x in indicesList)
                {
                    if (EdgesEqual(t, x))
                    {
                        found = true;
                        break;
                    }
                }
                if (!found)
                    breakpointGraph.Add(t);
            }
            breakpointGraph.Add(Tuple.Create(i1, i3));
            breakpointGraph.Add(Tuple.Create(i2, i4));

            foreach (var s in breakpointGraph)
            {
                result += "(" + s.Item1 + ", " + s.Item2 + "), ";
            }

            result = result.Substring(0, result.Length - 2);

            return result;
        }

        public static bool EdgesEqual(Tuple<int, int> t1, Tuple<int, int> t2)
        {
            bool equal = false;

            if (t1.Item1 == t2.Item1 && t1.Item2 == t2.Item2)
                equal = true;

            return equal;
        }
        public static int GetCycles(List<Tuple<int, int>> bp)
        {
            bool[] visited = new bool[bp.Count()];
            int c = 0;

            for (int i = 1; i < bp.Count(); i++)
            {
                if (!visited[i])
                {
                    int curr = i;
                    while (!visited[curr])
                    {
                        visited[curr] = true;
                        int val = bp[curr].Item1;
                        int val2 = bp[curr].Item2;

                        if (!visited[bp[curr].Item1])
                        {
                            curr = bp[curr].Item1;
                        }
                        else
                        {
                            curr = bp[curr].Item2;
                        }   
                    }
                    c++;
                }
            }
            return c;
        }

        public static List<Tuple<int, int>> getCycles(List<Tuple<int, int>> edges)
        {
            List<Tuple<int, int>> cycles = new List<Tuple<int, int>>();

            List<Tuple<int, int>> edgesPCopy = new List<Tuple<int, int>>(edges);

            bool isNewCycle = true;

            int first, second;
            int firstElementOfCycle = -1;
            int lastNode = -1;

            while (edgesPCopy.Count() > 0)
            {
                if (isNewCycle)
                {
                    first = edgesPCopy[0].Item1;
                    second = edgesPCopy[0].Item2;
                    cycles.Add(Tuple.Create(first, second));

                    edgesPCopy.RemoveAt(0);

                    firstElementOfCycle = first;
                    if (second % 2 == 0)
                    {
                        lastNode = second - 1;
                    }
                    else
                    {
                        lastNode = second + 1;
                    }
                    isNewCycle = false;
                }
                else
                {
                    //indexList = new ArrayList<>();
                    List<Tuple<int, int>> indexList = new List<Tuple<int, int>>();

                    int index = 0;
                    int dd = 0;
                    foreach (var v in edgesPCopy)
                    {
                        indexList.Add(v);
                        dd++;
                        if (v.Item2 == lastNode)
                            index = dd - 1;
                    }

                    first = edgesPCopy[index].Item1;
                    second = edgesPCopy[index].Item2;
                    edgesPCopy.RemoveAt(index);

                    if (lastNode == first)
                    {
                        cycles.Add(Tuple.Create(first, second));

                        if (second % 2 == 0)
                        {
                            lastNode = second - 1;
                        }
                        else
                        {
                            lastNode = second + 1;
                        }
                    }
                    else
                    {
                        cycles.Add(Tuple.Create(second, first));

                        if (first % 2 == 0)
                        {
                            lastNode = first - 1;
                        }
                        else
                        {
                            lastNode = first + 1;
                        }
                    }
                    if (lastNode == firstElementOfCycle)
                    {
                        //cycles = Rotate(cycles, 1);
                        int tmpFirst = cycles[cycles.Count() - 1].Item1;
                        int tmpSecond = cycles[cycles.Count() - 1].Item2;
                        cycles[cycles.Count - 1] = Tuple.Create(tmpSecond, tmpFirst);

                        isNewCycle = true;

                    }
                }
            }
            return cycles;
        }

        public static List<Tuple<int,int>> Rotate(List<Tuple<int,int>> list, int offset)
        {
            return list.Skip(offset).Concat(list.Take(offset)).ToList();
        }

        public static List<Tuple<int, int>> breakpointGraph(string g1, string g2)
        {
            string[] g1chroms = g1.Substring(1, g1.Length - 2).Split(")(");
            List<Tuple<int, int>> breakpointGraph = new List<Tuple<int, int>>();

            int[][] g1chromsInt = new int[g1chroms.Length][];

            int n = 0;
            for (int i = 0; i < g1chroms.Length; i++)
            {
                string[] parts = g1chroms[i].Split(" ");
                g1chromsInt[i] = new int[parts.Length];
                for (int j = 0; j < parts.Length; j++)
                {
                    g1chromsInt[i][j] = Int32.Parse(parts[j]);
                    if (Math.Abs(g1chromsInt[i][j]) > n)
                    {
                        n = Math.Abs(g1chromsInt[i][j]);
                    }
                }
            }

            string[] g2chroms = g2.Substring(1, g2.Length - 2).Split(")(");
            int[][] g2chromsInt = new int[g2chroms.Length][];
            for (int i = 0; i < g2chroms.Length; i++)
            {
                string[] parts = g2chroms[i].Split(" ");
                g2chromsInt[i] = new int[parts.Length];
                for (int j = 0; j < parts.Length; ++j)
                {
                    g2chromsInt[i][j] = Int32.Parse(parts[j]);
                }
            }

            int[,] bp = new int[2 * n + 1, 2]; // each arrow has 2 nodes (front and back) and each node has 2 edges. Front = 1 to n. Back = n+1 to 2n

            for (int i = 0; i < bp.GetLength(0); i++)
            {
                bp[i, 0] = -1;
                bp[i, 1] = -1;
            }

            for (int i = 0; i < g1chromsInt.Length; i++)
            {
                for (int j = 0; j < g1chromsInt[i].Length; j++)
                {
                    int curr = Math.Abs(g1chromsInt[i][j]);
                    int sign = g1chromsInt[i][j] / curr; // +1 for positive, -1 for negative
                    int next = -1;
                    int nextSign = 0;

                    if (j == g1chromsInt[i].Length - 1)
                    {
                        next = Math.Abs(g1chromsInt[i][0]);
                        nextSign = g1chromsInt[i][0] / next;
                    }
                    else
                    {
                        next = Math.Abs(g1chromsInt[i][j + 1]);
                        nextSign = g1chromsInt[i][j + 1] / next;
                    }
                    int prev = -1;
                    int prevSign = 0;
                    if (j == 0)
                    {
                        prev = Math.Abs(g1chromsInt[i][g1chromsInt[i].Length - 1]);
                        prevSign = g1chromsInt[i][g1chromsInt[i].Length - 1] / prev;
                    }
                    else
                    {
                        prev = Math.Abs(g1chromsInt[i][j - 1]);
                        prevSign = g1chromsInt[i][j - 1] / prev;
                    }

                    // If I'm negative and my next is negative, then my back goes to his forward
                    if (sign == -1 && nextSign == -1)
                    {
                        if (bp[frontToBack(n, curr), 0] == -1)
                        {
                            bp[frontToBack(n, curr), 0] = next;
                        }
                        else
                        {
                            bp[frontToBack(n, curr), 1] = next;
                        }
                    }
                    // If I'm negative and my next is positive, then my back goes to his back
                    else if (sign == -1 && nextSign == 1)
                    {
                        if (bp[frontToBack(n, curr), 0] == -1)
                        {
                            bp[frontToBack(n, curr), 0] = frontToBack(n, next);
                        }
                        else
                        {
                            bp[frontToBack(n, curr), 1] = frontToBack(n, next);
                        }
                    }
                    // If I'm positive and my next is negative, then my forward goes to his forward
                    else if (sign == 1 && nextSign == -1)
                    {
                        if (bp[curr, 0] == -1)
                        {
                            bp[curr, 0] = next;
                        }
                        else
                        {
                            bp[curr, 1] = next;
                        }
                    }
                    // If I'm positive and my next is positive, then my forward goes to his back
                    else if (sign == 1 && nextSign == 1)
                    {
                        if (bp[curr, 0] == -1)
                        {
                            bp[curr, 0] = frontToBack(n, next);
                        }
                        else
                        {
                            bp[curr, 1] = frontToBack(n, next);
                        }
                    }

                    // If I'm negative and my prev is negative, then my forward goes to his back
                    if (sign == -1 && prevSign == -1)
                    {
                        if (bp[curr, 0] == -1)
                        {
                            bp[curr, 0] = frontToBack(n, prev);
                        }
                        else
                        {
                            bp[curr, 1] = frontToBack(n, prev);
                        }
                    }
                    // If I'm negative and my prev is positive, then my forward goes to his forward
                    else if (sign == -1 && prevSign == 1)
                    {
                        if (bp[curr, 0] == -1)
                        {
                            bp[curr, 0] = prev;
                        }
                        else
                        {
                            bp[curr, 1] = prev;
                        }
                    }
                    // If I'm positive and my prev is negative, then my back goes to his back
                    else if (sign == 1 && prevSign == -1)
                    {
                        if (bp[frontToBack(n, curr), 0] == -1)
                        {
                            bp[frontToBack(n, curr), 0] = frontToBack(n, prev);
                        }
                        else
                        {
                            bp[frontToBack(n, curr), 1] = frontToBack(n, prev);
                        }
                    }
                    // If I'm positive and my prev is positive, then my back goes to his forward
                    else if (sign == 1 && prevSign == 1)
                    {
                        if (bp[frontToBack(n, curr), 0] == -1)
                        {
                            bp[frontToBack(n, curr), 0] = prev;
                        }
                        else
                        {
                            bp[frontToBack(n, curr), 1] = prev;
                        }
                    }
                }
            }

            for (int i = 0; i < g2chromsInt.Length; i++)
            {
                for (int j = 0; j < g2chromsInt[i].Length; j++)
                {
                    int curr = Math.Abs(g2chromsInt[i][j]);

                    int sign = g2chromsInt[i][j] / curr; // +1 for positive, -1 for negative
                    int next = -1;
                    int nextSign = 0;
                    if (j == g2chromsInt[i].Length - 1)
                    {
                        next = Math.Abs(g2chromsInt[i][0]);
                        nextSign = g2chromsInt[i][0] / next;
                    }
                    else
                    {
                        next = Math.Abs(g2chromsInt[i][j + 1]);
                        nextSign = g2chromsInt[i][j + 1] / next;
                    }
                    int prev = -1;
                    int prevSign = 0;
                    if (j == 0)
                    {
                        prev = Math.Abs(g2chromsInt[i][g2chromsInt[i].Length - 1]);
                        prevSign = g2chromsInt[i][g2chromsInt[i].Length - 1] / prev;
                    }
                    else
                    {
                        prev = Math.Abs(g2chromsInt[i][j - 1]);
                        prevSign = g2chromsInt[i][j - 1] / prev;
                    }

                    // If I'm negative and my next is negative, then my back goes to his forward
                    if (sign == -1 && nextSign == -1)
                    {
                        if (bp[frontToBack(n, curr), 0] == -1)
                        {
                            bp[frontToBack(n, curr), 0] = next;
                        }
                        else
                        {
                            bp[frontToBack(n, curr), 1] = next;
                        }
                    }
                    // If I'm negative and my next is positive, then my back goes to his back
                    else if (sign == -1 && nextSign == 1)
                    {
                        if (bp[frontToBack(n, curr), 0] == -1)
                        {
                            bp[frontToBack(n, curr), 0] = frontToBack(n, next);
                        }
                        else
                        {
                            bp[frontToBack(n, curr), 1] = frontToBack(n, next);
                        }
                    }
                    // If I'm positive and my next is negative, then my forward goes to his forward
                    else if (sign == 1 && nextSign == -1)
                    {
                        if (bp[curr, 0] == -1)
                        {
                            bp[curr, 0] = next;
                        }
                        else
                        {
                            bp[curr, 1] = next;
                        }
                    }
                    // If I'm positive and my next is positive, then my forward goes to his back
                    else if (sign == 1 && nextSign == 1)
                    {
                        if (bp[curr, 0] == -1)
                        {
                            bp[curr, 0] = frontToBack(n, next);
                        }
                        else
                        {
                            bp[curr, 1] = frontToBack(n, next);
                        }
                    }

                    // If I'm negative and my prev is negative, then my forward goes to his back
                    if (sign == -1 && prevSign == -1)
                    {
                        if (bp[curr, 0] == -1)
                        {
                            bp[curr, 0] = frontToBack(n, prev);
                        }
                        else
                        {
                            bp[curr, 1] = frontToBack(n, prev);
                        }
                    }
                    // If I'm negative and my prev is positive, then my forward goes to his forward
                    else if (sign == -1 && prevSign == 1)
                    {
                        if (bp[curr, 0] == -1)
                        {
                            bp[curr, 0] = prev;
                        }
                        else
                        {
                            bp[curr, 1] = prev;
                        }
                    }
                    // If I'm positive and my prev is negative, then my back goes to his back
                    else if (sign == 1 && prevSign == -1)
                    {
                        if (bp[frontToBack(n, curr), 0] == -1)
                        {
                            bp[frontToBack(n, curr), 0] = frontToBack(n, prev);
                        }
                        else
                        {
                            bp[frontToBack(n, curr), 1] = frontToBack(n, prev);
                        }
                    }
                    // If I'm positive and my prev is positive, then my back goes to his forward
                    else if (sign == 1 && prevSign == 1)
                    {
                        if (bp[frontToBack(n, curr), 0] == -1)
                        {
                            bp[frontToBack(n, curr), 0] = prev;
                        }
                        else
                        {
                            bp[frontToBack(n, curr), 1] = prev;
                        }
                    }
                }
            }

            for (int k = 0; k < bp.GetLength(0); k++)
            {
                int first = bp[k, 0];
                int second = bp[k, 1];
                breakpointGraph.Add(Tuple.Create(first, second));
            }

            return breakpointGraph;
        }

        public static int frontToBack(int n, int i)
        {
            return i + n;
        }

        public static int backToFront(int n, int i)
        {
            return i - n;
        }


        public static int TwoBreakDistance(string P, string Q)
        {
            int numBreaks = 0;

            List<Tuple<int, int>> combinedEdgeList = breakpointGraph(P, Q);

            int mycylces = GetCycles(combinedEdgeList);

            int n = (combinedEdgeList.Count() - 1) / 2;
            numBreaks =  n - mycylces;

            return numBreaks;
        }

        public static string GraphToGenome(string GenomeGraph)
        {
            string P = "";
            string[] cycle = GenomeGraph.Split("), ");
            List<Tuple<int, int>> coloredList = new List<Tuple<int, int>>();

            foreach (string str in cycle)
            {
                string val = str.Substring(1, str.Length - 1);

                if (val[val.Length - 1] == ')')
                    val = val.Substring(0, val.Length - 1);

                string [] items = val.Split(", ");
                coloredList.Add(Tuple.Create(Int32.Parse(items[0]), Int32.Parse(items[1])));
            }

            //            List<Tuple<int, int>> cycleList = new List<Tuple<int, int>>();

//            List<Tuple<int, int>> cycleList1 = new List<Tuple<int, int>>();

//           cycleList1 = getCycles(coloredList);
 
            coloredList.Add(coloredList[0]);
            coloredList.Insert(0, coloredList[coloredList.Count() - 2]);

            List<List<Tuple<int,int>>> cycleList = new List<List<Tuple<int,int>>>();
            List<Tuple<int, int>> valList = new List<Tuple<int, int>>();

            for (int i = 1; i < coloredList.Count - 1; i++)
            {
                int diff = Math.Abs(coloredList[i].Item2 - coloredList[i + 1].Item1);
                valList.Add(Tuple.Create(coloredList[i].Item1, coloredList[i].Item2));

                if (diff > 1)   // end of cycle
                {
                    valList.Insert(0, valList[valList.Count() - 1]);    // hack for connecting later on
                    List<Tuple<int, int>> tmpList = new List<Tuple<int, int>>(valList);

                    cycleList.Add(tmpList);

                    valList.Clear();
                }
            }

            for (int j = 0; j < cycleList.Count(); j++)
            {
                valList = cycleList[j];
                string val3 = "(";

                for (int k = 0; k < valList.Count() - 1; k++)
                {
                    val3 += valList[k].Item2.ToString() + " " + valList[k + 1].Item1.ToString() + " ";
                }
                val3 = val3.TrimEnd();

                val3 += ")";

                P += CycleToChromosome(val3) + " ";
            }   

            P = P.TrimEnd();

            return P;
        }
        public static string ColoredEdges(string input)
        {
            string Edges = "";
            string Nodes = "";
            string[] chromosome = input.Split(")");
            List<Tuple<string, string>> edgeList = new List<Tuple<string, string>>();

            foreach (string s in chromosome)
            {
                if (s.Length > 0)
                {
                    string val = s.Substring(1, s.Length - 1);
                    string[] val2 = val.Split(" ");          

                    Nodes = ChromosomeToCycle(val + " " + val2[0]);

                    Nodes = Nodes.Substring(1, Nodes.Length - 2);
                    string[] nodeArray = Nodes.Split(" ");

                    for (int j = 1; j < nodeArray.Length - 1; j+= 2)
                    {
                        edgeList.Add(Tuple.Create(nodeArray[j], nodeArray[j + 1]));
                    }   
                }
            }
            
            foreach (var ed in edgeList)
            {
                Edges += "(" + ed.Item1 + ", " + ed.Item2 + "), ";
            }

            Edges = Edges.Substring(0, Edges.Length - 2);

            return Edges;
        }

        public static string CycleToChromosome(string input)
        {
            string Chromosome = "";
            string Nodes = input.Substring(1, input.Length - 2);
            string[] permutations = Nodes.Split(" ");
            List<int> myList = new List<int>();

            for (int j = 1; j < permutations.Length; j+= 2)
            {
                int prevVal = Int32.Parse(permutations[j - 1]);
                int Val = Int32.Parse(permutations[j]);

                if (prevVal < Val)
                    myList.Add(Val / 2);
                else
                    myList.Add((prevVal / 2) * -1);
            }

            Chromosome = "(";
            foreach (int k in myList)
            {
                if (k > 0)
                    Chromosome += "+" + k.ToString() + " ";
                else
                    Chromosome += k.ToString() + " ";
            }

            Chromosome = Chromosome.TrimEnd();
            Chromosome += ")";

            return Chromosome;

        }

        public static string ChromosomeToCycle(string input)
        {
            string Nodes = "";
            string Chromosome = "";

            if (input[0].ToString() == "(")
                Chromosome = input.Substring(1, input.Length - 2);
            else
                Chromosome = input;
        
            string[] permutations = Chromosome.Split(" ");
            List<int> myList = new List<int>();

            for (int j = 0; j < permutations.Length; j++)
            {
                int value = Int32.Parse(permutations[j]);
                if (value > 0)
                {
                    myList.Add((value * 2) - 1);
                    myList.Add(value * 2);
                }
                else
                {
                    myList.Add(value * -2);
                    myList.Add((value * -2) - 1);
                }
            }

            Nodes = "(";
            foreach (int k in myList)
            {
                Nodes += k.ToString() + " ";
            }

            Nodes = Nodes.TrimEnd();
            Nodes += ")";

            return Nodes;
        }

        public static int NumberOfBreakpoints(string P)
        {
            int numBreakpoints = 0;
            string input = P.Substring(1, P.Length - 2);
            string[] permutations = input.Split(" ");
            List<int> myList = new List<int>();

            myList.Add(0);
            foreach (string s in permutations)
            {
                int value = Int32.Parse(s);
                myList.Add(value);
            }
            myList.Add(permutations.Length + 1);

            for (int i = 1; i < myList.Count(); i++)
            {
                int prevVal = myList[i - 1];
                int val = myList[i];

                if (prevVal > val || prevVal < val - 1)
                    numBreakpoints++;
            }

            return numBreakpoints;
        }

        public static string GreedySortingOutput(List<int> sortedList)
        {
            string output = "";
            output += "(";

            foreach(int i in sortedList)
            {
                if (i < 0)
                    output += i.ToString() + " ";
                else
                {
                    output += "+" + i.ToString() + " ";
                }
            }
            output = output.TrimEnd();
            output += ")";

            return output;
        }

        public static List<string> GreedySorting(string P)
        {
            List<string> result = new List<string>();
            int approxReversalDistance = 0;
            string input = P.Substring(1, P.Length - 2);
            string[] permutations = input.Split(" ");
            List<int> sortedList = new List<int>();
            string output = "";

            foreach (string s in permutations)
            {
                int value = Int32.Parse(s);
                sortedList.Add(value);
            }

            for (int k = 1; k <= sortedList.Count(); k++)
            {
                if (Math.Abs(sortedList[k-1]) == k)
                {
                    if (sortedList[k - 1] == k * -1)
                    {
                        sortedList[k - 1] *= -1;
                        approxReversalDistance++;
                        output = GreedySortingOutput(sortedList);
                        result.Add(output);
                    }
                }
                else
                {
                    bool reverseSign = false;
                    List<int> tmpList = new List<int>();

                    for (int j = k; j <= sortedList.Count(); j++)
                    {
                        int v = sortedList[j - 1];
                        tmpList.Add(v * - 1);
                        if (Math.Abs(v) == k)
                            break;
                    }
                    if (tmpList[tmpList.Count() - 1] == k * -1)
                        reverseSign = true;

                    tmpList.Reverse();
                    int tmpIndex = 0;
                    int nCount = k + tmpList.Count();
                    for (int t = k ; t < nCount; t++)
                    {
                        sortedList[t - 1] = tmpList[tmpIndex++];
                    }
                    approxReversalDistance++;
                    output = GreedySortingOutput(sortedList);
                    result.Add(output);


                    if (reverseSign)
                    {
                        sortedList[k - 1] *= -1;
                        approxReversalDistance++;
                        output = GreedySortingOutput(sortedList);
                        result.Add(output);
                    }
                }
            }


            return result;
        }

        public static List<string> LongestPathInDAG(int startingNode, int endingNode, List<string> edgeList)
        {
            List<string> result = new List<string>();
            HashSet<int> nodeSet = new HashSet<int>();
            HashSet<Tuple<int, int, int>> edgeSet = new HashSet<Tuple<int, int, int>>();
            int maxNodeValue = 0;
            Dictionary<int, List<Tuple<int, int>>> edgeDict = new Dictionary<int, List<Tuple<int, int>>>();

            // Parse and create Nodes and Edges HashSets
            foreach (string edge in edgeList)
            {
                string[] split1 = edge.Split("->");
                string[] split2 = split1[1].Split(":");
                int nodeVal = 0;
                int edgeVal = 0;
                int weight = 0;
                Int32.TryParse(split1[0], out nodeVal);
                Int32.TryParse(split2[0], out edgeVal);
                Int32.TryParse(split2[1], out weight);

                if (edgeVal >= startingNode && edgeVal <= endingNode)
                { 
                    if (edgeDict.ContainsKey(edgeVal))
                    {
                        List<Tuple<int, int>> newItem = new List<Tuple<int, int>>();
                        newItem = edgeDict[edgeVal];
                        newItem.Add(Tuple.Create(nodeVal, weight));
                        edgeDict[edgeVal] = newItem;

                    }
                    else
                    {
                        List<Tuple<int, int>> newItem = new List<Tuple<int, int>>();
                        newItem.Add(Tuple.Create(nodeVal, weight));
                        edgeDict.Add(edgeVal, newItem);
                    }
                }

                bool ret = false;

                if (edgeVal > maxNodeValue)
                    maxNodeValue = edgeVal;

                ret = nodeSet.Add(nodeVal);
                ret = edgeSet.Add(Tuple.Create(nodeVal, edgeVal, weight));
            }

            int nodeSize = maxNodeValue + 1;
            int[,] DAGGraph = new int[nodeSize, nodeSize];
            int[] dist = new int[DAGGraph.GetUpperBound(0) + 1];

            for (int i = 0; i < nodeSize; i++)
                for (int j = 0; j < nodeSize; j++)
                    DAGGraph[i, j] = INF;

            foreach(var s in edgeSet)
            {
                DAGGraph[s.Item1, s.Item2] = s.Item3;
            }

            int longestPath = LongestPathLength(startingNode, endingNode, DAGGraph, dist);

            result.Add(longestPath.ToString());


            List<string> pathList = new List<string>();
            int currentNode = endingNode;
            var Edge = edgeDict[currentNode];

            pathList.Add(currentNode.ToString());

            while (Edge.Count > 0)
            {
                List<Tuple<int, int>> connectedEdges = new List<Tuple<int, int>>();
                int prevNode = 0;
                int prevVal = -1;

                connectedEdges = Edge;

                foreach (var s in connectedEdges)
                {
                    string path = currentNode.ToString() + "->" + s.Item1.ToString();

                    if (dist[s.Item1] > prevVal && s.Item1 >= startingNode && s.Item1 <= endingNode)
                    {
                        prevVal = dist[s.Item1];
                        prevNode = s.Item1;

                    }
                }
                pathList.Add(prevNode.ToString());
                currentNode = prevNode;

                if (prevNode == startingNode)
                    break;
                else if (edgeDict.ContainsKey(prevNode))
                    Edge = edgeDict[prevNode];
                else
                {
                    if (prevNode == startingNode)
                        break;
                    else
                    {
                        currentNode = endingNode;
                        Edge = edgeDict[endingNode];    // start again and skip paths used
                        pathList.Clear();
                    }
                 }
            }
            pathList.Reverse();

            string ReversePath = "";
            for (int i = 0; i < pathList.Count(); i++)
                ReversePath += pathList[i] + "->";

            ReversePath = ReversePath.Substring(0, ReversePath.Length - 2);  // remove last ->

            result.Add(ReversePath);

            return result;
        }


        public static void TopoSort(int u, bool[] visited, Stack<int> stk, int[,] DAGGraph)
        {
            visited[u] = true;
            int nodeLength = DAGGraph.GetUpperBound(0);

            for (int v = 0; v < nodeLength; v++)
            {
                if (DAGGraph[u,v] != INF)
                {
                    if (!visited[v])
                        TopoSort(v, visited, stk, DAGGraph);
                }
            }

            stk.Push(u);
        }

        public static int LongestPathLength(int start, int end, int[,] DAGGraph, int[] dist)
        {
            Stack<int> stk = new Stack<int>();
            int nodeSize = DAGGraph.GetUpperBound(0) + 1;
            bool[] vis = new bool[nodeSize];
            int result = 0;

            for (int i = 0; i < nodeSize; i++)
                vis[i] = false;

            for (int i = 0; i < nodeSize; i++)
            {
                if (!vis[i])
                    TopoSort(i, vis, stk, DAGGraph);
            }

            for (int i = 0; i < nodeSize; i++)
            {
                dist[i] = INF;
            }

            dist[start] = 0;

            while (stk.Count() > 0)
            {
                int nextVert = stk.Pop();

                if (dist[nextVert] != INF)
                {
                    for (int v = 0; v < nodeSize; v++)
                    {
                        if (DAGGraph[nextVert,v] != INF)
                        {
                            if (dist[v] < dist[nextVert] + DAGGraph[nextVert, v])
                            {
                                dist[v] = dist[nextVert] + DAGGraph[nextVert, v];
                            }
                        }
                    }
                }
            }

            result = dist[end];

            return result;
        }


        /// <summary>
        /// Topological Sorting (Kahn's algorithm) 
        /// </summary>
        /// <remarks>https://en.wikipedia.org/wiki/Topological_sorting</remarks>
        /// <typeparam name="T"></typeparam>
        /// <param name="nodes">All nodes of directed acyclic graph.</param>
        /// <param name="edges">All edges of directed acyclic graph.</param>
        /// <returns>Sorted node in topological order.</returns>
        public static List<T> TopologicalSort<T>(HashSet<T> nodes, HashSet<Tuple<T, T, T>> edges) where T : IEquatable<T>
        {
            // Empty list that will contain the sorted elements
            var L = new List<T>();

            // Set of all nodes with no incoming edges
            var S = new HashSet<T>(nodes.Where(n => edges.All(e => e.Item2.Equals(n) == false)));

            // while S is non-empty do
            while (S.Any())
            {

                //  remove a node n from S
                var n = S.First();
                S.Remove(n);

                // add n to tail of L
                L.Add(n);

                // for each node m with an edge e from n to m do
                foreach (var e in edges.Where(e => e.Item1.Equals(n)).ToList())
                {
                    var m = e.Item2;

                    // remove edge e from the graph
                    edges.Remove(e);

                    // if m has no other incoming edges then
                    if (edges.All(me => me.Item2.Equals(m) == false))
                    {
                        // insert m into S
                        S.Add(m);
                    }
                }
            }

            // if graph has edges then
            if (edges.Any())
            {
                // return error (graph has at least one cycle)
                return null;
            }
            else
            {
                // return L (a topologically sorted order)
                return L;
            }
        }

        public static string OutputLCS(string [,] backtrack, string v, int i, int j)
        {
            string LCS = "";

            if (i == 0 || j == 0)
                return LCS;

            if (backtrack[i, j] == down_arrow)
                LCS += OutputLCS(backtrack, v, i - 1, j);
            else if (backtrack[i, j] == right_arrow)
                LCS += OutputLCS(backtrack, v, i, j - 1);
            else
            {
                LCS += OutputLCS(backtrack, v, i - 1, j - 1);
                LCS += v[i-1].ToString();
            }

            return LCS;
        }

        public static string[,] LCSBackTrack(string v, string w)
        {
            int downLength = v.Length + 1;
            int rightLength = w.Length + 1;
            int[,] pathArray = new int[downLength, rightLength];

            string[,] Backtrack = new string[downLength, rightLength];

            for (int i = 0; i < downLength; i++)
            {
                for (int k = 0; k < rightLength; k++)
                {
                    Backtrack[i, k] = blank;
                    pathArray[i, k] = 0;
                }
            }
            
            for (int i = 1; i < downLength; i++)
            {
                for (int j = 1; j < rightLength; j++)
                {
                    int dwn = pathArray[i - 1, j];
                    int rt = pathArray[i, j - 1];
                    int eq = pathArray[i - 1, j - 1];
                    if (v[i-1] == w[j-1])
                    {
                        eq += 1;
                    }
                    pathArray[i,j] = Max(dwn, rt, eq);

                    if (pathArray[i, j] == pathArray[i - 1, j])
                        Backtrack[i, j] = down_arrow;
                    else if (pathArray[i, j] == pathArray[i, j - 1])
                        Backtrack[i, j] = right_arrow;
                    else if (pathArray[i, j] == (pathArray[i - 1, j - 1] + 1) && v[i-1] == w[j-1])
                        Backtrack[i, j] = diag_arrow;
                }
            }

            return Backtrack;
        }

        public static int ManhattanTourist(int n, int m, int[,] DownList, int[,] RightList)
        {
            int result = 0;
            int[,] pathArray = new int[n + 1, m + 1];

            pathArray[0, 0] = 0;


            for (int i = 1; i <= n; i++)
                pathArray[i, 0] = pathArray[i - 1, 0] + DownList[i-1, 0];

            for (int j = 1; j <= m; j++)
                pathArray[0, j] = pathArray[0, j - 1] + RightList[0, j-1];

            for (int i = 1; i <= n; i++)
            {
                for (int j = 1; j <= m; j++)
                {
                    int downVal = 0;
                    int rightVal = 0;

                    downVal = pathArray[i - 1, j] + DownList[i-1, j];
                    rightVal = pathArray[i, j - 1] + RightList[i, j-1];

                    if (downVal >= rightVal)
                        pathArray[i, j] = downVal;
                    else
                        pathArray[i, j] = rightVal;
                }
            }

            result = pathArray[n, m];

            return result;
        }

        public static int DPChange(int money, List<int> Coins)
        {
            int numCoins = 0;
            int[] MinNumCoins = new int[money + 1];

            MinNumCoins[0] = 0;

            for (int m = 1; m <= money; m++)
            {
                MinNumCoins[m] = int.MaxValue;
                for (int i = 0; i < Coins.Count(); i++)
                {
                    int coin = Coins[i];
                    if (m >= coin)
                    {
                        if (MinNumCoins[m - coin] + 1 < MinNumCoins[m])
                            MinNumCoins[m] = MinNumCoins[m - coin] + 1;
                    }
                }
            }

            numCoins = MinNumCoins[money];

            return numCoins;
        }

        public static string SpectralConvolution(string spectrum)
        {
            string result = "";
            string[] array = spectrum.Split(" ");
            List<int> spectrumList = new List<int>();
            int val = 0;
            int totalMass = 0;
            int convolution = 0;

            foreach (string s in array)
            {
                Int32.TryParse(s, out val);
                totalMass += val;
                spectrumList.Add(val);
            }
            spectrumList.Add(totalMass);

            spectrumList.Sort();

            for (int i = 0; i < spectrumList.Count() - 1; i++)
            {
                for (int j = 0; j < spectrumList.Count - 1; j++)
                {
                    convolution = spectrumList[i] - spectrumList[j];
                    if (convolution > 0)
                        result += convolution.ToString() + " ";
                }
            }

            result = result.Trim();

            return result;
        }

        public static int LinearScore(string peptide, string spectrum, bool cyclic = false)
        {
            int result = 0;
            List<int> peptideList = new List<int>();

            string[] spectrumArray = spectrum.Split(" ");
            List<int> spectrumList = new List<int>();

            if (cyclic)
                peptideList = CyclicSpectrum(peptide);
            else
                peptideList = LinearSpectrum(peptide);

            foreach (string s in spectrumArray)
            {
                int val = 0;
                Int32.TryParse(s, out val);

                spectrumList.Add(val);
            }

            for (int i = 0; i < peptideList.Count(); i++)
            {
                for (int j = 0; j < spectrumList.Count(); j++)
                {
                    if (peptideList[i] == spectrumList[j])
                    {
                        result += 1;
                        spectrumList.RemoveAt(j);
                        break;
                    }
                }
            }
        
            return result;
        }
        
        public static List<string> TrimLeaderboard(List <string> leaderboard, string spectrum, int N)
        {
            List<string> newLeaderboard = new List<string>();
            List<Peptide> tmpBoard = new List<Peptide>();
            string[] spectrumArray = spectrum.Split(" ");
            Dictionary<string, int> scoreLookup = new Dictionary<string, int>();

            int currentScore = 0;
            string currentPeptide = "";

            for (int i = 0; i < leaderboard.Count(); i++)
            {
                currentPeptide = leaderboard[i];

                if (scoreLookup.ContainsKey(currentPeptide))
                    currentScore = scoreLookup[currentPeptide];
                else
                    currentScore = LinearScore(currentPeptide, spectrum, true);

                Peptide tmpPeptide = new Peptide(currentPeptide, currentScore);
                tmpBoard.Add(tmpPeptide);
            }

            MySortingClass sort1 = new MySortingClass();
            tmpBoard.Sort(sort1);
            tmpBoard.Reverse();

            for (int j = 0; j < N && j < tmpBoard.Count(); j++)
            {
                currentScore = tmpBoard[j].score;
                newLeaderboard.Add(tmpBoard[j].code);
            }

            for (int k = N; k < tmpBoard.Count; k++)
            {
                if (tmpBoard[k].score == currentScore)
                    newLeaderboard.Add(tmpBoard[k].code);
                else if (k > N)
                    break;      // no need to check further
            }
            
            return newLeaderboard;
        }

        public static string LeaderboardCyclopeptideSequencing(string spectrum, int N, int M = 0)
        {
//            string result = "";
            string LeaderPeptide = "";
            List<string> leaderboard = new List<string>();
            List<int> spectrumList = new List<int>();
            Dictionary<string, int> scoreLookup = new Dictionary<string, int>();

            int currentScore = 0;
            int leaderScore = 0;
            List<Int64> massArray = new List<Int64>();

            string[] spectrumArray = spectrum.Split(" ");
            foreach (string s in spectrumArray)
            {
                Int32 val = 0;
                Int32.TryParse(s, out val);
                spectrumList.Add(val);
            }

            if (M > 0)  // Use Convolution
            {
                string convolution = SpectralConvolution(spectrum);
                string[] array = convolution.Split(" ");
                int val = 0;
                List<Peptide> convolutionList = new List<Peptide>();
                Dictionary<string, int> dict = new Dictionary<string, int>();

                foreach (string s in array)
                {
                    Int32.TryParse(s, out val);

                    if (val >= 57 && val <= 200)
                    {
                        if (dict.ContainsKey(s))
                            dict[s] += 1;
                        else
                            dict.Add(s, 1);
                    }
                }

                foreach (var s in dict)
                {
                    Peptide pep1 = new Peptide(s.Key, s.Value);
                    convolutionList.Add(pep1);
                }

                MySortingClass sort1 = new MySortingClass();
                convolutionList.Sort(sort1);
                convolutionList.Reverse();

                for (int j = 0; j < M && j < convolutionList.Count(); j++)
                {
                    currentScore = convolutionList[j].score;
                    Int32.TryParse(convolutionList[j].code, out val);

                    massArray.Add(val);
                }

                for (int k = M; k < convolutionList.Count; k++)
                {
                    if (convolutionList[k].score == currentScore)
                    {
                        Int32.TryParse(convolutionList[k].code, out val);
                        massArray.Add(val);
                    }
                    else if (k > M)
                        break;      // no need to check further
                }

                massArray.Sort();
            }
            else
                massArray = peptideMass;


            int parentMass = ParentMass(spectrum);

            leaderboard.Add("");

            while (leaderboard.Count > 0)
            {
                leaderboard = PurgePeptideArray(leaderboard);

                leaderboard = ExpandPeptides(leaderboard, massArray);
                for (int i = 0; i < leaderboard.Count; i++)
                {
                    string peptide = leaderboard[i];

                    Int32 peptideMass = PeptideMass(peptide);

                    if (peptideMass == parentMass)
                    {
                        if (scoreLookup.ContainsKey(peptide))
                            currentScore = scoreLookup[peptide];
                        else
                        {
                            currentScore = LinearScore(peptide, spectrum, true);
                            scoreLookup.Add(peptide, currentScore);
                        }

                        if (scoreLookup.ContainsKey(LeaderPeptide))
                            leaderScore = scoreLookup[LeaderPeptide];
                        else
                        {
                            leaderScore = LinearScore(LeaderPeptide, spectrum, true);
                            scoreLookup.Add(LeaderPeptide, leaderScore);
                        }

                        if (currentScore > leaderScore)
                        {
                            LeaderPeptide = peptide;
                        }

                    }
                    else if (peptideMass > parentMass)
                    {
                        leaderboard[i] = "-1";

                    }
                }

                leaderboard = PurgePeptideArray(leaderboard);       // remove -1 entrees

                leaderboard = TrimLeaderboard(leaderboard, spectrum, N);
            }

            return LeaderPeptide;
        }

        public static string CyclopeptideSequencing(string spectrum)
        {
            string result = "";
            List<string> peptideArray = new List<string>();
            List<int> spectrumList = new List<int>();
            string[] spectrumArray = spectrum.Split(" ");
            foreach (string s in spectrumArray)
            {
                Int32 val = 0;
                Int32.TryParse(s, out val);
                spectrumList.Add(val);
            }

            int parentMass = ParentMass(spectrum);

            peptideArray.Add("");

            while (peptideArray.Count > 0)
            {
                peptideArray = PurgePeptideArray(peptideArray);

                peptideArray = ExpandPeptides(peptideArray, peptideMass);
                for (int i = 0; i < peptideArray.Count; i++)
                {
                    string peptide = peptideArray[i];

                    Int32 peptideMass = PeptideMass(peptide);

                    if (peptideMass == parentMass)
                    {
                        List<int> cycloList = CyclicSpectrum(peptide);

                        if (CompareSpectrum(cycloList, spectrumList))
                        {
                            result += peptide + " ";
                        }

                        peptideArray[i] = "-1";
                    }
                    else if (!IsPeptideConsistentWithSpectrum(peptide, spectrumList))
                    {
                        peptideArray[i] = "-1";

                    }
                }
            }

            return result;
        }

        public static List<string> PurgePeptideArray(List<string> inputList)
        {
            List<string> outputList = new List<string>();

            foreach (string s in inputList)
            {
                if (s != "-1")
                    outputList.Add(s);
            }

            return outputList;
        }

        public static bool IsPeptideConsistentWithSpectrum(string peptide, List<int> spectrum)
        {
            bool IsConsistent = false;
            List<int> theoriticalList = new List<int>();
            Dictionary<int, int> counts = new Dictionary<int, int>();
            List<int> peptideList = new List<int>();

            peptideList = LinearSpectrum(peptide);

            foreach (int p in peptideList)
            {
                if (counts.ContainsKey(p))
                    counts[p] += 1;
                else
                    counts.Add(p, 1);
            }

            foreach (var val in counts)
            {
                int count = val.Value;
                int countsFound = 0;
                bool match = false;

                foreach(int sp in spectrum)
                {
                    if (sp == val.Key)
                    {
                        match = true;
                        countsFound++;
                    }
                }

                if (match)
                {
                    if (countsFound >= count)
                        IsConsistent = true;
                    else
                    {
                        IsConsistent = false;
                        break;
                    }
                }
                else
                {
                    IsConsistent = false;
                    break;
                }
            }
            return IsConsistent;
        }



        public static bool CompareSpectrum(List<int> spectrum1, List<int> spectrum2)
        {
            bool equal = false;

            if (spectrum1.Count() == spectrum2.Count())
            {
                for (int i = 0; i < spectrum1.Count(); i++)
                {
                    if (spectrum1[i] != spectrum2[i])
                    {
                        equal = false;
                        break;
                    }
                    else
                        equal = true;
                }
            }

            return equal;
        }

        public static int ParentMass(string spectrum)
        {
            int parentMass = 0;
            string[] spectrumArray = spectrum.Split(" ");
            List<int> parent = new List<int>();

            foreach (string s in spectrumArray)
            {
                Int32 m = 0;
                Int32.TryParse(s, out m);
                parent.Add(m);
            }

            parent.Sort();
            parentMass = parent[parent.Count - 1];

            return parentMass;
        }

        public static int PeptideMass(string peptide)
        {
            string[] split = peptide.Split("-");
            Int32 peptideMass = 0;
            Int32 mass = 0;

            foreach (string s in split)
            {
                Int32.TryParse(s, out mass);

                peptideMass += mass;
            }

            return peptideMass;
        }


        public static List<string> ExpandPeptides(List<string> peptideArray, List<Int64> massArray)
        {
            List<string> expandedArray = new List<string>();

            foreach (string p in peptideArray)
            {
                foreach (int m in massArray)
                {
                    string newVal = "";
                    if (p.Length > 0)
                        newVal = p + "-" + m.ToString();
                    else
                        newVal = m.ToString();

                    expandedArray.Add(newVal);
                }
            }

            return expandedArray;
        }

        public static Int64 BFCyclopeptideSequencing(int spectrum)
        { 
              Int64 result = 0; 
//              List<int> arr = new List<int> { 1, 2, 3 }; 
  
              result = MassCounting(spectrum, peptideMass); 
  
              return result; 
          } 
  
 
          public static Int64 MassCounting(Int64 x, List<Int64> c)
          { 
              Int64[] d = new Int64[x + 1]; 
              d[0] = 1; 
              Int64 k = c.Count(); 
  
 
              for (int i = 1; i <= x; i++) 
              { 
                  for (int j = 0; j<k; j++) 
                  { 
                      if (i - c[j] < 0) 
                          continue; 
                      d[i] += d[i - c[j]]; 
                  } 
              } 
              return d[x]; 
          } 

        public static List<int> TheoreticalSpectrum(string peptide)
        {
            List<int> result = new List<int>();
            List<string> subPeptides = new List<string>();

            result.Add(0);

            subPeptides = CreateSubpeptides(peptide);

            foreach (string p in subPeptides)
            {
                List<int> mass = new List<int>();

                mass = LinearSpectrum(p);

                result.Add(mass[mass.Count - 1]);
            }

            Int64 total = MassPeptide(peptide);

//            foreach (char c in peptide)
//                total += integerMass[c.ToString()];

            result.Add((int)total);

            result.Sort();

            return result;
        }

        public static Int64 MassPeptide(string peptide)
        {
            Int64 mass = 0;

            foreach (char c in peptide)
                mass += integerMass[c.ToString()];

            return mass;
        }

        public static List<string> CreateSubpeptides(string peptide)
        {
            List<string> subPeptides = new List<string>();
            int k = peptide.Length;
            string doublePeptid = peptide + peptide;    // help with wraparound
//            int length = 1;

            for (int i = 1; i < k; i++)
            {
                for (int j = 0; j < k; j++)
                {
                    subPeptides.Add(doublePeptid.Substring(j, i));
                }
            }

            return subPeptides;
        }

        public static string ConvertNumbericPeptide(string peptide)
        {
            string newPeptide = "";

            return newPeptide;
        }

         public static List<int> CyclicSpectrum(string peptide)
        {
            List<int> cycloSpectrum = new List<int>();
            int peptideMass = 0;
            List<int> prefixMassList = new List<int>();
            int index = 0;
            int dummy = 0;

            string[] peptideArray = peptide.Split("-");

            prefixMassList.Insert(index++, 0);

            if (peptideArray.Length > 0 && int.TryParse(peptideArray[0], out dummy))      // numeric peptide
            {
                for (int i = 0; i < peptideArray.Length; i++)
                {
                    int tmp = 0;
                    Int32.TryParse(peptideArray[i], out tmp);

                    int mass = prefixMassList[i] + tmp;
                    prefixMassList.Insert(index++, mass);
                }
            }
            else   // regular peptide
            {
                for (int i = 0; i < peptide.Length; i++)
                {
                    char p = peptide[i];

                    foreach (string j in integerMass.Keys)
                    {
                        if (j == peptide[i].ToString())
                        {
                            int mass = prefixMassList[i] + integerMass[j];
                            prefixMassList.Insert(index++, mass);
                            break;
                        }
                    }
                }
            }

            peptideMass = prefixMassList[prefixMassList.Count - 1];

            index = 0;
            cycloSpectrum.Insert(index++, 0);

            for (int i = 0; i < prefixMassList.Count; i++)
            {
                for (int j = i + 1; j < prefixMassList.Count; j++)
                {
                    int newMass = prefixMassList[j] - prefixMassList[i];
                    cycloSpectrum.Insert(index++, newMass);

                    if (i > 0 && (j < prefixMassList.Count - 1))
                    {
                        int v = 0;
                        v = peptideMass - newMass;
                        cycloSpectrum.Insert(index++, v);
                    }
                }
            }

            cycloSpectrum.Sort();

            return cycloSpectrum;
        }
        public static List<int> LinearSpectrum(string peptide)
        {
            List<int> linearSpectrum = new List<int>();
            List<int> prefixMassList = new List<int>();
            int index = 0;
            int dummy = 0;

            string[] peptideArray = peptide.Split("-");

            prefixMassList.Insert(index++, 0);

            if (peptideArray.Length > 0 && int.TryParse(peptideArray[0], out dummy))      // numeric peptide
            {
                for (int i = 0; i < peptideArray.Length; i++)
                {
                    int tmp = 0;
                    Int32.TryParse(peptideArray[i], out tmp);

                    int mass = prefixMassList[i] + tmp;
                    prefixMassList.Insert(index++, mass);
                }
            }
            else   // regular peptide
            {
                for (int i = 0; i < peptide.Length; i++)
                {
                    char p = peptide[i];

                    foreach (string j in integerMass.Keys)
                    {
                        if (j == peptide[i].ToString())
                        {
                            int mass = prefixMassList[i] + integerMass[j];
                            prefixMassList.Insert(index++, mass);
                            break;
                        }
                    }
                }
            }

            index = 0;
            linearSpectrum.Insert(index++, 0);

            for (int i = 0; i < prefixMassList.Count; i++)
            {
                for (int j = i + 1; j < prefixMassList.Count; j++)
                {
                    int newMass = prefixMassList[j] - prefixMassList[i];
                    linearSpectrum.Insert(index++, newMass);
                }
            }

            linearSpectrum.Sort();

            return linearSpectrum;
        }

        public static List<string> PeptideEncoding(string DNA, string peptide)
        {
            List<string> output = new List<string>();
            string complement = "";
            string subString = "";
            string complementRNA = "";
            string subStringRNA = "";
            int offset = peptide.Length;
            string acid1 = "";
            string acid2 = "";

            for (int i = 0; i < DNA.Length; i++)
            {
                if ((i + (3 * offset) - 1) < DNA.Length)
                {
                    subString = DNA.Substring(i, 3 * offset);
                    complement = ReverseComplement(subString);

                    subStringRNA = EncodeDNAtoRNA(subString);
                    complementRNA = EncodeDNAtoRNA(complement);

                    acid1 = "";
                    acid2 = "";

                    for (int x = 0; x < subString.Length; x += 3)
                    {
                        acid1 += ProteinTranslation(subStringRNA.Substring(x, 3));
                        acid2 += ProteinTranslation(complementRNA.Substring(x, 3));
                    }

                    if (acid1 == peptide)
                        output.Add(subString);
                    else if (acid2 == peptide)
                        output.Add(subString);
                }
            }

            return output;
        }

        public static string EncodeDNAtoRNA(string DNA)
        {
            string RNA = DNA.Replace('T', 'U');

            return RNA;
        }

        public static string ProteinTranslation(string RNA)
        {
            string output = "";
            string aminoacid = "";
            string key = "";

            for (int i = 0; i < RNA.Length; i += 3)
            {
                key = RNA.Substring(i, 3);
                aminoacid = condonTable[key];

                if (aminoacid != "STOP")
                    output += aminoacid;
            }
            return output;
        }

        public static List<string> ContigGeneration(List<string> kmers)
        {
            List<string> result = new List<string>();
            List<string> directedGraph = new List<string>();
            List<string> contigList = new List<string>();

            directedGraph = DeBruijnGraph(kmers);

            contigList = MaximalNonBranchingPaths(directedGraph);

            for (int i = 0; i < contigList.Count; i++)
            {
                string[] strJoin = contigList[i].Split(" -> ");
                string contig = strJoin[0];

                for (int x = 1; x < strJoin.Length; x++)
                {
                    contig += strJoin[x].Substring(strJoin[x].Length - 1, 1);
                }
                result.Add(contig);
            }

            return result;
        }

        public static List<string> MaximalNonBranchingPaths(List<string> directedGraph)
        {
            List<string> output = new List<string>();
            Dictionary<string, Node> dict = new Dictionary<string, Node>();
            Node startingNode = new Node();

            foreach (string str in directedGraph)
            {
                string[] nodes = str.Split(" -> ");

                if (nodes.Length == 2)
                {
                    string[] connectedNodes = nodes[1].Split(",");
                    string key = nodes[0];
                    Node child = new Node(key);

                    foreach (string childNode in connectedNodes)
                    {
                        child.connectedNodes.Add(childNode);
                    }
                    dict.Add(key, child);
                }
            }

            // now Determine in & out degrees of each node
            foreach (Node nd in dict.Values)
            {
                nd.outDegrees += nd.connectedNodes.Count;

                foreach (string sn in nd.connectedNodes)
                {
                    if (dict.ContainsKey(sn))
                        dict[sn].inDegrees += 1;
                }
            }

            foreach (Node nd in dict.Values)
            {
                if (nd.inDegrees == 1 && nd.outDegrees == 1)
                {
                    int dummy = 33;
                }
                else    // non 1-in-1-out nodes
                {
                    if (nd.outDegrees > 0)
                    {
                        string key = nd.node;
                        string nonBranchPath;
                        nd.included = true;

                        foreach (string edge in nd.connectedNodes)
                        {
                            nonBranchPath = key;

                            if (dict.ContainsKey(edge))
                            {
                                Node neighbor = dict[edge];
                                while (neighbor.inDegrees == 1 && neighbor.outDegrees == 1)
                                {
                                    nonBranchPath += " -> " + neighbor.node;
                                    neighbor.included = true;

                                    if (neighbor.connectedNodes.Count == 1)
                                    {
                                        if (dict.ContainsKey(neighbor.connectedNodes[0]))
                                            neighbor = dict[neighbor.connectedNodes[0]];
                                        else
                                        {
                                            nonBranchPath += " -> " + neighbor.connectedNodes[0];
                                            break;
                                        }
                                    }
                                }
                                if (neighbor.inDegrees + neighbor.outDegrees > 2)       // end node
                                {
                                    nonBranchPath += " -> " + neighbor.node;
                                    neighbor.included = true;
                                }
                            }
                            else
                                nonBranchPath += " -> " + edge;

                            output.Add(nonBranchPath);
                        }
                    }
                }
            }

            // now loop to find isolated cycles
            foreach (Node nd in dict.Values)
            {
                if (!nd.included)
                {
                    if (nd.inDegrees == 1 && nd.outDegrees == 1)
                    {
                        // check for isolated cycle where all nodes have in/out degrees of 1
                        Dictionary<string, Node> isolatedCycle = new Dictionary<string, Node>();

                        string key = nd.node;
                        string cycle = "";
                        bool allNodes1 = true;
                        nd.included = true;

                        foreach (string edge in nd.connectedNodes)
                        {
                            cycle = key;

                            if (dict.ContainsKey(edge))
                            {
                                Node neighbor = dict[edge];

                                if (neighbor.inDegrees != 1 || neighbor.outDegrees != 1)
                                    allNodes1 = false;

                                while (neighbor.inDegrees == 1 && neighbor.outDegrees == 1)
                                {
                                    cycle += " -> " + neighbor.node;
                                    neighbor.included = true;

                                    if (neighbor.node == key)   // complete cycle
                                        break;

                                    if (neighbor.connectedNodes.Count == 1)
                                    {
                                        if (dict.ContainsKey(neighbor.connectedNodes[0]))
                                        {
                                            neighbor = dict[neighbor.connectedNodes[0]];
                                            if (neighbor.inDegrees != 1 || neighbor.outDegrees != 1)
                                                allNodes1 = false;
                                        }
                                        else
                                        {
                                            cycle += " -> " + neighbor.connectedNodes[0];
                                            break;
                                        }
                                    }
                                }
                            }

                            if (allNodes1)
                            {
                                output.Add(cycle);
                            }
                        }
                    }
                }
            }

            return output;
        }


        public static string StringSpelledByGappedPatterns(List<string> GappedPatterns, int k, int d)
        {
            string result = "";
            List<string> pattern1 = new List<string>();
            List<string> pattern2 = new List<string>();
            string[] split;
            string prefixPattern = "";
            string suffixPattern = "";
            bool perfectMatch = true;

            foreach (string str in GappedPatterns)
            {
                split = str.Split("|");
                pattern1.Add(split[0]);
                pattern2.Add(split[1]);
            }

            prefixPattern = StringSpelledByGenomePath(pattern1);
            suffixPattern = StringSpelledByGenomePath(pattern2);

            for (int i = k + d; i < prefixPattern.Length; i++)
            {
                string p1 = prefixPattern[i].ToString();
                string p2 = suffixPattern[i - k - d].ToString();
                if (prefixPattern[i] != suffixPattern[i - k - d])
                {
                    perfectMatch = false;
                    break;
                }
            }

            if (perfectMatch)
            {
                string p1 = suffixPattern.Substring(suffixPattern.Length - (k + d), k + d);
                result = prefixPattern + suffixPattern.Substring(suffixPattern.Length - (k + d), k + d);
            }
            else
                result = "there is no string spelled by the gapped patterns";

            return result;
        }

        static public List<string> ConvertToReadPairs(string sequence)
        {
            List<string> readPairList = new List<string>();
            string[] kmers;
            string prefix = "";
            string suffix = "";

            kmers = sequence.Split("->");
            int k = kmers[0].Length / 2;    // two nodes included

            for (int i = 1; i < kmers.Length; i++)
            {
                prefix = kmers[i - 1].Substring(0, k) + kmers[i].Substring(k - 1, 1);
                suffix = kmers[i - 1].Substring(k, k) + kmers[i].Substring(kmers[i].Length - 1, 1);

                readPairList.Add(prefix + "|" + suffix);
            }

            return readPairList;
        }

        static public string StringReconstructionFromReadPairs(int k, int d, List<string> GappedPatterns)
        {
            string sequence = "";
            string result = "";
            List<string> directedGraph = new List<string>();
            List<string> orderedKmers = new List<string>();

            directedGraph = PairedDeBruijnGraph(GappedPatterns);

            sequence = EulerianCycle(directedGraph);

            orderedKmers = ConvertToReadPairs(sequence);

            result = StringSpelledByGappedPatterns(orderedKmers, k, d);


            return result;
        }

        static public string StringReconstruction(int k, List<string> kmers)
        {
            string sequence = "";
            List<string> directedGraph = new List<string>();
            List<string> orderedKmers = new List<string>();

            directedGraph = DeBruijnGraph(kmers);

            sequence = EulerianCycle(directedGraph);

            string[] tmpKmers = sequence.Split("->");

            foreach (string s in tmpKmers)
                orderedKmers.Add(s);

            sequence = StringSpelledByGenomePath(orderedKmers);

            return sequence;
        }


        static public string UniversalCircularString(int k)
        {
            string result = "";
            List<string> kmers = new List<string>();

            kmers = UniversalGenerator(k);

            kmers.Sort();

            result = StringReconstruction(k, kmers);

            string prefix = result.Substring(0, k - 1);
            string suffix = result.Substring(result.Length - k + 1, k - 1);

            if (prefix == suffix)
            {
                result = result.Substring(0, result.Length - k + 1);
            }

            return result;
        }

        static public List<string> UniversalGenerator(int kmerSize)
        {
            List<string> outputList = new List<string>();

            for (int i = 0; i < kmerSize; i++)
            {
                outputList = GenerateTreeNodes(zero_ones, outputList);
            }

            outputList.Sort();

            return outputList;
        }

        static public List<string> GenerateTreeNodes(char[] pattern, List<string> inputList)
        {
            List<string> outputList = new List<string>();

            if (inputList.Count > 0)
            {
                foreach (string str in inputList)
                {
                    foreach (char c in pattern)
                    {
                        outputList.Add(str + c.ToString());
                    }
                }
            }
            else
            {
                foreach (char c in pattern)
                {
                    outputList.Add(c.ToString());
                }
            }

            return outputList;
        }

        static public string EulerianCycle(List<string> directedGraph)
        {
            string cycleOutput = "";
            Dictionary<string, Node> dict = new Dictionary<string, Node>();
            Node startingNode = new Node();
            string startKey = "";
            Random random = new Random();
            int ranIndex = random.Next(0, directedGraph.Count - 1);
            int count = 0;

            foreach (string str in directedGraph)
            {
                string[] nodes = str.Split(" -> ");

                if (nodes.Length == 2)
                {
                    string[] connectedNodes = nodes[1].Split(",");
                    string key = nodes[0];
                    Node child = new Node(key);

                    foreach (string childNode in connectedNodes)
                    {
                        child.connectedNodes.Add(childNode);
                    }
                    dict.Add(key, child);

                    // hack to get random starting node
                    if (ranIndex == count)
                        startKey = key;

                    count++;
                }
            }

            // now Determine in & out degrees of each node
            foreach (Node nd in dict.Values)
            {
                nd.outDegrees += nd.connectedNodes.Count;

                foreach (string sn in nd.connectedNodes)
                {
                    if (dict.ContainsKey(sn))
                        dict[sn].inDegrees += 1;
                }
            }

            startingNode = dict[startKey];
            foreach (Node nd in dict.Values)
            {
                if (nd.inDegrees != nd.outDegrees)
                {
                    if (nd.outDegrees == nd.inDegrees + 1)
                    {
                        startingNode = nd;
                        break;
                    }
                }
            }

            Circuit circuit = new Circuit();

            FindCircuit(dict, startingNode, circuit);

            for (int k = circuit.circuit.Count - 1; k >= 0; k--)
            {
                cycleOutput += circuit.circuit[k] + "->";
            }
            // remove last arrow
            cycleOutput = cycleOutput.Substring(0, cycleOutput.Length - 2);


            return cycleOutput;
        }

        public static void FindCircuit(Dictionary<string, Node> nodeDict, Node currentNode, Circuit circuit)
        {
            if (currentNode.connectedNodes.Count == 0)
            {
                circuit.AddNode(currentNode);
            }
            else
            {
                while (currentNode.connectedNodes.Count > 0)
                {
                    string nextNode = currentNode.connectedNodes[0];  // pick next one 

                    currentNode.connectedNodes.RemoveAt(0);

                    if (nodeDict.ContainsKey(nextNode))
                    {
                        Node neighbor = nodeDict[nextNode];

                        FindCircuit(nodeDict, neighbor, circuit);
                    }
                    else // end of Eulerian path
                    {
                        Node endNode = new Node(nextNode);
                        circuit.AddNode(endNode);
                    }
                }
                circuit.AddNode(currentNode);
            }
        }

        public static List<string> DeBruijnGraph(List<string> kmers)
        {
            List<string> output = new List<string>();
            Dictionary<string, List<string>> dict = new Dictionary<string, List<string>>();
            string prefix = "";
            int k = kmers[0].Length;
            string currentNode = "";
            string connectedNode = "";
            string suffix = "";
            bool found = false;
            bool nodeFound = false;

            for (int i = 0; i < kmers.Count; i++)
            {
                currentNode = kmers[i];
                suffix = currentNode.Substring(1, k - 1);
                found = false;
                List<string> Nodes = new List<string>();

                for (int x = 0; x < kmers.Count; x++)
                {
                    connectedNode = kmers[x];
                    prefix = connectedNode.Substring(0, k - 1);

                    nodeFound = false;

                    if (suffix == prefix)
                    {
                        // check if already in list
                        foreach (string s in Nodes)
                        {
                            if (s == suffix)
                            {
                                nodeFound = true;
                                break;
                            }
                        }

                        if (!nodeFound)
                        {
                            Nodes.Add(suffix);
                            found = true;
                        }
                    }
                    else if (!nodeFound && x == kmers.Count - 1)    // scanned all nodes, check if not found, if so at the last node, so add
                    {
                        if (Nodes.Count == 0)
                        {
                            Nodes.Add(suffix);
                            found = true;
                        }
                    }
                }

                if (found)
                {
                    Nodes.Sort();
                    string key = currentNode.Substring(0, k - 1);

                    if (!dict.ContainsKey(key))
                    {
                        dict.Add(key, Nodes);
                    }
                    else
                    {
                        List<string> tmpNodes = new List<string>();

                        tmpNodes = dict[key];

                        foreach (string s in Nodes)
                            tmpNodes.Add(s);

                        tmpNodes.Sort();

                        dict[key] = tmpNodes;
                    }
                }
            }

            string outputNode = "";

            foreach (string key in dict.Keys)
            {
                List<string> row = dict[key];

                outputNode = key + " -> ";

                foreach (string item in row)
                    outputNode += item + ",";

                outputNode = outputNode.Substring(0, outputNode.Length - 1);        // remove last comma

                output.Add(outputNode);
            }
            output.Sort();

            WriteListToFile("C:\\Temp\\output.txt", output);

            return output;
        }

        public static List<string> PairedDeBruijnGraph(List<string> GappedPatterns)
        {
            List<string> output = new List<string>();
            Dictionary<string, List<string>> dict = new Dictionary<string, List<string>>();
            List<string> kmers = new List<string>();
            List<string> pairedKmers = new List<string>();
            string prefix = "";
            string currentNode = "";
            string connectedNode = "";
            string suffix = "";
            bool found = false;
            bool nodeFound = false;
            string pairedPrefix = "";
            string pairedSuffix = "";
            string combinedSuffix = "";
            string combinedPrefix = "";
            string[] split;

            foreach (string str in GappedPatterns)
            {
                split = str.Split("|");
                kmers.Add(split[0]);
                pairedKmers.Add(split[1]);
            }

            int k = kmers[0].Length;

            for (int i = 0; i < kmers.Count; i++)
            {
                currentNode = kmers[i];
                suffix = currentNode.Substring(1, k - 1);
                pairedSuffix = pairedKmers[i].Substring(1, k - 1);

                combinedSuffix = suffix + pairedSuffix;

                found = false;
                List<string> Nodes = new List<string>();

                for (int x = 0; x < kmers.Count; x++)
                {
                    connectedNode = kmers[x];
                    prefix = connectedNode.Substring(0, k - 1);
                    pairedPrefix = pairedKmers[x].Substring(0, k - 1);

                    combinedPrefix = prefix + pairedPrefix;

                    nodeFound = false;

                    if (combinedSuffix == combinedPrefix)
                    {
                        // check if already in list
                        foreach (string s in Nodes)
                        {
                            if (s == combinedSuffix)
                            {
                                nodeFound = true;
                                break;
                            }
                        }

                        if (!nodeFound)
                        {
                            Nodes.Add(combinedSuffix);
                            found = true;
                        }
                    }
                    else if (!nodeFound && x == kmers.Count - 1)    // scanned all nodes, check if not found, if so at the last node, so add
                    {
                        if (Nodes.Count == 0)
                        {
                            Nodes.Add(combinedSuffix);
                            found = true;
                        }
                    }
                }

                if (found)
                {
                    Nodes.Sort();

                    string key = kmers[i].Substring(0, k - 1) + pairedKmers[i].Substring(0, k - 1);

                    if (!dict.ContainsKey(key))
                    {
                        dict.Add(key, Nodes);
                    }
                    else
                    {
                        List<string> tmpEdges = new List<string>();

                        tmpEdges = dict[key];

                        foreach (string s in Nodes)
                            tmpEdges.Add(s);

                        tmpEdges.Sort();

                        dict[key] = tmpEdges;
                    }
                }
            }

            string outputNode = "";

            foreach (string key in dict.Keys)
            {
                List<string> row = dict[key];

                outputNode = key + " -> ";

                foreach (string item in row)
                {
                    outputNode += item + ",";
                }

                outputNode = outputNode.Substring(0, outputNode.Length - 1);        // remove last comma

                output.Add(outputNode);
            }
            output.Sort();

            WriteListToFile("C:\\Temp\\output.txt", output);

            return output;
        }

        public static List<string> DeBruijnGraph(int k, string Text)
        {
            List<string> output = new List<string>();
            List<string> kmers = new List<string>();

            kmers = Composition(Text, k);

            output = DeBruijnGraph(kmers);

            return output;
        }

        static public int NumberOfOccurances(string kmer, List<string> kmerList)
        {
            int count = 0;

            foreach (string str in kmerList)
            {
                if (str == kmer)
                    count += 1;
            }
            return count;
        }

        static public List<string> OverlapGraph(List<string> kmers)
        {
            List<string> output = new List<string>();
            Dictionary<string, List<string>> dict = new Dictionary<string, List<string>>();
            string prefix = "";
            int k = kmers[0].Length;
            string currentNode = "";
            string connectedNode = "";
            string suffix = "";
            bool found = false;
            bool nodeFound = false;

            for (int i = 0; i < kmers.Count; i++)
            {
                currentNode = kmers[i];
                suffix = currentNode.Substring(1, k - 1);
                found = false;
                List<string> Nodes = new List<string>();

                for (int x = 0; x < kmers.Count; x++)
                {
                    if (x != i)
                    {
                        connectedNode = kmers[x];
                        prefix = connectedNode.Substring(0, k - 1);

                        if (suffix == prefix)
                        {
                            // check if already in list
                            nodeFound = false;
                            foreach (string s in Nodes)
                            {
                                if (s == connectedNode)
                                {
                                    nodeFound = true;
                                    break;
                                }
                            }

                            if (!nodeFound)
                            {
                                Nodes.Add(connectedNode);
                                found = true;
                            }
                        }
                    }
                }

                if (found)
                {
                    Nodes.Sort();
                    if (!dict.ContainsKey(currentNode))
                        dict.Add(currentNode, Nodes);
                    else
                        Console.WriteLine("Key exists");
                }
            }

            string outputNode = "";

            foreach (string key in dict.Keys)
            {
                List<string> row = dict[key];

                outputNode = key + " -> ";

                foreach (string item in row)

                    outputNode += item + ",";

                outputNode = outputNode.Substring(0, outputNode.Length - 1);        // remove last comma

                output.Add(outputNode);
            }

            WriteListToFile("C:\\Temp\\output.txt", output);


            return output;
        }


        static public string StringSpelledByGenomePathOld(List<string> kmers)
        {
            string sequence = "";
            string prefix = "";
            int k = kmers[0].Length;
            string prevVal = "";

            sequence = kmers[0];

            for (int i = 1; i < kmers.Count; i++)
            {
                prefix = kmers[i].Substring(0, k - 1);

                prevVal = sequence.Substring(sequence.Length - k + 1, k - 1);

                if (prefix == prevVal)
                {
                    sequence += kmers[i].Substring(k - 1, 1);
                }
            }

            return sequence;
        }

        static public string StringSpelledByGenomePath(List<string> kmers)
        {
            string sequence = "";
            string prefix = "";
            string suffix = "";
            int k = kmers[0].Length;

            sequence = kmers[0];

            for (int i = 1; i < kmers.Count; i++)
            {
                prefix = kmers[i].Substring(0, k - 1);
                suffix = kmers[i - 1].Substring(1, k - 1);

                if (prefix == suffix)
                {
                    sequence += kmers[i].Substring(kmers[i].Length - 1, 1);       // add last character of current kmer
                }
            }

            return sequence;
        }


        static public string Reconstruction(int k, List<string> kmerArray)
        {
            string sequence = "";
            string suffix = "";
            string prefix = "";

            foreach (string kmer in kmerArray)
            {
                suffix = kmer.Substring(1, kmer.Length - 1);

                for (int i = 0; i < kmerArray.Count; i++)
                {
                    prefix = kmerArray[i].Substring(0, kmer.Length - 1);
                    if (prefix == suffix)
                    {
                        break;
                    }
                }
                sequence += prefix;
            }

            return sequence;
        }

        static public List<string> Composition(string Text, int k)
        {
            List<string> kmerArray = new List<string>();
            string str = "";

            for (int i = 0; i < Text.Length - k + 1; i++)
            {
                str = Text.Substring(i, k);
                kmerArray.Add(Text.Substring(i, k));
            }

            return kmerArray;
        }

        //===================================================================================================
        //  RandomizedMotifSearch(Dna, k, t)
        //      randomly select k-mers Motifs = (Motif1, …, Motift) in each string from Dna
        //      BestMotifs ← Motifs
        //      while forever
        //          Profile ← Profile(Motifs)
        //          Motifs ← Motifs(Profile, Dna)
        //          if Score(Motifs) < Score(BestMotifs)
        //              BestMotifs ← Motifs
        //          else
        //              return BestMotifs
        //===================================================================================================
        static public List<string> RandominzedMotifSearch(string[] Dna, int k, int t, int N = 1000)
        {
            List<string> BestMotifs = new List<string>();
            List<string> tmpMotifs = new List<string>();
            List<string> lowestMotifs = new List<string>();
            int tmpScore = 0;
            int bestScore = int.MaxValue;
            int lowestScore = int.MaxValue;

            Random random = new Random();

            for (int a = 0; a < N; a++)
            {
                for (int i = 0; i < Dna.Length; i++)
                {
                    int length = Dna[i].Length;
                    int ran = random.Next(0, length - k);
                    string kmer = Dna[i].Substring(ran, k);

                    BestMotifs.Add(kmer);
                    tmpMotifs.Add(kmer);
                }

                while (true)
                {
                    double[,] profileMatrix = Profile(tmpMotifs, false);

                    tmpMotifs = Motif(profileMatrix, Dna, k, true);
                    tmpScore = Score(tmpMotifs);
                    bestScore = Score(BestMotifs);

                    if (tmpScore < bestScore)
                    {
                        BestMotifs.Clear();
                        foreach (string s in tmpMotifs)
                            BestMotifs.Add(s);

                        if (tmpScore < lowestScore)
                        {
                            lowestScore = tmpScore;
                            lowestMotifs.Clear();
                            foreach (string s in tmpMotifs)
                                lowestMotifs.Add(s);
                        }
                    }
                    else
                        break;
                }

                BestMotifs.Clear();
                foreach (string s in lowestMotifs)
                    BestMotifs.Add(s);
            }

            string tmp = "";

            foreach (string s in BestMotifs)
                tmp += s + "  ";

            return BestMotifs;
        }

        static public int Random(double[] array)
        {
            int index = array.Length - 1;     // Set to last index by default
            double total = 0;
            double[] distributionArray = new double[array.Length];
            Random random = new Random();

            foreach (double d in array)
                total += d;

            for (int i = 0; i < distributionArray.Length; i++)
            {
                if (i == 0)
                    distributionArray[i] = (array[i] / total);
                else
                    distributionArray[i] = ((array[i] / total) + distributionArray[i - 1]);
            }

            double dRan = random.NextDouble();
            for (int a = 0; a < distributionArray.Length; a++)
            {
                if (dRan <= distributionArray[a])
                {
                    index = a;
                    break;
                }
            }

            return index;
        }

        static public List<string> GibbsSampler(string[] Dna, int k, int t, int N)
        {
            List<string> BestMotifs = new List<string>();
            List<string> tmpMotifs = new List<string>();
            List<string> lowestMotifs = new List<string>();
            List<string> seedMotifs = new List<string>();

            string selectedMotif = "";
            int tmpScore = 0;
            int ranMotif = 0;
            int bestScore = int.MaxValue;
            int lowestScore = int.MaxValue;

            Random random = new Random();

            seedMotifs = RandominzedMotifSearch(Dna, k, t, 30);

            for (int a = 0; a < N; a++)
            {
                for (int b = 0; b < 100; b++)
                {
                    if (a > 0)
                    {
                        seedMotifs.Clear();

                        for (int c = 0; c < Dna.Length; c++)
                        {
                            int length = Dna[c].Length;
                            int ran = random.Next(0, length - k);
                            string kmer = Dna[c].Substring(ran, k);

                            seedMotifs.Add(kmer);
                        }
                    }

                    if (b == 0)
                    {
                        BestMotifs.Clear();
                        foreach (string s in seedMotifs)
                            BestMotifs.Add(s);
                    }

                    if (Score(seedMotifs) < Score(BestMotifs))
                    {
                        BestMotifs.Clear();
                        foreach (string s in seedMotifs)
                            BestMotifs.Add(s);
                    }
                }

                if (a == 0)
                {
                    tmpMotifs.Clear();
                    foreach (string s in BestMotifs)
                        tmpMotifs.Add(s);

                }

                while (true)
                {
                    ranMotif = random.Next(0, t - 1);

                    selectedMotif = tmpMotifs[ranMotif];

                    tmpMotifs.RemoveAt(ranMotif);

                    double[,] profileMatrix = Profile(tmpMotifs, true);

                    selectedMotif = ProfileRandom(Dna[ranMotif], k, profileMatrix);

                    tmpMotifs.Insert(ranMotif, selectedMotif);

                    tmpScore = Score(tmpMotifs);
                    bestScore = Score(BestMotifs);

                    if (tmpScore < bestScore)
                    {
                        BestMotifs.Clear();
                        foreach (string s in tmpMotifs)
                            BestMotifs.Add(s);

                        if (tmpScore < lowestScore)
                        {
                            lowestScore = tmpScore;
                            lowestMotifs.Clear();
                            foreach (string s in tmpMotifs)
                                lowestMotifs.Add(s);
                        }
                    }
                    else
                        break;
                }

                if (lowestMotifs.Count > 0)
                {
                    BestMotifs.Clear();
                    foreach (string s in lowestMotifs)
                        BestMotifs.Add(s);
                }
            }


            WriteListToFile("C:\\Temp\\Output.txt", BestMotifs);

            string tmp = "";

            foreach (string s in BestMotifs)
                tmp += s + "  ";

            return BestMotifs;
        }

        static public List<string> GreedyMotifSearch(string[] Dna, int k, int t)
        {
            List<string> BestMotifs = new List<string>();
            List<string> tmpMotifs = new List<string>();
            string kmer = "";
            string nextMotif = "";
            string strand = "";
            string baseStrand = "";

            for (int i = 0; i < Dna.Length; i++)
                BestMotifs.Add(Dna[i].Substring(0, k));

            baseStrand = Dna[0];

            for (int m = 0; m < baseStrand.Length - k + 1; m++)
            {
                tmpMotifs.Clear();
                kmer = baseStrand.Substring(m, k);
                tmpMotifs.Add(kmer);

                for (int n = 1; n < t; n++)
                {
                    double[,] profileMatrix = Profile(tmpMotifs, true);

                    strand = Dna[n];

                    nextMotif = ProfileMostProbableKmer(strand, k, profileMatrix);

                    tmpMotifs.Add(nextMotif);
                }

                if (Score(tmpMotifs) < Score(BestMotifs))
                {
                    BestMotifs.Clear();
                    foreach (string s in tmpMotifs)
                        BestMotifs.Add(s);
                }
            }

            string tmp = "";

            foreach (string s in BestMotifs)
                tmp += s + "  ";

            return BestMotifs;
        }

        static public int Score(List<string> Motifs)
        {
            int motifLength = 1;
            int t = 0;
            int maxCount = 0;
            int totalScore = 0;
            string consensusMotif = "";
            char[] Nucleotides = { 'A', 'C', 'G', 'T' };

            if (Motifs.Count > 0)
            {
                motifLength = Motifs[0].Length;
                t = Motifs.Count;
            }

            int[,] countArray = new int[NucleotideSize, motifLength];
            int[] scoreArray = new int[motifLength];
            //            string str = "";
            //            string c = "";
            string currentMotif = "";
            char currentChar = ' ';

            for (int i = 0; i < motifLength; i++)
            {
                for (int x = 0; x < t; x++)     // row of Motifs
                {
                    currentMotif = Motifs[x];
                    currentChar = currentMotif[i];

                    for (int m = 0; m < Nucleotides.Length; m++)
                    {
                        if (currentChar == Nucleotides[m])
                        {
                            countArray[m, i] += 1;
                            break;
                        }
                    }
                }

                maxCount = 0;
                int index = 0;
                for (int b = 0; b < Nucleotides.Length; b++)
                {
                    if (countArray[b, i] > maxCount)
                    {
                        maxCount = countArray[b, i];
                        index = b;
                    }
                }
                currentChar = Nucleotides[index];
                consensusMotif += currentChar.ToString();
            }

            //            Console.WriteLine("Score - Consensu Mofit = {0}", consensusMotif);

            totalScore = 0;

            foreach (string s in Motifs)
                totalScore += HammingDistance(consensusMotif, s);

            return totalScore;
        }

        static public double[,] Profile(List<string> Motifs, bool bIncludePseudocounts = false)
        {
            int motifLength = 1;
            int t = 0;
            int denominator = 0;

            if (Motifs.Count > 0)
            {
                motifLength = Motifs[0].Length;
                t = Motifs.Count;
                if (bIncludePseudocounts)
                    denominator = t * 2;        // Strange Denominator effect
                else
                    denominator = t;
            }

            double[,] profileArray = new double[NucleotideSize, motifLength];
            double[,] countArray = new double[NucleotideSize, motifLength];
            string str = "";
            string c = "";

            if (bIncludePseudocounts == true)
            {
                for (int a = 0; a < NucleotideSize; a++)
                    for (int b = 0; b < motifLength; b++)
                        countArray[a, b] = 1;
            }


            for (int i = 0; i < Motifs.Count; i++)
            {
                str = Motifs[i];
                for (int m = 0; m < motifLength; m++)
                {
                    c = str.Substring(m, 1);

                    if (c == "A")
                    {
                        countArray[0, m] += 1;
                    }
                    else if (c == "C")
                    {
                        countArray[1, m] += 1;
                    }
                    else if (c == "G")
                    {
                        countArray[2, m] += 1;
                    }
                    else if (c == "T")
                    {
                        countArray[3, m] += 1;
                    }
                }
            }

            for (int x = 0; x < NucleotideSize; x++)
                for (int y = 0; y < motifLength; y++)
                    profileArray[x, y] = countArray[x, y] / denominator;

            return profileArray;
        }

        static public List<string> Motif(double[,] profile, string[] Dna, int k, bool bRandom = false)
        {
            List<string> motifs = new List<string>();
            string mostCommon = "";

            foreach (string str in Dna)
            {
                mostCommon = ProfileMostProbableKmer(str, k, profile, bRandom);
                motifs.Add(mostCommon);
            }

            return motifs;
        }

        static public string ProfileMostProbableKmer(string Text, int k, double[,] matrix, bool bUseRandom = false)
        {
            string result = "";
            string pattern = "";
            string c = "";
            double probability = 1;
            double highestProbabilty = -1;
            double[] probArray = new double[Text.Length - k + 1];

            for (int i = 0; i < Text.Length - k + 1; i++)
            {
                pattern = Text.Substring(i, k);
                probability = 1;
                for (int m = 0; m < pattern.Length; m++)
                {
                    c = pattern.Substring(m, 1);

                    if (c == "A")
                        probability *= matrix[0, m];
                    else if (c == "C")
                        probability *= matrix[1, m];
                    else if (c == "G")
                        probability *= matrix[2, m];
                    else if (c == "T")
                        probability *= matrix[3, m];
                }

                probArray[i] = probability;

                if (probability > highestProbabilty)
                {
                    highestProbabilty = probability;
                    result = pattern;
                }
            }

            if (bUseRandom)
            {
                int ranIndex = Random(probArray);
                result = Text.Substring(ranIndex, k);
            }

            //            Console.WriteLine("ProfileMostProbableKmer:  {0}", result);
            return result;
        }

        static public string ProfileRandom(string Text, int k, double[,] matrix)
        {
            string result = "";
            string origPattern = "";
            int origKmerIndex = 0;
            string pattern = "";
            string c = "";
            double probability = 1;
            double highestProbabilty = -1;
            double[] probArray = new double[Text.Length - k + 1];

            for (int i = 0; i < Text.Length - k + 1; i++)
            {
                pattern = Text.Substring(i, k);
                probability = 1;
                for (int m = 0; m < pattern.Length; m++)
                {
                    c = pattern.Substring(m, 1);

                    if (c == "A")
                        probability *= matrix[0, m];
                    else if (c == "C")
                        probability *= matrix[1, m];
                    else if (c == "G")
                        probability *= matrix[2, m];
                    else if (c == "T")
                        probability *= matrix[3, m];
                }

                probArray[i] = probability;

                if (probability > highestProbabilty)
                {
                    highestProbabilty = probability;
                    origPattern = pattern;
                    origKmerIndex = i;
                }
            }

            int ranIndex = Random(probArray);
            result = Text.Substring(ranIndex, k);

            //            Console.WriteLine("ProfileMostProbableKmer:  {0}", result);
            return result;
        }

        // MedianString(Dna, k)
        //      distance ← ∞
        //      for i ←0 to 4k −1
        //          Pattern ← NumberToPattern(i, k)
        //          if distance > DistanceBetweenPatternAndStrings(Pattern, Dna)
        //              distance ← DistanceBetweenPatternAndStrings(Pattern, Dna)
        //              Median ← Pattern
        //      return Median
        static public string MedianString(string[] Dna, int k)
        {
            string Median = "";
            int distance = int.MaxValue;
            Int64 arraySize = (Int64)Math.Pow(4, k);
            string pattern = "";
            int tmpDistance = 0;

            for (Int64 i = 0; i < arraySize - 1; i++)
            {
                pattern = NumberToPattern2(i, k);
                tmpDistance = DistanceBetweenPatternAndStrings(pattern, Dna);
                if (distance > tmpDistance)
                {
                    distance = tmpDistance;
                    Median = pattern;
                }
            }

            Console.WriteLine("MedianString:  {0}", Median);
            return Median;
        }

        //===================================================================================================
        //        DistanceBetweenPatternAndStrings(Pattern, Dna)
        //              k ← |Pattern|
        //              distance ← 0
        //              for each string Text in Dna
        //                  HammingDistance ← ∞
        //                  for each k-mer Pattern’ in Text
        //                      if HammingDistance > HammingDistance(Pattern, Pattern’)
        //                          HammingDistance ← HammingDistance(Pattern, Pattern’)
        //                  distance ← distance + HammingDistance
        //              return distance
        //===================================================================================================
        static public int DistanceBetweenPatternAndStrings(string Pattern, string[] Dna)
        {
            int distance = 0;
            int k = Pattern.Length;
            string str = "";
            int tmpDistance = 0;

            foreach (string Text in Dna)
            {
                int HamDistance = int.MaxValue;

                for (int m = 0; m < Text.Length - Pattern.Length + 1; m++)
                {
                    str = Text.Substring(m, k);
                    tmpDistance = HammingDistance(Pattern, str);
                    if (HamDistance > tmpDistance)
                        HamDistance = tmpDistance;
                }
                distance += HamDistance;
            }

            //            Console.WriteLine("DistanceBetweenPatternAndStrings:  {0}", distance);
            return distance;
        }

        //===================================================================================================
        //        Patterns ← an empty set
        //        
        //        for each k-mer Pattern in Dna
        //            for each k-mer Pattern’ differing from Pattern by at most d mismatches
        //               if Pattern' appears in each string from Dna with at most d ﻿mismatches
        //                    add Pattern' to Patterns
        //        
        //    remove duplicates from Patterns
        //===================================================================================================
        static public string MotifEnumeration(string[] Dna, int k, int d)
        {
            string Patterns = "";
            List<string> motifList = new List<string>();

            foreach (string kmer in Dna)
            {
                for (int m = 0; m < kmer.Length - k + 1; m++)        // All kmers in DNA
                {
                    string str = kmer.Substring(m, k);
                    List<string> neighbors = new List<string>();    // All kmers with d mismatches

                    neighbors = Neighbors(str, d);

                    foreach (string n in neighbors)
                    {
                        string strMatch = "";
                        bool bInDNA = true;

                        foreach (string inner in Dna)
                        {
                            strMatch = ApproximatePatternMatching(n, inner, d);  // Check each DNA for pattern with d mismatch
                            if (strMatch.Length == 0)
                            {
                                bInDNA = false;
                                break;
                            }
                        }
                        if (bInDNA) // str in all DNA
                            motifList.Add(n);
                    }
                }
            }

            motifList = RemoveDuplicates(motifList);

            foreach (string s in motifList)
                Patterns += s + " ";

            Console.WriteLine("MotifEnumeration:  {0}", Patterns);
            return Patterns;
        }

        static public List<string> Neighbors(string Pattern, int d)
        {
            List<string> Neighborhood = new List<string>();
            List<string> suffixNeighbors = new List<string>();
            string suffix;
            string prefix;
            char[] nucleotide = { 'A', 'C', 'T', 'G' };
            string kmer = "";

            if (0 == d)
            {
                Neighborhood.Add(Pattern);
                return Neighborhood;
            }

            if (Pattern.Length == 1)
            {
                Neighborhood.Add("A");
                Neighborhood.Add("C");
                Neighborhood.Add("G");
                Neighborhood.Add("T");

                return Neighborhood;
            }

            Neighborhood.Clear();

            suffix = Pattern.Substring(1, Pattern.Length - 1);
            prefix = Pattern[0].ToString();

            suffixNeighbors = Neighbors(suffix, d);
            foreach (string str in suffixNeighbors)
            {
                if (HammingDistance(suffix, str) < d)
                {
                    foreach (char c in nucleotide)
                    {
                        kmer = c.ToString() + str;
                        Neighborhood.Add(kmer);
                    }
                }
                else
                {
                    kmer = prefix + str;
                    Neighborhood.Add(kmer);
                }
            }

            //            Console.WriteLine("Neighbors:  {0}", Neighborhood);
            return Neighborhood;
        }

        static public string ImmediateNeighbors(string Pattern)
        {
            string Neighborhood = Pattern;
            char[] nucleotide = { 'A', 'C', 'T', 'G' };
            char symbol;
            string kmer = "";


            for (int i = 0; i < Pattern.Length; i++)
            {
                symbol = Pattern[i];
                foreach (char c in nucleotide)
                {
                    if (symbol != c)
                    {
                        kmer = ReplaceStringAt(Pattern, c, i);
                        Neighborhood += " " + kmer;

                    }
                }
            }

            Console.WriteLine("ImmediateNeighbors:  {0}", Neighborhood);
            return Neighborhood;
        }

        static public string ClumpFinder(string Genome, int k, int L, int t)
        {
            string result = "";
            string subPattern = "";
            List<string> FrequentPatterns = new List<string>();
            List<string> tmpPatterns = new List<string>();

            for (Int64 i = 0; i < Genome.Length - L + 1; i++)
            {
                subPattern = Genome.Substring((int)i, L);

                tmpPatterns = FrequentWords2(subPattern, k, t);

                foreach (string s in tmpPatterns)
                {
                    FrequentPatterns.Add(s);
                }
            }

            FrequentPatterns = RemoveDuplicates(FrequentPatterns);


            foreach (string str in FrequentPatterns)
                result += str + " ";

            Console.WriteLine("ClumpFinder:  {0}", result);
            return result;
        }

        static public int PatternCount(string Text, string Pattern)
        {
            int Count = 0;
            string substring = "";

            for (int i = 0; i < Text.Length - Pattern.Length + 1; i++)
            {
                substring = Text.Substring(i, Pattern.Length);
                if (substring == Pattern)
                    Count += 1;
            }

            //            Console.WriteLine("PatternCount:  {0} = {1}", Pattern, Count);

            return Count;
        }

        static public int[] ComputingFrequenciesWithMismatches(string Text, int k, int d, bool IncludeReverse = false)
        {
            Int64 arraySize = (Int64)Math.Pow(4, k);
            int[] frequncyArray = new int[arraySize];
            List<string> Neighborhood = new List<string>();
            string pattern = "";
            string complement = "";
            Int64 j = 0;

            for (int i = 0; i < Text.Length - k + 1; i++)
            {
                pattern = Text.Substring(i, k);
                Neighborhood = Neighbors(pattern, d);

                foreach (string p in Neighborhood)
                {
                    j = PatternToNumber2(p);
                    frequncyArray[j] += 1;
                }

                if (IncludeReverse)
                {
                    complement = ReverseComplement(pattern);

                    Neighborhood = Neighbors(complement, d);

                    foreach (string p in Neighborhood)
                    {
                        j = PatternToNumber2(p);
                        frequncyArray[j] += 1;
                    }
                }

            }

            // Casues exception, look into different implementation.
            //            for (Int64 i = 0; i < arraySize; i++)
            //                result += frequncyArray[i].ToString() + " ";

            //            Console.WriteLine("ComnputingFrequencieswithMatchs:  {0}", result);
            return frequncyArray;
        }

        static public string ComputingFrequencies(string Text, int k)
        {
            string result = "";
            Int64 arraySize = (Int64)Math.Pow(4, k);
            Int64[] frequncyArray = new Int64[arraySize];
            string pattern = "";
            Int64 j = 0;

            for (Int64 i = 0; i < arraySize; i++)
                frequncyArray[i] = 0;

            for (Int64 m = 0; m < Text.Length - k + 1; m++)
            {
                pattern = Text.Substring((int)m, k);
                j = PatternToNumber2(pattern);
                frequncyArray[j] += 1;
            }

            for (Int64 i = 0; i < arraySize; i++)
                result += frequncyArray[i].ToString() + " ";

            Console.WriteLine("ComputingFrequencies:  {0}", result);
            return result;
        }

        static public List<string> FrequentWords(string Text, int k)
        {
            List<string> FrequentPatterns = new List<string>();
            List<int> Count = new List<int>();
            List<string> Kmers = new List<string>();
            int MaxCount = 0;
            int localCount = 0;

            string pattern = "";

            for (int i = 0; i < Text.Length - k; i++)
            {
                pattern = Text.Substring(i, k);

                localCount = PatternCount(Text, pattern);
                if (localCount > MaxCount)
                    MaxCount = localCount;

                Count.Add(localCount);
                Kmers.Add(pattern);
            }

            List<string> tmpKmers = new List<string>();

            for (k = 0; k < Count.Count; k++)
            {
                if (Count[k] == MaxCount)
                    tmpKmers.Add(Kmers[k]);
            }

            FrequentPatterns = RemoveDuplicates(tmpKmers);

            /*            foreach (string s in tmpKmers)
                        {
                            FrequentPatterns += s + " ";
                        } */

            //            Console.WriteLine("FrequentWords:  {0}", FrequentPatterns);

            return FrequentPatterns;

        }

        static public List<string> FrequentWords2(string Text, int k, int t)
        {
            List<string> FrequentPatterns = new List<string>();
            List<int> Count = new List<int>();
            List<string> Kmers = new List<string>();
            int MaxCount = 0;
            int localCount = 0;

            string pattern = "";

            for (int i = 0; i < Text.Length - k; i++)
            {
                pattern = Text.Substring(i, k);

                localCount = PatternCount(Text, pattern);
                if (localCount > MaxCount)
                    MaxCount = localCount;

                Count.Add(localCount);
                Kmers.Add(pattern);
            }

            List<string> tmpKmers = new List<string>();

            for (k = 0; k < Count.Count; k++)
            {
                if (Count[k] == MaxCount && Count[k] >= t)
                    tmpKmers.Add(Kmers[k]);
            }

            FrequentPatterns = RemoveDuplicates(tmpKmers);

            return FrequentPatterns;

        }

        static public List<string> FrequentWordsWithMismatches(string Text, int k, int d, bool IncludeReverse = false)
        {
            List<string> FrequentPatterns = new List<string>();
            int[] localCount;
            int MaxCount = 0;
            string pattern = "";
            string str = "";

            localCount = ComputingFrequenciesWithMismatches(Text, k, d, IncludeReverse);
            for (int m = 0; m < localCount.Length; m++)
            {
                if (localCount[m] > MaxCount)
                    MaxCount = localCount[m];
            }

            List<string> tmpKmers = new List<string>();

            for (int i = 0; i < localCount.Length; i++)
            {
                if (localCount[i] == MaxCount)
                {
                    pattern = NumberToPattern2(i, k);
                    tmpKmers.Add(pattern);
                }
            }

            FrequentPatterns = RemoveDuplicates(tmpKmers);

            foreach (string s in FrequentPatterns)
            {
                str += s + " ";
            }

            Console.WriteLine("FrequentWordsWithMismatches:  {0}", str);

            return FrequentPatterns;
        }



        static public string ReverseComplement(string input)
        {
            string output = "";

            for (int i = input.Length - 1; i >= 0; i--)
            {
                if (input[i] == 'A')
                    output += "T";
                else if (input[i] == 'T')
                    output += "A";
                else if (input[i] == 'G')
                    output += "C";
                else if (input[i] == 'C')
                    output += "G";
                else
                    Console.WriteLine("ReverseComplement:  Invalid Value");
            }

            //            Console.WriteLine("ReverseComplment:  {0}", output);
            return output;
        }

        static public string PatternMatchIndexes(string pattern, string text)
        {
            string result = "";

            for (int i = 0; i < text.Length - pattern.Length + 1; i++)
            {
                if (text.Substring(i, pattern.Length) == pattern)
                    result += i.ToString() + " ";
            }
            Console.WriteLine("PatternMatchIndexes:  {0}", result);
            return result;
        }

        static public List<string> RemoveDuplicates(List<string> inputList)
        {
            List<string> outputList = new List<string>();

            for (int i = 0; i < inputList.Count; i++)
            {
                if (i == 0)
                    outputList.Add(inputList[i]);
                else
                {
                    if (!outputList.Contains(inputList[i]))
                        outputList.Add(inputList[i]);
                }
            }

            return outputList;
        }

        private static string MinimumSkew(String strInput)
        {
            string strResult = "";
            int count = 0;
            int skewValue = 0;
            int minIndex = 1000;
            string tmp = "";
            int[] arrayOffset = new int[strInput.Length + 1];

            // set first postion to 0;
            arrayOffset[0] = skewValue;
            Console.WriteLine($"Position {count}: Value: skewIndex:  {skewValue}");
            tmp += skewValue.ToString() + " ";

            foreach (char strKmer in strInput)
            {

                if (strKmer == 'G')
                    skewValue += 1;
                else if (strKmer == 'C')
                    skewValue -= 1;

                if (skewValue < minIndex)
                    minIndex = skewValue;

                Console.WriteLine($"Position {count}: Value: {strKmer}   skewIndex:  {skewValue}");

                tmp += skewValue.ToString() + " ";

                count++;

                arrayOffset[count] = skewValue;
            }

            for (int i = 0; i < arrayOffset.Length; i++)
            {
                if (arrayOffset[i] == minIndex)
                    strResult += i.ToString() + " ";
            }
            Console.WriteLine($"Minimum SkewValue:  {minIndex}");

            return strResult;
        }

        // Return index within text that matches pattern with "n" mismatches.
        static public int ApproximatePatternCount(string strPattern, string strText, int nDistance)
        {
            int nResult = 0;
            int nSubDistance = 0;
            string strTmp;
            string strOutput = "";

            for (int i = 0; i < strText.Length - strPattern.Length + 1; i++)
            {
                strTmp = strText.Substring(i, strPattern.Length);
                nSubDistance = HammingDistance(strPattern, strText.Substring(i, strPattern.Length));
                if (nSubDistance <= nDistance)
                {
                    strOutput += strTmp + " ";
                }
            }

            string[] strKmers = strOutput.Split(" ");
            nResult = strKmers.Length - 1;

            return nResult;
        }

        // Return number of kmers within text that match pattern with "n" mismatches
        static public string ApproximatePatternMatching(string strPattern, string strText, int nDistance)
        {
            string strOutput = "";
            int nSubDistance = 0;
            string strTmp;

            for (int i = 0; i < strText.Length - strPattern.Length + 1; i++)
            {
                nSubDistance = HammingDistance(strPattern, strText.Substring(i, strPattern.Length));
                if (nSubDistance <= nDistance)
                {
                    strTmp = i.ToString();
                    strOutput += strTmp + " ";
                }
            }

            return strOutput;
        }

        static public int HammingDistance(string str1, string str2)
        {
            int nDistance = 0;  // assume strings equal;

            if (str1.Length == str2.Length)
            {
                for (int i = 0; i < str1.Length; i++)
                {
                    if (str1[i] != str2[i])
                        nDistance += 1;
                }
            }
            else
            {
                Console.WriteLine("Hamming Distance:  Strings different lengths");
                nDistance = -1;
            }

            return nDistance;
        }

        static public string NumberToPattern2(Int64 index, int k)
        {
            string result = "";
            Int64 prefixIndex = 0;
            int remainder = 0;
            string symbol = "";
            string prefixPattern = "";

            if (k == 1)
                return NumberToSymbol(index);

            prefixIndex = index / 4;
            remainder = (int)index % 4;
            symbol = NumberToSymbol(remainder);
            prefixPattern = NumberToPattern2(prefixIndex, k - 1);

            result = prefixPattern + symbol;

            //            Console.WriteLine("NumberToPattern2:  {0}", result);
            return result;
        }

        static public string NumberToPattern(int index, int kmerSize)
        {
            string result = "";

            List<string> kmerList = new List<string>();

            kmerList = BuildLexicographicallyList(kmerSize);

            result = kmerList[index];

            Console.WriteLine("NumberToPattern:  {0}", result);
            return result;
        }

        static public int SymbolToNumber(string symbol)
        {
            int result = -1;
            if (symbol == "A")
                result = 0;
            else if (symbol == "C")
                result = 1;
            else if (symbol == "G")
                result = 2;
            else if (symbol == "T")
                result = 3;

            return result;
        }

        static public string NumberToSymbol(Int64 index)
        {
            string result = "";

            if (index == 0)
                result = "A";
            else if (index == 1)
                result = "C";
            else if (index == 2)
                result = "G";
            else if (index == 3)
                result = "T";

            return result;
        }

        static public Int64 PatternToNumber2(string pattern)
        {
            Int64 result = 0;
            string symbol = "";
            string prefix = "";

            if (pattern.Length == 0)
                return result;

            symbol = pattern.Substring(pattern.Length - 1, 1);
            prefix = pattern.Substring(0, pattern.Length - 1);

            result = 4 * PatternToNumber2(prefix) + SymbolToNumber(symbol);


            //            Console.WriteLine("PatternToNumber2:  {0}", result);
            return result;
        }

        static public int PatternToNumber(string pattern)
        {
            int result = 0;
            List<string> kmerList = new List<string>();

            kmerList = BuildLexicographicallyList(pattern.Length);

            result = kmerList.BinarySearch(pattern);

            Console.WriteLine("PatternToNumber:  {0}", result);
            return result;
        }

        static public List<string> BuildLexicographicallyList(int kmerSize)
        {
            List<string> outputList = new List<string>();

            for (int i = 0; i < kmerSize; i++)
            {
                outputList = AddTreeLevel(outputList);
            }

            outputList.Sort();

            return outputList;
        }
        static public List<string> AddTreeLevel(List<string> inputList)
        {
            char[] nucleotide = { 'A', 'C', 'T', 'G' };

            List<string> outputList = new List<string>();

            if (inputList.Count > 0)
            {
                foreach (string str in inputList)
                {
                    foreach (char c in nucleotide)
                    {
                        outputList.Add(str + c.ToString());
                    }

                }

            }
            else
            {
                foreach (char c in nucleotide)
                {
                    outputList.Add(c.ToString());
                }
            }

            return outputList;
        }

        static public string[] MyReadFile(string inputFile)
        {
            string[] readText = { };

            try
            {
                readText = File.ReadAllLines(inputFile);

            }
            catch (Exception e)
            {
                Console.WriteLine("{0} Exception caught.", e);
            }

            return readText;
        }

        static public void WriteListToFile(string fileName, List<string> strList)
        {
            if (File.Exists(fileName))
            {
                File.Delete(fileName);
            }

            using (System.IO.StreamWriter myfile = new System.IO.StreamWriter(fileName, true))
            {
                try
                {
                    //                    myfile.WriteLine("Length:  {0}", strList.Count);
                    foreach (string s in strList)

                        myfile.WriteLine(s);

                }
                catch (Exception e)
                {
                    Console.WriteLine("{0} Exception caught.", e);
                }
            }

        }

        static public string ReplaceStringAt(string Text, char c, int index)
        {
            StringBuilder sb = new StringBuilder(Text);
            sb[index] = c;
            string newStr = sb.ToString();
            return newStr;
        }
        private static void Main(string[] args)
        {
            string[] fileText = MyReadFile(inputFile);

            //            string[] dnaArray = new string[fileText.Length - 2];
            string Text = "";
            int k = 0;
            int d = 0;
            int t = 0;
            int N = 1;
            string strResult = "";

            // Set Method to run


            if ("StringSpelledByGenomePath" == method)
            {
                //                string[] strLine = new string[fileText.Length];
                Int32.TryParse(fileText[0], out k);

                //                StringSpelledByGenomePath(fileText);



            }


            if ("Composition" == method)
            {
                //                string[] strLine = new string[fileText.Length];
                Int32.TryParse(fileText[0], out k);
                List<string> searchResult = new List<string>();

                searchResult = Composition(fileText[1], k);

                WriteListToFile("C:\\Temp\\output.txt", searchResult);

            }


            if ("Random" == method)
            {
                double[] array = { .1, .2, .3, .2 };
                int index = Random(array);
            }

            if ("GreedyMotifSearch" == method)
            {
                string[] strLine = new string[fileText.Length - 2];
                Int32.TryParse(fileText[0], out k);
                Int32.TryParse(fileText[1], out t);
                List<string> searchResult = new List<string>();

                for (int i = 2; i < fileText.Length; i++)
                {
                    strLine[i - 2] = fileText[i];
                }

                searchResult = GreedyMotifSearch(strLine, k, t);
            }

            if ("RandominzedMotifSearch" == method)
            {
                string[] strLine = new string[fileText.Length - 2];
                Int32.TryParse(fileText[0], out k);
                Int32.TryParse(fileText[1], out t);
                List<string> searchResult = new List<string>();

                for (int i = 2; i < fileText.Length; i++)
                {
                    strLine[i - 2] = fileText[i];
                }

                searchResult = RandominzedMotifSearch(strLine, k, t);
            }

            if ("GibbsSampler" == method)
            {
                string[] strLine = new string[fileText.Length - 3];
                Int32.TryParse(fileText[0], out k);
                Int32.TryParse(fileText[1], out t);
                List<string> searchResult = new List<string>();
                Int32.TryParse(fileText[2], out N);



                for (int i = 3; i < fileText.Length; i++)
                {
                    strLine[i - 3] = fileText[i];
                }

                searchResult = GibbsSampler(strLine, k, t, N);
            }

            if ("MedianString" == method)
            {
                string[] strLine = new string[fileText.Length - 1];

                Int32.TryParse(fileText[0], out k);
                List<string> searchResult = new List<string>();

                for (int i = 1; i < fileText.Length; i++)
                {
                    strLine[i - 1] = fileText[i];
                }

                MedianString(strLine, k);
            }


            if ("ProfileMostProbableKmer" == method)
            {
                string[] strLine;
                Text = fileText[0];
                Int32.TryParse(fileText[1], out k);
                double[,] array = new double[NucleotideSize, k];
                double value = 0;

                for (int i = 2; i < fileText.Length; i++)
                {
                    strLine = fileText[i].Split(" ");
                    for (int n = 0; n < k; n++)
                    {
                        Double.TryParse(strLine[n], out value);
                        array[i - 2, n] = value;
                    }
                }

                strResult = ProfileMostProbableKmer(Text, k, array);
            }


            if ("Motif" == method)
            {
                string[] strLine;
                string[] strLine2 = new string[5];
                double[,] array = new double[NucleotideSize, 4];
                double value = 0;

                for (int i = 0; i < 4; i++)
                {
                    strLine = fileText[i].Split(" ");
                    for (int n = 0; n < 4; n++)
                    {
                        Double.TryParse(strLine[n], out value);
                        array[i, n] = value;
                    }
                }

                int count = 0;
                for (int i = 4; i < fileText.Length; i++)
                {
                    strLine2[count] = fileText[i];
                    count += 1;
                }


                Motif(array, strLine2, 4);
            }



            //            strResult = MedianString(dnaArray, d);

            //            int dis = DistanceBetweenPatternAndStrings(str2, dnaArray);

            //            strResult = MotifEnumeration(dnaArray, k, d);


            // strResult = ImmediateNeighbors("ATG");

            // strResult = MinimumSkew(str1);

            //  strResult = PatternToNumber2(str1);

            //  ReverseComplement(str1);

            //  dist = HammingDistance(str1, str2);

            //nb = Neighbors(str1, nDistance);

            //                FrequentWords(str1, nDistance);

            //                result = NumberToPattern2(index, nDistance);

            // result = ComputingFrequencies(str1, nDistance);

            //                PatternToNumber("ATGCAA");

            //              NumberToPattern(5437, 8);

            //            PatternMatchIndexes(str1, str2);


            if ("OverlapGraph" == method)
            {
                List<string> output = new List<string>();
                List<string> input = new List<string>();

                foreach (string s in fileText)
                    input.Add(s);

                output = OverlapGraph(input);
            }

            if ("DeBruijnGraph" == method)
            {
                string[] strLine = new string[fileText.Length];
                List<string> searchResult = new List<string>();
                Int32.TryParse(fileText[0], out k);

                searchResult = DeBruijnGraph(k, fileText[1]);

                WriteListToFile("C:\\Temp\\output.txt", searchResult);

            }

            if ("DeBruijnGraph2" == method)
            {
                List<string> strLine = new List<string>();
                List<string> searchResult = new List<string>();

                foreach (string s in fileText)
                    strLine.Add(s);

                searchResult = DeBruijnGraph(strLine);

                WriteListToFile("C:\\Temp\\output.txt", searchResult);

            }
            if ("EulerianCycle" == method)
            {
                List<string> strLine = new List<string>();
                string result = "";

                foreach (string s in fileText)
                    strLine.Add(s);

                result = EulerianCycle(strLine);

                //                WriteListToFile("C:\\Temp\\output.txt", result);

            }
            if ("StringReconstruction" == method)
            {
                List<string> strLine = new List<string>();
                string result = "";
                Int32.TryParse(fileText[0], out k);

                for (int i = 1; i < fileText.Length; i++)
                    strLine.Add(fileText[i]);

                result = StringReconstruction(k, strLine);

                //                WriteListToFile("C:\\Temp\\output.txt", result);

            }


            if ("UniversalCircularString" == method)
            {
                List<string> strLine = new List<string>();
                string result = "";
                Int32.TryParse(fileText[0], out k);

                result = UniversalCircularString(k);
            }

            if ("StringSpelledByGappedPatterns" == method)
            {
                List<string> strLine = new List<string>();
                string result = "";
                Int32.TryParse(fileText[0], out k);
                Int32.TryParse(fileText[1], out d);

                for (int i = 2; i < fileText.Length; i++)
                    strLine.Add(fileText[i]);

                result = StringSpelledByGappedPatterns(strLine, k, d);

                //                WriteListToFile("C:\\Temp\\output.txt", result);

            }

            if ("StringReconstructionFromReadPairs" == method)
            {
                List<string> strLine = new List<string>();
                string result = "";
                Int32.TryParse(fileText[0], out k);
                Int32.TryParse(fileText[1], out d);

                for (int i = 2; i < fileText.Length; i++)
                    strLine.Add(fileText[i]);

                result = StringReconstructionFromReadPairs(k, d, strLine);

                //                WriteListToFile("C:\\Temp\\output.txt", result);

            }

            if ("MaximalNonBranchingPaths" == method)
            {
                List<string> strLine = new List<string>();
                List<string> result = new List<string>();

                for (int i = 0; i < fileText.Length; i++)
                    strLine.Add(fileText[i]);

                result = MaximalNonBranchingPaths(strLine);

                WriteListToFile("C:\\Temp\\output.txt", result);

            }

            if ("ContigGeneration" == method)
            {
                List<string> result = new List<string>();
                List<string> strLine = new List<string>();

                for (int i = 0; i < fileText.Length; i++)
                    strLine.Add(fileText[i]);

                result = ContigGeneration(strLine);

                WriteListToFile("C:\\Temp\\output.txt", result);
            }

            if ("ProteinTranslation" == method)
            {
                string result = "";

                result = ProteinTranslation(fileText[0]);

                //                WriteListToFile("C:\\Temp\\output.txt", result);
            }

            if ("PeptideEncoding" == method)
            {
                List<string> result = new List<string>();

                result = PeptideEncoding(fileText[0], fileText[1]);

                WriteListToFile("C:\\Temp\\output.txt", result);
            }

            if ("LinearSpectrum" == method)
            {
                List<int> result = new List<int>();
                List<string> result2 = new List<string>();
                string output = "";
                int index = 0;

                result = LinearSpectrum(fileText[0]);

                for (int i = 0; i < result.Count; i++)
                {
                    result2.Insert(i, result[i].ToString());
                    output += result[i].ToString() + " ";
                }

                WriteListToFile("C:\\Temp\\output.txt", result2);
            }

            if ("CyclicSpectrum" == method)
            {
                List<int> result = new List<int>();
                List<string> result2 = new List<string>();
                string output = "";
                int index = 0;

                result = CyclicSpectrum(fileText[0]);

                WriteListToFile("C:\\Temp\\output.txt", result2);
            }

            if ("CreateSubpeptides" == method)
            {
                List<int> result = new List<int>();
                List<string> result2 = new List<string>();
                string output = "";
                int index = 0;

                CreateSubpeptides(fileText[0]);

                for (int i = 0; i < result.Count; i++)
                {
                    result2.Insert(i, result[i].ToString());
                    output += result[i].ToString() + " ";
                }

                WriteListToFile("C:\\Temp\\output.txt", result2);
            }

            if ("TheoreticalSpectrum" == method)
            {
                List<int> result = new List<int>();
                List<string> result2 = new List<string>();
                string output = "";
                int index = 0;

                result = TheoreticalSpectrum(fileText[0]);

                for (int i = 0; i < result.Count; i++)
                {
                    result2.Insert(i, result[i].ToString());
                    output += result[i].ToString() + " ";
                }

                WriteListToFile("C:\\Temp\\output.txt", result2);

            }

            if ("BFCyclopeptideSequencing" == method)
            {
                Int64 result = 0;
                Int32.TryParse(fileText[0], out k);

                result = BFCyclopeptideSequencing(k);
            }

            if ("CyclopeptideSequencing" == method)
            {
                string result = "";

                result = CyclopeptideSequencing(fileText[0]);
            }

            if ("LinearScore" == method)
            {
                int result = 0;

                result = LinearScore(fileText[0], fileText[1], false);
            }

            if ("LeaderboardCyclopeptideSequencing" == method)
            {
                string result = "";
                Int32.TryParse(fileText[1], out d);

                result = LeaderboardCyclopeptideSequencing(fileText[1], d);
            }

            if ("ConvolutionCyclopeptideSequencing" == method)
            {
                string result = "";
                Int32.TryParse(fileText[1], out d);
                int M = 0;
                Int32.TryParse(fileText[0], out M);

                result = LeaderboardCyclopeptideSequencing(fileText[2], d, M);
            }

            if ("TrimLeaderboard" == method)
            {
                Int32.TryParse(fileText[2], out d);
                List<string> result = new List<string>();
                List<string> leaderboard = new List<string>();
                string[] array = fileText[0].Split(" ");
                foreach (string s in array)
                    leaderboard.Add(s);

                result = TrimLeaderboard(leaderboard, fileText[1], d);

                WriteListToFile("C:\\Temp\\output.txt", result);

            }

            if ("SpectralConvolution" == method)
            {
                string result = "";

                result = SpectralConvolution(fileText[0]);

//                WriteListToFile("C:\\Temp\\output.txt", result);

            }

            if ("DPChange" == method)
            {
                int result = 0;
                int money = 0;
                int val = 0;
                List<int> Coins = new List<int>();
                string[] tmpCoins = fileText[1].Split(",");

                foreach (string s in tmpCoins)
                {
                    Int32.TryParse(s, out val);
                    Coins.Add(val);

                }

                Int32.TryParse(fileText[0], out money);

                result = DPChange(money, Coins);
            }


            if ("ManhattanTourist" == method)
            {
                int result = 0;
                string[] array = fileText[0].Split(" ");
                int n = 0;
                int m = 0;
                Int32.TryParse(array[0], out n);
                Int32.TryParse(array[0], out m);
                int[,] rightArray = new int[n + 1, m];
                int[,] downArray = new int[n, m + 1];

                bool fillRight = false;
                int row = 0;

                for (int i = 1; i < fileText.Length; i++)
                {
                    if (fileText[i] == "-")
                    {
                        fillRight = true;
                        row = 0;
                    }
                    else
                    {
                        int val = 0;
                        string[] intValues = fileText[i].Split(" ");

                        for (int s = 0; s < intValues.Length; s++)
                        {
                            Int32.TryParse(intValues[s], out val);
                            if (fillRight)
                                rightArray[row, s] = val;
                            else
                                downArray[row, s] = val;
                        }
                        row++;
                    }
                }

                result = ManhattanTourist(n, m, downArray, rightArray);
            }

            if ("LCSBackTrack" == method)
            {
                string[,] result;

                result = LCSBackTrack(fileText[0], fileText[1]);

                int test = 0;
            }

            if ("OutputLCS" == method)
            {
                string result;
                string v = fileText[0];
                string w = fileText[1];

                string[,] backtrack = LCSBackTrack(v, w);

                result = OutputLCS(backtrack, v, v.Length, w.Length);

                int test = 0;
            }

            if ("MassCounting" == method)
            {
                Int64 result;
                string v = fileText[0];
                string w = fileText[1];

                string[,] backtrack = LCSBackTrack(v, w);
                Int64 weight = 22;
                List<Int64> masses = new List<Int64>();
                masses.Add(2);
                masses.Add(3);

                result = MassCounting(weight, masses);

                int test = 0;
            }

            /*            if ("TopologicalSort" == method)
                        {
                            //        List<T> TopologicalSort<T>(HashSet<T> nodes, HashSet<Tuple<T, T>> edges) where T : IEquatable<T>

                            HashSet<int> nodes = new HashSet<int>{ 7, 5, 3, 8, 11, 2, 9, 10 };
                            HashSet<Tuple<int, int>> edges = new HashSet<Tuple<int, int>> { Tuple.Create(7, 11),
                                    Tuple.Create(7, 8),
                                    Tuple.Create(5, 11),
                                    Tuple.Create(3, 8),
                                    Tuple.Create(3, 10),
                                    Tuple.Create(11, 2),
                                    Tuple.Create(11, 9),
                                    Tuple.Create(11, 10),
                                    Tuple.Create(8, 9)};


            //                var ret = TopologicalSort(nodes, edges);

            //                System.Diagnostics.Debug.Assert(ret.SequenceEqual(new[] { 7, 5, 11, 2, 3, 8, 9, 10 }));
                        }   */

            if ("TopologicalSort" == method)
            {
                List<string> result = new List<string>();
                int startingNode = 0;
                int endingNode = 0;
                List<string> edges = new List<string>();

                Int32.TryParse(fileText[0], out startingNode);
                Int32.TryParse(fileText[1], out endingNode);

                for (int i = 2; i < fileText.Length; i++)
                {
                    edges.Add(fileText[i]);
                }

                result = LongestPathInDAG(startingNode, endingNode, edges);

                int test = 0;
            }

            if ("LongestPathInDAG" == method)
            {
                List<string> result = new List<string>();
                int startingNode = 0;
                List<string> result2 = new List<string>();

                int endingNode = 0;
                List<string> edges = new List<string>();

                for (int i = 2; i < fileText.Length; i++)
                {
                    edges.Add(fileText[i]);
                }

                Int32.TryParse(fileText[0], out startingNode);
                Int32.TryParse(fileText[1], out endingNode);

                result = LongestPathInDAG(startingNode, endingNode, edges);

                int test = 0;
            }

            if ("SequenceAlignment" == method)
            {
                SequenceAlignment seq = new SequenceAlignment();

                seq.SequenceAlign(fileText[0], fileText[1]);
            }

            if ("EditDistance" == method)
            {
                SequenceAlignment seq = new SequenceAlignment();
                int distance = 0;

                distance = seq.EditDistance(fileText[0], fileText[1]);
            }

            if ("FittingAlignment" == method)
            {
                SequenceAlignment seq = new SequenceAlignment();
                List<string> result = new List<string>();

                result = seq.FittingAlignment(fileText[0], fileText[1]);
            }

            if ("OverlapAlignment" == method)
            {
                SequenceAlignment seq = new SequenceAlignment();
                List<string> result = new List<string>();

                result = seq.OverlapAlignment(fileText[0], fileText[1]);
            }

            if ("AlignmentWithAffineGapPenalties" == method)
            {
                SequenceAlignment seq = new SequenceAlignment();
                List<string> result = new List<string>();

                result = seq.AlignmentWithAffineGapPenalties(fileText[0], fileText[1]);
            }

            if ("FindMiddleEdge" == method)
            {
                SequenceAlignment seq = new SequenceAlignment();
                string result = "";

                result = seq.FindMiddleEdge(fileText[0], fileText[1]);
            }

            if ("multipleSequenceAlignment" == method)
            {
                SequenceAlignment seq = new SequenceAlignment();
                string[] result;

                result = seq.multipleSequenceAlignment(fileText[0], fileText[1], fileText[2]);
            }

            if ("GreedySorting" == method)
            {
                List<string> result = new List<string>();

                result = GreedySorting(fileText[0]);

                WriteListToFile("C:\\Temp\\output.txt", result);

            }

            if ("NumberOfBreakpoints" == method)
            {
                int result = 0;

                result = NumberOfBreakpoints(fileText[0]);
            }

            if ("ChromosomeToCycle" == method)
            {
                string result = "";

                result = ChromosomeToCycle(fileText[0]);
            }

            if ("CycleToChromosome" == method)
            {
                string result = "";

                result = CycleToChromosome(fileText[0]);
            }

            if ("ColoredEdges" == method)
            {
                string result = "";

                result = ColoredEdges(fileText[0]);
            }

            if ("GraphToGenome" == method)
            {
                string result = "";

                result = GraphToGenome(fileText[0]);
            }

            if ("TwoBreakDistance" == method)
            {
                int result = 0;

                result = TwoBreakDistance(fileText[0], fileText[1]);
            }

            if ("TwoBreakOnGenomeGraph" == method)
            {
                string result = "";
                string[] myInput = fileText[1].Split(",");
                int i1, i2, i3, i4;

                i1 = Int32.Parse(myInput[0]);
                i2 = Int32.Parse(myInput[1]);
                i3 = Int32.Parse(myInput[2]);
                i4 = Int32.Parse(myInput[3]);

                result = TwoBreakOnGenomeGraph(fileText[0], i1, i2, i3, i4);
            }

            if ("TwoBreakOnGenome" == method)
            {
                string result = "";
                string[] myInput = fileText[1].Split(",");
                int i1, i2, i3, i4;

                i1 = Int32.Parse(myInput[0]);
                i2 = Int32.Parse(myInput[1]);
                i3 = Int32.Parse(myInput[2]);
                i4 = Int32.Parse(myInput[3]);

                result = TwoBreakOnGenome(fileText[0], i1, i2, i3, i4);
            }

            if ("sharedKmers" == method)
            {
                List<string> result = new List<string>();
                k = Int32.Parse(fileText[0]);

                result = sharedKmers(k, fileText[1], fileText[2]);

                WriteListToFile("C:\\Temp\\output.txt", result);
            }

        }


        private static int arrayComparator(int[] o1, int[] o2)
        {
            if (o1[0] < o2[0])
            {
                return -1;
            }
            else if (o1[0] > o2[0])
            {
                return 1;
            }
            else
            {
                if (o1[1] < o2[1])
                {
                    return 1;
                }
                else if (o1[1] > o2[1])
                {
                    return -1;
                }
                else
                    return 0;
            }
         }

    public static int Max(params int[] numbers)
        {
            return numbers.Max();
        }

        public class MySortingClass : IComparer<Peptide>
        {
            public int Compare(Peptide x, Peptide y)
            {
                int compareScore = x.score.CompareTo(y.score);
                if (compareScore == 0)
                {
                    return x.code.CompareTo(y.code);
                }
                return compareScore;
            }
        }
        public class Peptide
        {
            public string code { get; set; }
            public int mass { get; set; }
            public int score { get; set; }

            public bool cyclic { get; set; }

            public int count { get; set; }

            Peptide()
            {
                this.code = "";
                this.mass = 0;
                this.score = 0;
                this.cyclic = false;
                this.count = 0;
            }

            public Peptide(string code, int score)
            {
                this.code = code;
                this.score = score;
                this.mass = 0;
                this.cyclic = false;
                this.count = 0;
            }
        }

        public class Node
        {
            public string node { get; set; }
            public List<string> connectedNodes { get; set; }
            public int inDegrees { get; set; }
            public int outDegrees { get; set; }
            public bool included { get; set; }

            public Node()
            {
                this.node = node;

                connectedNodes = new List<string>();
                inDegrees = outDegrees = 0;
                included = false;
            }

            public Node(string node)
            {
                this.node = node;
                connectedNodes = new List<string>();
                inDegrees = outDegrees = 0;
                included = false;
            }
        }

        public class Circuit
        {
            public List<string> circuit;
            public int currentPos;

            public Circuit()
            {
                circuit = new List<string>();
                currentPos = 0;
            }

            public void AddNode(Node newNode)
            {
                circuit.Add(newNode.node);
                currentPos++;
            }
        }

        public static void PL(object o) { Console.WriteLine(o); } //alias
        public static void PL() { Console.WriteLine(); } //alias
        public static void P(object o) { Console.Write(o); } //alias

    }
}