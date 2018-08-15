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
        private static char[] zero_ones = { '0', '1' };

        //private const string inputFile = "..\\..\\..\\Data Files\\dataset_6207_2.txt";
        private const string inputFile = "..\\..\\..\\Data Files\\MyData.txt";
        private const string method = "MaximalNonBranchingPaths";


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
            foreach(Node nd in dict.Values)
            {
                nd.outDegrees += nd.connectedNodes.Count;

                foreach(string sn in nd.connectedNodes)
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
                    sequence += kmers[i].Substring(kmers[i].Length-1, 1);       // add last character of current kmer
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
            int index = array.Length-1;     // Set to last index by default
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
                    distributionArray[i] = ((array[i] / total) + distributionArray[i-1]);
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
        static public int DistanceBetweenPatternAndStrings(string Pattern, string [] Dna)
        {
            int distance = 0;
            int k = Pattern.Length;
            string str = "";
            int tmpDistance = 0;

            foreach (string Text in Dna)
            {
                int HamDistance = int.MaxValue;

                for (int m = 0; m < Text.Length - Pattern.Length +1; m++)
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
                for (int m = 0; m < kmer.Length - k +1; m++)        // All kmers in DNA
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
            int [] localCount;
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

            Console.WriteLine("ReverseComplment:  {0}", output);
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
            //            circuit[currentPos] = newNode.node;
            circuit.Add(newNode.node);
            currentPos++;
        }
    }


}