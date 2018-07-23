using System;
using System.IO;
using System.Collections.Generic;
using System.Text;


namespace BioInformaticsConsoleApp
{
    class Program
    {
        private const int NucleotideSize = 4;

        private static void Main(string[] args)
        {
            string inputFile = "..\\..\\..\\Data Files\\dataset_159_3.txt";

            string[] fileText = MyReadFile(inputFile);
            string[] dnaArray = new string[fileText.Length - 2];
            string Text = "";
            string str1 = "";
            string str2 = "";
            string str3 = "";
            string str4 = "";
            int nResult = 0;
            int k = 0;
            int d = 0;
            int t = 0;
            string strResult = "";

            // Set Method to run
            string method = "ProfileMostProbableKmer";

            if ("ProfileMostProbableKmer" == method)
            {
                string [] strLine;
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



            if (fileText.Length == 1)
            {
                str1 = fileText[0];
            }
            else if (fileText.Length == 2)
            {
                str1 = fileText[0];
                str2 = fileText[1];
                Int32.TryParse(str2, out k);
            }
            else if (fileText.Length == 3)
            {
                str1 = fileText[0];
                str2 = fileText[1];
                str3 = fileText[2];
                Int32.TryParse(str2, out k);
                Int32.TryParse(str3, out d);
            }
            else if (fileText.Length == 4)
            {
                str1 = fileText[0];
                str2 = fileText[1];
                str3 = fileText[2];
                str4 = fileText[3];
                Int32.TryParse(str2, out k);
                Int32.TryParse(str3, out d);
                Int32.TryParse(str4, out t);
            }
            else if (fileText.Length > 4)
            {
                str1 = fileText[0];
                str2 = fileText[1];
                Int32.TryParse(str1, out k);
                Int32.TryParse(str2, out d);

                for (int i = 0; i < fileText.Length - 2; i++) // offset of k, d
                    dnaArray[i] = fileText[i + 2];
           
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

        }


        static public string ProfileMostProbableKmer(string Text, int k, double[,] matrix)
        {
            string result = "";
            string pattern = "";
            string c = "";
            double probability = 1;
            double highestProbabilty = 0;

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

                if (probability > highestProbabilty)
                {
                    highestProbabilty = probability;
                    result = pattern;
                }
            }


            Console.WriteLine("ProfileMostProbableKmer:  {0}", result);
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
            using (System.IO.StreamWriter myfile = new System.IO.StreamWriter(fileName, true))
            {
                try
                {
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

    }

}