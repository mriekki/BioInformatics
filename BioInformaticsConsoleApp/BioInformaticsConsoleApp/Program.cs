using System;
using System.IO;
using System.Collections.Generic;


namespace BioInformaticsConsoleApp
{
    class Program
    {
        static void Main(string[] args)
        {
            string inputFile = "c:\\Temp\\large.txt";

            string[] fileText = MyReadFile(inputFile);

            // Minimum Skew
            if (fileText.Length == 1)
            {
                ReverseComplement(fileText[0]);

                Console.WriteLine(MinimumSkew(fileText[0]));
            }

            // Hamming Distance Test
            /*
            if (fileText.Length == 2)
            {
                string str1 = fileText[0];
                string str2 = fileText[1];
                int nDistance = HammingDistance(str1, str2);

                Console.WriteLine("Hamming Distance equals {0}", nDistance);

            }       */


            if (fileText.Length == 2)
            {
                string str1 = fileText[0];
                string str2 = fileText[1];
                int nDistance = 0;

                //              Int32.TryParse(str2, out nDistance);


                //                FrequentWords(str1, nDistance);
                PatternToNumber("ATGCAA");

                NumberToPattern(5437, 8);

                PatternMatchIndexes(str1, str2);
            }

            if (fileText.Length == 2)
            {
                string str1 = fileText[0];
                string str2 = fileText[1];
                int nDistance = 0;
                string strNeighborhood = "";

                Int32.TryParse(str2, out nDistance);
//                strNeighborhood = ImmediateNeighbors(str1);

                strNeighborhood = Neighbors(str1, nDistance);
            }


            if (fileText.Length == 4)
            {
                string str1 = fileText[0];
                string str2 = fileText[1];
                string str3 = fileText[2];
                string str4 = fileText[3];
                string result = "";

                int nDistance = 0;
                int k = 0;
                int L = 0;
                int t = 0;

                Int32.TryParse(str2, out k);
                Int32.TryParse(str3, out L);
                Int32.TryParse(str4, out t);

                result = ClumpFinder(str1, k, L, t);

                k = 1;
            }

        }

        static public string ClumpFinder(string Genome, int k, int L, int t)
        {
            string result = "";
            string pattern = "";
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
            int[] arrayOffset = new int[strInput.Length + 1];

            foreach (char strKmer in strInput)
            {
                if (strKmer.ToString().Length > 0)
                {
                    if (count > 0)
                    {
                        if (strKmer == 'G')
                            skewValue++;
                        else if (strKmer == 'C')
                            skewValue--;
                    }
                    if (skewValue < minIndex)
                        minIndex = skewValue;

                    Console.WriteLine($"Position {count}: Value: {strKmer}   skewIndex:  {skewValue}");

                    count++;

                    arrayOffset[count] = skewValue;

                }
            }

            for (int i = 0; i < arrayOffset.Length; i++)
            {
                if (arrayOffset[i] == minIndex)
                    strResult += i.ToString() + " ";
            }
            Console.WriteLine($"Minimum SkewValue:  {minIndex}");

            return strResult;
        }

        static public string ComputingFrequenciesWithMismatches(string strInputText, int nKmer, int nDistance)
        {
            string strResult = "";
            int[] freqArray = new int[(int)Math.Pow(4, nKmer)];
            string strPattern = "";

            for (int i = 0; i < freqArray.Length; i++)
                freqArray[i] = 0;

            for (int i = 0; i < strInputText.Length - nKmer; i++)
            {
                strPattern = strInputText.Substring(i, nKmer);
            }

            return strResult;
        }

        static public string Neighbors(string strPattern, int d)
        {
            string strResult = strPattern;
            string strSuffixNeighbors = "";
            char[] nucleotied = { 'A', 'C', 'T', 'G' };

            if (d == 0)
                return strResult;

            if (strPattern.Length == 1)
            {
                strResult = "A, C, G, T";
                
            }

            strResult = "";

            strSuffixNeighbors = Neighbors(strPattern.Substring(1), d);

            foreach (char s in strSuffixNeighbors)
            {
                if (HammingDistance(strPattern.Substring(1), s.ToString()) < d)
                {
                    foreach (char c in nucleotied)
                    {
//                        strResult += ", " + c.ToString + 
                    }

                }
            }

   
            return strResult;
        }

        static public string IterativeNeighbors(string strPattern, int d)
        {
            string strNeighborhood = strPattern;

            for (int i = 1; i <= d; i++)
            {

            }
            return strNeighborhood;
        }

        static public string ImmediateNeighbors(string strPattern)
        {
            string strNeighborhood = strPattern;
            char symbol;
            char[] nucleotied = { 'A', 'C', 'T', 'G' };
            string tmpStr = "";

            for (int i = 0; i < strPattern.Length; i++)
            {
                symbol = strPattern[i];
                foreach (char c in nucleotied)
                {
                    if (symbol != c)
                    {
                        tmpStr = ReplaceString(strPattern, i, c.ToString());
                        
                        strNeighborhood += ", " + tmpStr;
                    }
                }
            }

            Console.WriteLine("ImmediateNeighbors - Neighborhood = {0}", strNeighborhood);

            return strNeighborhood;
        }

        static public string ReplaceString(string strPattern, int index, string strNew)
        {
            string strResult = "";

            for (int i = 0; i < strPattern.Length; i+= strNew.Length)
            {
                if (i == index)
                    strResult += strNew;
                else
                    strResult += strPattern[i];
            }
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

        static public string NumberToPattern(int index, int kmerSize)
        {
            string result = "";

            List<string> kmerList = new List<string>();

            kmerList = BuildLexicographicallyList(kmerSize);

            result = kmerList[index];

            Console.WriteLine("NumberToPattern:  {0}", result);
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

    }

}
