﻿using System;
using System.IO;
using System.Collections.Generic;


namespace BioInformaticsConsoleApp
{
    class Program
    {
        static void Main(string[] args)
        {
            string inputFile = "D:\\Temp\\dataset_9_6.txt";

            string[] fileText = MyReadFile(inputFile);

            // Minimum Skew
            if (fileText.Length == 1)
            {
                string str1 = fileText[0];
                Int64 result = 0;

                result = PatternToNumber2(str1);
                int k = 0;
//                ReverseComplement(fileText[0]);

//                Console.WriteLine(MinimumSkew(fileText[0]));
            }


            if (fileText.Length == 2)
            {
                string str1 = fileText[0];
                string str2 = fileText[1];
                int nDistance = 0;
                int index = 0;
                string result = "";

                Int32.TryParse(str1, out index);
                Int32.TryParse(str2, out nDistance);


                //                FrequentWords(str1, nDistance);

//                result = NumberToPattern2(index, nDistance);

                  result = ComputingFrequencies(str1, nDistance);

                //                PatternToNumber("ATGCAA");

                //              NumberToPattern(5437, 8);

                //            PatternMatchIndexes(str1, str2);

                int k = 11;
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

            if (fileText.Length == 3)
            {
                string str1 = fileText[0];
                string str2 = fileText[1];
                string str3 = fileText[2];
                int nDistance = 0;
                int nKmer = 0;

                Int32.TryParse(str3, out nDistance);
                Int32.TryParse(str2, out nKmer);

                int n = ApproximatePatternCount(str1, str2, nDistance);

                Console.WriteLine("ApproximatePatternCount {0}", n);

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

            Console.WriteLine("NumberToPattern2:  {0}", result);
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

    }

}
