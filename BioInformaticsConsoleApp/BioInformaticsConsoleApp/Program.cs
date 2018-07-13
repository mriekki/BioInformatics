using System;
using System.IO;

namespace BioInformaticsConsoleApp
{
    class Program
    {
        static void Main(string[] args)
        {
            string inputFile = "D:\\Temp\\Neighbors.txt";

            string[] fileText = MyReadFile(inputFile);

            // Minimum Skew
            if (fileText.Length == 1)
            {
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

                ComputingFrequenciesWithMismatches(str1, nKmer, nDistance);

//                int n = ApproximatePatternCount(str1, str2, nDistance);

//                Console.WriteLine("ApproximatePatternCount {0}", n);

            }


/*            foreach (string s in fileText)
            {
                Console.WriteLine(s);
            }       */
          }

         private static string MinimumSkew(String strInput)
        {
            string strResult = "";
            int count = 0;
            int skewValue = 0;
            int minIndex = 1000;
            int[] arrayOffset = new int[strInput.Length + 1];

            foreach(char strKmer in strInput)
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
