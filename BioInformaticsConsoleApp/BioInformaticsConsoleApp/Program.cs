﻿using System;
using System.IO;

namespace BioInformaticsConsoleApp
{
    class Program
    {
        static void Main(string[] args)
        {
            string inputFile = "D:\\Temp\\dataset_9_4.txt";

            string[] fileText = MyReadFile(inputFile);

            // Hamming Distance Test
            /*
            if (fileText.Length == 2)
            {
                string str1 = fileText[0];
                string str2 = fileText[1];
                int nDistance = HammingDistance(str1, str2);

                Console.WriteLine("Hamming Distance equals {0}", nDistance);

            }       */


            if (fileText.Length == 3)
            {
                string str1 = fileText[0];
                string str2 = fileText[1];
                string str3 = fileText[2];
                int nDistance = 0;

                Int32.TryParse(str3, out nDistance);

                string strIndex = ApproximatePatternMatching(str1, str2, nDistance);

                Console.WriteLine("ApproximatePatternMatching {0}", strIndex);

            }


/*            foreach (string s in fileText)
            {
                Console.WriteLine(s);
            }       */
          }

         private static string MinimumSkew(String strInput)
        {
            string strResult;
            int count = 0;
            int skewValue = 0;
            int minIndex = 1000;

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
                }
            }

            Console.WriteLine($"Minimum SkewValue:  {minIndex}");

            strResult = "hello world";

            return strResult;
        }

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
