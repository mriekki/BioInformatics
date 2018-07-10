using System;
using System.IO;

namespace BioInformaticsConsoleApp
{
    class Program
    {
        static void Main(string[] args)
        {
            string inputFile = "D:\\Temp\\dataset_9_3.txt";

            string[] fileText = MyReadFile(inputFile);


            if (fileText.Length == 2)
            {
                string str1 = fileText[0];
                string str2 = fileText[1];
                int nDistance = HammingDistance(str1, str2);

                Console.WriteLine("Hamming Distance equals {0}", nDistance);

            }

/*            foreach (string s in fileText)
            {
                Console.WriteLine(s);
            }       */
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
