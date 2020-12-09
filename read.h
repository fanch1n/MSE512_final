#include <iostream>
#include <iterator>
#include <fstream>
#include <vector>
#include <algorithm> // for std::copy
bool getFileContent(std::string fileName, std::vector<std::string> & vecOfStrs);
std::vector<double> getData(std::string fileName);

// int main(){
//       // Get the contents of file in a vector
//     std::vector<std::string> vecOfStr;
//     std::cout<< "read bcc.txt" << std::endl;
//     bool result = getFileContent("bcc.txt", vecOfStr);
//     std::vector<double> vecData;
//     vecData = getData("bcc.txt");
//     for(int i=0; i< vecData.size(); i++){
//         std::cout << vecData[i] << std::endl;
//     }
//
  //
  //  if(result)
  //  {
  //        // Print the vector contents
  //      for(std::string & line : vecOfStr)
  //      {
  //          std::cout<<line<<std::endl;
  //      }
  //   }
  //  return 0;
// }
/*----------------------------------------*/
bool getFileContent(std::string fileName, std::vector<std::string> & vecOfStrs)
{
    // Open the File
    std::ifstream in(fileName.c_str());
    // Check if object is valid
    if(!in)
    {
        std::cerr << "Cannot open the File : "<<fileName<<std::endl;
        return false;
    }
    std::string str;
    // Read the next line from File untill it reaches the end.
    while (std::getline(in, str))
    {
        // Line contains string of length > 0 then save it in vector
        if(str.size() > 0)
            vecOfStrs.push_back(str);
    }
    //Close The File
    in.close();
    return true;
}

/*----------------------------------------*/
std::vector<double> getData(std::string fileName)
{
    std::vector<double> data;
    // Open the File
    std::ifstream infile(fileName, std::ios::in);
    // Check if object is valid
    if(!infile.is_open())
    {
        std::cerr << "Cannot open the File : "<<fileName<<std::endl;
        exit(1);
    }

    double num = 0.0;
    //keep storing values from the text file as long as the data exists;
    while (infile >> num )
    {
        data.push_back(num);
    }

    //Close The File
    infile.close();
    return data;
}
