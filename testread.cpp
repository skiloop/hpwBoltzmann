#include<iostream>
#include<cstdlib>
#include<fstream>
#include<string>
#include<exception>
#include<iomanip>
using namespace std;

int main(int argc, char* argv[]) {

    double f;
    double E0;
    double phi;
    int NumOfWaveLength;
    int NumOfCellPerWaveLen;
    int m;
    int nbound;
    int scatwidth;
    int PlotStep;
    int CStep;

    string fname = "hpw.param";
    string line;
    string freq = "frequency = ";
    ifstream infile("hpw.param");
    int width = 40;
    int vwidth = 22;
    if (!infile.is_open()) {
        cerr << "file open failed!\t" << fname << endl;
        exit(-1);
    }
    try {

        // read frequency        
        infile >> f;
        cout << setiosflags(ios_base::right) << setw(width) << "incident wave frequency:" << resetiosflags(ios_base::left) << setw(vwidth) << f << endl;
        // read amplitude of electric component
        infile >> E0;
        cout << setiosflags(ios_base::right) << setw(width) << "electric component amplitude:" << resetiosflags(ios_base::left) << setw(vwidth) << E0 << endl;
        // read phrase angle
        infile >> phi;
        cout << setiosflags(ios_base::right) << setw(width) << "incident phrase angle:" << resetiosflags(ios_base::left) << setw(vwidth) << phi << endl;
        infile >> NumOfWaveLength;
        cout << setiosflags(ios_base::right) << setw(width) << "Number of wavelength in domain:" << resetiosflags(ios_base::left) << setw(vwidth) << NumOfWaveLength << endl;
        infile >> NumOfCellPerWaveLen;
        cout << setiosflags(ios_base::right) << setw(width) << "Number of coarse cells per wavelength:" << resetiosflags(ios_base::left) << setw(vwidth) << NumOfCellPerWaveLen << endl;
        infile >> m; //number of fine cells per coarse cell
        cout << setiosflags(ios_base::right) << setw(width) << "Number of fine cells per coarse cell:" << resetiosflags(ios_base::left) << setw(vwidth) << m << endl;
        infile >> nbound; //number of coarse cells in PML boundary
        cout << setiosflags(ios_base::right) << setw(width) << "PML width:" << resetiosflags(ios_base::left) << setw(vwidth) << nbound << endl;
        infile >> scatwidth; // number of coarse cells between connecting interface and PML boundary 
        cout << setiosflags(ios_base::right) << setw(width) << "scatter field width:" << resetiosflags(ios_base::left) << setw(vwidth) << scatwidth << endl;
        infile >> PlotStep;
        cout << setiosflags(ios_base::right) << setw(width) << "PlotStep:" << resetiosflags(ios_base::left) << setw(vwidth) << PlotStep << endl;
        infile >> CStep;
        cout << setiosflags(ios_base::right) << setw(width) << "Capture Step:" << resetiosflags(ios_base::left) << setw(vwidth) << CStep << endl;
        infile.close();
    } catch (exception &e) {
        cerr << e.what() << endl;
        exit(-1);
    }

}
