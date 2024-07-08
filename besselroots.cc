#include "iitpave2.h"
MatrixXd besselroots(int o, int N, int k) {
    std::string filename;
    if (o == 0 && k == 1) {
        filename = "B0r.csv";
    } else if (o == 1 && k == 1) {
        filename = "B1r.csv";  
        cout << "Successfully read Bessel roots." << endl;
    } else {
        std::cerr << "Unsupported combination of parameters." << std::endl;
        return MatrixXd(0, 0); // Return an empty matrix on error
    }

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Unable to open file: " << filename << std::endl;
        return MatrixXd(0, 0); // Return an empty matrix on error
    }

    std::vector<double> roots;
    std::string line;
    while (getline(file, line)) {
        std::stringstream ss(line);
        std::string value;
        while (getline(ss, value, ',')) {
            roots.push_back(std::stod(value));
        }
    }
    file.close();

    MatrixXd x(1, N);
    for (int i = 0; i < N; ++i) {
        x(0, i) = roots[i];
    }

    return x;
}