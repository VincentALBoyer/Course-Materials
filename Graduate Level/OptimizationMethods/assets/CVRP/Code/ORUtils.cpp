#include "ORUtils.hpp"

double ORUtils::coordtogeo(double val) {
    double pi = 3.141592;
    double deg = (int)val;
    double min = val - deg;
    return pi * (deg + 5.0 * min / 3.0) / 180;
}

// Trim from the start (left)

std::string ORUtils::ltrim(std::string s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
        }));
    return s;
}

// Trim from the end (right)

std::string ORUtils::rtrim(std::string s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
        }).base(), s.end());
    return s;
}

std::string ORUtils::getData(ifstream& input, string keyword) {
    input.clear();
    input.seekg(0, ios::beg);
    string line;
    while (getline(input, line)) {
        size_t pos = line.find(keyword);
        if (pos != string::npos) {
            pos = line.find(':', pos);
            if (pos == string::npos) {
                return "n/a";
            }
            string s = line.substr(pos + 1);
            return trim(s);
        }
    }
    return "n/a";
}

bool ORUtils::gotoSection(ifstream& input, string keyword) {
    input.clear();
    input.seekg(0, ios::beg);
    string line;
    while (getline(input, line)) {
        if (line.find(keyword) != string::npos) {
            return true;
        }
    }
    return false;
}

bool ORUtils::gotoKeyword(ifstream& input, string keyword) {
    input.clear();
    input.seekg(0, ios::beg);
    string buffer = "";
    while (buffer != keyword) {
        input >> buffer;
        if (input.eof()) {
            return false;
        }
    }
    return true;

}

std::string ORUtils::extractfilename(std::string filepath, bool withext) {
    size_t pos = filepath.find_last_of("/\\");
    size_t posext = filepath.find_last_of(".");
    if (withext)
        return filepath.substr(pos + 1);
    else
        return filepath.substr(pos + 1, posext - pos - 1);
}

std::string ORUtils::extractpath(std::string filepath) {
    size_t pos = filepath.find_last_of("/\\");
    return (pos == string::npos) ? "" : filepath.substr(0, pos);
}

double ORUtils::dotproduct(const std::vector<double>& a, const std::vector<double>& b)
{
    if (a.size() != b.size())
        throw std::invalid_argument("Vectors must be of the same size for dot product.");

    // Kahan summation for improved numerical stability
    double sum = 0.0;
    double c = 0.0; // Compensation for lost low-order bits
    for (size_t i = 0; i < a.size(); ++i) {
        double prod = a[i] * b[i];
        double y = prod - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    return sum;
}



// Initialize the static member variable
ORRandom* ORRandom::_instance = nullptr;
